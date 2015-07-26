'''Reads in output files produced by fiber-optic distributed temperature sensing machines and plots them'''
__author__ = "Aaron Curtis <aarongc@nmt.edu>"
__version__ = '0.1'

# -*- coding: utf-8 -*-
import os, datetime
import pandas
import numpy
import glob2
from matplotlib import pyplot
from scipy.io import loadmat

pyplot.ion() #set up interactive plotting

def read_ddf_col(datadir, ddf_column=0, channel=1):
    '''
    Reads a column (temp, Stokes, AntiStokes) from all the .ddf files in all datadirs
    '''
    for datadir in glob2.glob(datadir + '**/'):
        for filepath in os.listdir(datadir):
            if filepath.endswith('.ddf'):
                timestamp = datetime.datetime.strptime(filepath, 'channel ' + str(channel) + ' %Y%m%d %H%M%S 00001.ddf')
                data_for_timestamp = pandas.read_csv(
                    datadir + filepath,
                    skiprows = 25,
                    delimiter ='\t')
                data_for_timestamp.set_index('length (m)', inplace=True)
                try:
                    data[timestamp] = data_for_timestamp.ix[:,ddf_column] # would be better to do data_for_timestamp['temperature (<deg>C)'] but the degree symbol is a problem
                except NameError:
                    # first iteration
                    data = pandas.DataFrame(data_for_timestamp.ix[:,ddf_column],columns=[timestamp])
                    ref_temps = pandas.DataFrame
    data.sort(axis=1, inplace=True)
    return data

def read_trefs(datadir):
    '''
    Reads all .tcd files in datadirs and returns a DataFrame containing all thermistor data
    '''
    trefs = []
    for tcd_filepath in glob2.glob(datadir + '**/*.tcd'):
        #index_col=False is a workaround to avoid getting
        #*** NotImplementedError: file structure not yet supported
        data = pandas.read_csv(
            tcd_filepath,
            skiprows = 6,
            parse_dates = [[0,1]],
            delimiter ='\t',
            index_col = False)
        data = data.set_index(data.columns[0])
        data.columns=['internal','ext1','ext2'] #necessary because the header row starts with two tab characters and pandas deals with that differently depending on version :-(
        trefs.append(data)
    return pandas.concat(trefs)

def read_ctemps_mat(filepath):
    m = loadmat(filepath)
    times = [datetime.datetime.fromordinal(int(dt))
             + datetime.timedelta(days=dt%1)
             - datetime.timedelta(days=366)
             - datetime.timedelta(hours=5)
             for dt in m['datetime'][0]]
    return pandas.DataFrame(m['calTemp'], index=m['distance'][:,0], columns=times)

def plot_dts(dts_dataframe, min_dist=None, max_dist=None, min_time=None, max_time=None):
    pyplot.figure()
    plotax = pyplot.axes()
    pltdata = dts_dataframe.loc[min_dist:max_dist,min_time:max_time]
    myplot = plotax.pcolorfast(pltdata.astype(dtype=float))
    pyplot.colorbar(myplot, ax=plotax)
    xlocs, xlabels = pyplot.xticks()
    #For some reason an extra tick is created beyond the end of the data. Remove it using [:-1].
    xlocs, xlabels = xlocs[:-1], xlabels[:-1]
    xdates = pltdata.iloc[0,xlocs].index
    xdates = [ts.strftime('%Y-%m-%d %H:%M') for ts in xdates]
    ylocs, ylabels = pyplot.yticks()
    ylocs, ylabels = ylocs[:-1], ylabels[:-1]
    ydists = pltdata.iloc[ylocs,0].index
    pyplot.xticks(xlocs, xdates, rotation=45)
    pyplot.yticks(ylocs, ydists)
    return pltdata

def read_and_plot_pwrrat(dtsdir, **kwargs):
    '''
    Read Stokes and AntiStokes columns, then plot 1/(ln(S/AS)).
    Works on a single directory only as currently written.
    '''
    stokes = read_ddf_col(dtsdir, ddf_column=1).astype(dtype=float)
    antistokes = read_ddf_col(dtsdir, ddf_column=2).astype(dtype=float)
    ratio = 1/numpy.log(stokes/antistokes)
    ratio_trimmed = plot_dts(ratio, **kwargs)
    return ratio_trimmed 
