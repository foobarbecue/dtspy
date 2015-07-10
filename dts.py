'''Reads in output files produced by fiber-optic distributed temperature sensing machines and plots them'''
__author__ = "Aaron Curtis <aarongc@nmt.edu>"
__version__ = '0.1'

# -*- coding: utf-8 -*-
import os, datetime
import pandas
import numpy
from matplotlib import pyplot
import ipdb

pyplot.ion() #set up interactive plotting

def read_dts_dirs(datadirs, ddf_column=0, channel=1):
    for datadir in datadirs:
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

def plot_dts(dts_dataframe, min_dist=None, max_dist=None, min_time=None, max_time=None):
    plotax = pyplot.axes()
    myplot = plotax.pcolorfast(dts_dataframe.loc[min_dist:max_dist,min_time:max_time].astype(dtype=float))
    pyplot.colorbar(myplot, ax=plotax)
    xlocs, xlabels = pyplot.xticks()
    #For some reason an extra tick is created beyond the end of the data. Remove it using [:-1].
    xlocs, xlabels = xlocs[:-1], xlabels[:-1]
    xdates = dts_dataframe.iloc[0,xlocs].index
    ylocs, ylabels = pyplot.yticks()
    ylocs, ylabels = ylocs[:-1], ylabels[:-1]
    ydists = dts_dataframe.iloc[ylocs,0].index
    pyplot.xticks(xlocs, xdates, rotation=45)
    pyplot.yticks(ylocs, ydists)