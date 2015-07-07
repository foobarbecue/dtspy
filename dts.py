# -*- coding: utf-8 -*-
import os, datetime
import pandas
import numpy
from matplotlib import pyplot

pyplot.ion() #set up interactive plotting
datadirs=['/home/aaron/work/phd/3ddts/SN411049/data/full data set/MammothKnitting/channel 4/2013/jan/',
          '/home/aaron/work/phd/3ddts/SN411049/data/full data set/MammothKnitting/channel 4/2012/dec/']

def read_dts_dirs(datadirs=datadirs, ddf_column=0, channel=4):
    for datadir in datadirs:
        for filepath in os.listdir(datadir):
            if filepath.endswith('.ddf') and 'channel ' + str(channel) in filepath:
                timestamp = datetime.datetime.strptime(filepath, 'channel ' + str(channel) + ' %Y%m%d %H%M%S 00001.ddf')
                data_for_timestamp = pandas.read_csv(
                    datadir + filepath,
                    skiprows = 25,
                    delimiter ='\t')
                data_for_timestamp.set_index('length (m)', inplace=True)
                try:
                    data[timestamp] = data_for_timestamp.ix[:,ddf_column] # would be better to do data_for_timestamp['temperature (°C)'] but the degree symbol is a problem
                except NameError:
                    # first iteration
                    data = pandas.DataFrame(data_for_timestamp.ix[:,ddf_column],columns=[timestamp])
                    ref_temps = pandas.DataFrame
                #print("finished reading in " + filepath)
    data.sort(axis=1, inplace=True)
    return data

def plot_dts(dts_dataframe, min_dist=None, max_dist=None, min_time=None, max_time=None):
    plotax = pyplot.axes()
    myplot = plotax.pcolorfast(dts_dataframe.ix[min_dist:max_dist,min_time:max_time])
    pyplot.colorbar(myplot, ax=plotax)
    locs, labels = pyplot.xticks()
    #For some reason an extra tick is created beyond the end of the data. Remove it using [:-1].
    locs, labels = locs[:-1], labels[:-1]
    xdates = dts_dataframe.iloc[0,locs].index
    pyplot.xticks(locs, xdates, rotation=45)

temperature = read_dts_dirs(ddf_column=0, channel=4)
#stokes = read_dts_dirs(ddf_column='Stokes', channel=4)
#anti_stokes = read_dts_dirs(ddf_column='anti-Stokes', channel=4)


plot_dts(dts_dataframe = temperature,
        min_dist = 10,
        max_dist = 690,
        min_time = 0,
        max_time = -60)