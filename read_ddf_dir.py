# -*- coding: utf-8 -*-
import os, datetime
import pandas
from matplotlib import pyplot

pyplot.ion() #set up interactive plotting
datadir='/home/aaron/work/phd/3ddts/SN411049/data/full data set/MammothKnitting/channel 4/2013/jan/'

def read_dts_dir(datadir, channel=4):
    for filepath in os.listdir(datadir):
        if filepath.endswith('.ddf'):
            timestamp = datetime.datetime.strptime(filepath, 'channel ' + str(channel) + ' %Y%m%d %H%M%S 00001.ddf')
            data_for_timestamp = pandas.read_csv(
                datadir + filepath,
                skiprows = 25,
                delimiter ='\t')
            data_for_timestamp.set_index('length (m)', inplace=True)
            # would be better to do data_for_timestamp['temperature (°C)'] but the degree symbol is a problem
            try:
                data[timestamp] = data_for_timestamp.ix[:,0]
            except NameError:
                # first iteration
                data = pandas.DataFrame(data_for_timestamp.ix[:,0],columns=[timestamp])
            print("finished reading in " + filepath)
    data.sort(axis=1, inplace=True)
    return data

def plot_dts(dts_dataframe, min_dist=None, max_dist=None, min_time=None, max_time=None):
    plotax = pyplot.axes()
    myplot = plotax.pcolorfast(dts_dataframe.ix[min_dist:max_dist,min_time:max_time])
    pyplot.colorbar(myplot, ax=plotax)
    pyplot.show()

plot_dts(dts_dataframe = read_dts_dir(datadir),
         min_dist = 300,
         max_dist = 600,
         min_time = '2013-01-03 13:31:07',
         max_time = '2013-01-03 13:31:07')
