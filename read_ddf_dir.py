# -*- coding: utf-8 -*-
import pandas, os, datetime
 
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
    return data

