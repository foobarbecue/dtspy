"""Cable Utilities: Tools for integrating DTS and TLS data """

__author__ = "Aaron Curtis <aarongc@nmt.edu>"
__date__ = "2015-07-09"
__version__ = '0.1'

import pandas, numpy

def find_loc_btw_2pts(a, b, dist):
    a = numpy.array(a)
    b = numpy.array(b)
    ba = b - a
    ba_normalized = ba/numpy.linalg.norm(ba)
    return a + ba_normalized * dist

def read_imsurvey_polyline(filepath):
    '''
    Read in the output of InnovMetric IMSurvey's "export polyline to text"
    '''
    return pandas.read_csv(filepath, sep=" ", comment="#", names=['x','y','z'])

def points_to_vectors(xyz):
    '''
    Convert a matrix of x, y, z points to a matrix of normalized vectors with
    rows of [xdiff, ydiff, yzdiff, distance_from_last_point, cumulative distance]
    '''
    ijk = xyz.diff()
    ijk.iloc[0] = xyz.iloc[0] #diff leaves NaNs in the first row, we want the position vector
    
    #Calculate euclidean distance from previous point
    dists = numpy.sqrt(numpy.sum(numpy.square(ijk), axis=1))
    cum_dists = numpy.cumsum(dists)
    
    #Normalize the vectors
    ijk_normd = ijk.divide(dists,axis='rows')
    
    #Store distance and cumulative distance from first point along cable
    ijk['dist'] = dists
    ijk['cum_dists'] = numpy.cumsum(dists)
    
    return ijk