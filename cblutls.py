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
    return pandas.read_csv(filepath, sep=" ", comment="#", names=['x','y','z']) #TODO import in segments

def points_to_vectors(xyz):
    '''
    Convert a matrix of x, y, z points to a matrix of normalized vectors with
    rows of [xdiff, ydiff, yzdiff, distance_from_last_point, cumulative distance]
    '''
    ijk = xyz.diff()
    ijk.iloc[0] = xyz.iloc[0] #diff leaves NaNs in the first row, we want the position vector
    
    #Calculate euclidean distance from previous point
    dists = numpy.sqrt(numpy.sum(ijk**2, axis='columns'))
    cum_dists = numpy.cumsum(dists)
    
    #Normalize the vectors and rename the columns
    ijk_normd = ijk.divide(dists,axis='rows').columns(['i','j','k'])
    
    #Store distance and cumulative distance from first point along cable
    ijk['dist'] = dists
    ijk['cum_dists'] = numpy.cumsum(dists)
    
    #Return a rather redundant matrix containing:
    ## xyz, ijk, distances, cumulative distances
    return pandas.concat([xyz, ijk], axis='columns')

def find_closest_point(xyz_array, xyz_point):
    #Calculate distances between each point in xyz_array and xyz_point
    dists_squared = numpy.sum((xyz_array-xyz_point)^2)
    return numpy.argmin(dists_squared)

def add_dist_control_points(xyz, control_file):
    '''
    control_file: a list of points that represent known cable distance / location pairs,
    in the format x, y, z, dist
    '''
    closest_point_indices = xyz.apply(find_closest_point)

