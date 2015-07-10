"""Cable Utilities: Tools for integrating DTS and TLS data """

__author__ = "Aaron Curtis <aarongc@nmt.edu>"
__date__ = "2015-07-09"
__version__ = '0.1'

import pandas, numpy, ipdb

def pts_to_vectors(xyz):
    '''
    Convert a matrix of x, y, z points to a matrix of normalized vectors with
    rows of [xdiff, ydiff, yzdiff, distance_from_last_point, cumulative distance]
    '''
    ijk = xyz[['x','y','z']].diff()
    #ijk.iloc[0] = xyz.iloc[0] #diff leaves NaNs in the first row, we want the position vector
    
    #Calculate euclidean distance from previous point
    dists = numpy.sqrt(numpy.sum(ijk**2, axis='columns'))
    cum_dists = numpy.cumsum(dists)
    
    #Normalize the vectors and rename the columns
    ijk = ijk.divide(dists,axis='rows')
    ijk.columns = ['i','j','k']
    
    #Store distance and cumulative distance from first point along cable
    ijk['euc_dist'] = dists

    #Return a rather redundant matrix containing:
    ## xyz, ijk, distances, cumulative distances
    return pandas.concat([xyz, ijk], axis='columns')
    
def find_2closest_points(xyz_array, xyz_point):
    #Calculate distances between each point in xyz_array and xyz_point
    dists_squared = numpy.sum((xyz_array-xyz_point[:3])**2, axis='columns')
    return xyz_array.ix[numpy.argpartition(dists_squared,2)[:2]]

def find_loc_btw_2pts(a, b, dist):
    a = numpy.array(a)
    b = numpy.array(b)
    ba = b - a
    ba_normalized = ba/numpy.linalg.norm(ba)
    return a + ba_normalized * dist

def add_dist_control_pts(xyz, control_file):
    '''
    control_file: a list of points that represent known cable distance / location pairs,
    in the format x, y, z, dist
    '''
    closest_point_indices = xyz.apply(find_closest_point)

class CableSection():
    '''Represents a section of DTS cable in 3D space'''
    def __init__(self, polyline_filepath, dist_ref_pts_filepath):
        #Read in the output of InnovMetric IMSurvey's "export polyline to text"
        xyz = pandas.read_csv(polyline_filepath, sep=" ", comment="#", names=['x','y','z'])
        
        #Read in cable distance reference points
        dist_ref_pts = pandas.read_csv(dist_ref_pts_filepath, sep=" ", names=['x', 'y', 'z', 'ref_dist'], index_col=False)
        
        for pt in dist_ref_pts.values: #TODO rewrite using apply()
            closest2 = find_2closest_points(xyz, pt)
            #insert point between those two
            upper = xyz.ix[closest2.index[0]:]
            lower = xyz.ix[:closest2.index[1]]
            pt = pandas.DataFrame([pt], columns=['x','y','z','dist']) #TODO this is bad
            pt['is_distref']=True
            xyz = pandas.concat([upper, pt, lower]).reset_index(drop=True)
        self.xyzijk = pts_to_vectors(xyz)