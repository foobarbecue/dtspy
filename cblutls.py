"""Cable Utilities: Tools for integrating DTS and TLS data """

__author__ = "Aaron Curtis <aarongc@nmt.edu>"
__date__ = "2015-07-09"
__version__ = '0.1'

import pandas, numpy, ipdb
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

def pts_to_vectors(xyz):
    '''
    Convert a matrix of x, y, z points to a matrix of normalized vectors with
    rows of [xdiff, ydiff, yzdiff, distance_from_last_point, cumulative distance]
    '''
    ijk = xyz[['x','y','z']].diff()
    #ijk.iloc[0] = xyz.iloc[0] #diff leaves NaNs in the first row, we want the position vector
    
    #Calculate euclidean distance from previous point
    dists = numpy.sqrt(numpy.sum(ijk**2, axis='columns'))
    
    #Normalize the vectors and rename the columns
    ijk = ijk.divide(dists,axis='rows')
    ijk.columns = ['i','j','k']
    
    #Store distance and cumulative distance from first point along cable
    ijk['euc_dist'] = dists
    ijk['cum_euc_dist'] = numpy.cumsum(dists)
    ijk.cum_euc_dist[0] = 0.
    
    #Return a rather redundant matrix containing:
    ## xyz, ijk, distances, cumulative distances
    return pandas.concat([xyz, ijk], axis='columns')
    
def find_2closest_points(xyz_array, xyz_point):
    #Calculate distances between each point in xyz_array and xyz_point
    dists_squared = numpy.sum((xyz_array[['x','y','z']]-xyz_point[:3])**2, axis='columns')
    return xyz_array.ix[numpy.argpartition(dists_squared,2)[:2]]

def find_loc_btw_2pts(a, b, dist):
    a = numpy.array(a)
    b = numpy.array(b)
    ba = b - a
    ba_normalized = ba/numpy.linalg.norm(ba)
    return a + ba_normalized * dist

class CableSection():
    '''Represents a section of DTS cable in 3D space'''
    def __init__(self, polyline_filepath, dist_ref_pts_filepath, extrapolate=False, dts_data=None):
        #Read in the output of InnovMetric IMSurvey's "export polyline to text"
        xyz = pandas.read_csv(polyline_filepath, sep=" ", comment="#", names=['x','y','z'])
        xyz['is_distref'] = False
        
        #Read in cable distance reference points
        dist_ref_pts = pandas.read_csv(dist_ref_pts_filepath, sep=" ", names=['x', 'y', 'z', 'ref_dist'], index_col=False)
        
        for pt in dist_ref_pts.values: #TODO rewrite using apply()
            closest2 = find_2closest_points(xyz, pt)
            #insert point between those two
            upper = xyz.ix[:closest2.index[1]]
            lower = xyz.ix[closest2.index[0]:]
            pt = pandas.DataFrame([pt], columns=['x','y','z','cable_dist']) #TODO this is bad
            pt['is_distref']=True
            xyz = pandas.concat([upper, pt, lower]).reset_index(drop=True)
        self.data = pts_to_vectors(xyz)
        self.data.set_index('cum_euc_dist', inplace=True)
        self.interp_dists()
        if extrapolate:
            self.extrap_dists()

    def plot(self):
        f = pyplot.figure()
        ax = Axes3D(f)
        ax.plot(self.data.x.values, self.data.y.values, self.data.z.values)
        dr = self.get_distrefs()
        ax.plot(dr.x.values, dr.y.values, dr.z.values, 'r.')
        pyplot.show()
    
    def get_distrefs(self):
        return self.data[self.data.is_distref==True]
    
    def interp_dists(self):
        '''
        Add a cable distance for each polyline point, interpolated from the reference distances        
        TODO: solution for data that's not between two reference sections
        '''
        distrefs = self.get_distrefs()
        #We're using the euclidean distance from TLS polylines as the index
        self.data.cable_dist = numpy.interp(
            x = self.data.index.values,
            xp = distrefs.index.values,
            fp=distrefs.cable_dist.values,
            left=numpy.nan, right=numpy.nan)