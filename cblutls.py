"""Cable Utilities: Tools for integrating DTS and TLS data """

__author__ = "Aaron Curtis <aarongc@nmt.edu>"
__date__ = "2015-07-09"
__version__ = '0.1'

import pandas, numpy
import glob2
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import pyplot, colors
from mpl_toolkits.mplot3d import Axes3D
try:
    from mayavi import mlab
except:
    print("Mayavi not found. Fancy plotting will be unavailable.")

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

class Cable():
    '''
    Represents a DTS cable, made up of CableSections.
    Polyline files should end in .txt, and distance reference files should end in .dtr.
    See cblutls.CableSection class for details of file formats.
    '''
    def __init__(self, polyline_dirpath, dts_data_dirpath=None, **kwargs):
        self.sections=[]
        for polyline_filepath in glob2.glob(polyline_dirpath + '*.txt'):
            self.sections.append(CableSection(
                polyline_filepath,
                polyline_filepath.replace('.txt','.dtr'),
                **kwargs))

    def plot_w_mpl(self, ax=None):
        for section in self.sections:
            if not ax:
                f = pyplot.figure()
                ax = Axes3D(f)
            section.plot_w_mpl(ax)
            
    def plot_w_mayavi(self):
        for section in self.sections:
            section.plot_w_mayavi(colorbar=False)
        mlab.colorbar()
    
    def get_dts_data(self):
        '''
        Combines the dts data from all the cable sections into one dataframe
        '''
        dd = self.sections[0].dts_data
        for section in self.sections[1:]:
            dd = dd.combine_first(section.dts_data)
        return dd

class CableSection():
    '''
    Represents a section of DTS cable in 3D space
    
    Polyline file should be space-delimited, with one header row like:
        x y z
        
    Distance reference file should be space-delimited, NaNs marked as NaN with one header row like:
        x y z cable_dist fiber_dist description
    
    dts_data is the output of dts.read_dts_dirs
    
    cbl_fbr_b and cbl_fbr_m are the parameters for a linear equation <fbr>=m<cbl>+b that relates cable distance to fiber distance
    
    coord_cbl_m is the ratio of the coordinate system units (meters if UTM) to the cable length units (often feet)
    '''
    def __init__(self, polyline_filepath, dist_ref_pts_filepath, dts_data=None, cbl_fbr_m=None, cbl_fbr_b=None, coord_cbl_m=0.3048):
        self.dts_data = dts_data.copy()
        self.polyline_filepath = polyline_filepath
        self.dist_ref_pts_filepath = dist_ref_pts_filepath
        
        #Read in the output of InnovMetric IMSurvey's "export polyline to text"
        xyz = pandas.read_csv(polyline_filepath, sep=" ", comment="#")
        xyz['is_distref'] = False
        #Remove duplicates -- they cause problems later
        xyz.drop_duplicates(inplace=True)
        
        #Read in cable distance reference points
        dist_ref_pts = pandas.read_csv(dist_ref_pts_filepath, sep=" ", comment="#")
        
        #Insert each distance reference point between the 2 nearest cable polyline points
        for pt in dist_ref_pts.values: #TODO rewrite using apply()
            closest2 = find_2closest_points(xyz, pt)
            upper = xyz.ix[:closest2.index[1]]
            lower = xyz.ix[closest2.index[0]:]
            pt = pandas.DataFrame([pt], columns=dist_ref_pts.columns)
            pt['is_distref']=True
            xyz = pandas.concat([upper, pt, lower]).reset_index(drop=True)

        #Add i, j, k, euc_dist and cum_euc_dist columns
        self.data = pts_to_vectors(xyz)
        self.data.set_index('cum_euc_dist', inplace=True)
        
        #Add fiber distance for each polyline point
        self.interp_dists(coord_cbl_m)
        if cbl_fbr_m and cbl_fbr_b:
            self.set_fiber_d_frm_cable_d(cbl_fbr_m, cbl_fbr_b)
            if self.dts_data is not None:
                #calculate locations of dts points
                for dim in ['x','y','z']:
                    #dts_data.index is the fiber lengths
                    self.data.sort('fiber_dist',inplace=True)
                    spl = InterpolatedUnivariateSpline(self.data.fiber_dist.values, self.data[dim].values, k=1)
                    eval_at = self.dts_data[self.data.fiber_dist.min(): self.data.fiber_dist.max()].index.values
                    self.dts_data[dim] = numpy.nan
                    self.dts_data.ix[eval_at,dim] = spl(eval_at)

    def set_fiber_d_frm_cable_d(self, m, b):
        #TODO check that this doesn't wipe out the fiber dist ref points
        self.data.fiber_dist = self.data.cable_dist*m+b

    def plot_w_mpl(self, ax=None, cbl_dst_lbls=False):
        if not ax:
            f = pyplot.figure()
            ax = Axes3D(f)
        ax.plot(self.data.x.values, self.data.y.values, self.data.z.values)
        dr = self.get_distrefs()
        ax.plot(dr.x.values, dr.y.values, dr.z.values, 'r.')
        #find the non-nan dts values; i.e. the ones we've located on the polyline
        dtsnn = self.dts_data[self.dts_data.x.notnull()]
        #Plot average DTS temperatures
        time_averaged_dts = dtsnn.drop(['x','y','z'],axis='columns').mean(axis='columns').values
        time_averaged_dts_normalized = colors.Normalize(
            min(time_averaged_dts),
            max(time_averaged_dts))(time_averaged_dts)
        ax.scatter(
            dtsnn.x.values,
            dtsnn.y.values,
            dtsnn.z.values,
            marker='.',
            s=160,
            edgecolors='face',
            c=pyplot.cm.jet(time_averaged_dts_normalized))
        if cbl_dst_lbls:
            for idx, row in dtsnn.iterrows():    
                ax.text(
                    row.x,
                    row.y,
                    row.z,
                    idx
                )
        pyplot.show()

    def plot_w_mayavi(self, colorbar=True):
        #Plot average DTS temperatures
        dtsnn = self.dts_data[self.dts_data.x.notnull()]
        time_averaged_dts = dtsnn.drop(['x','y','z'],axis='columns').mean(axis='columns').values
        dts_plot = mlab.plot3d(dtsnn.x.values, dtsnn.y.values, dtsnn.z.values, time_averaged_dts, tube_radius=0.1)
        if colorbar:
            mlab.colorbar(dts_plot)
        dr = self.get_distrefs()
        mlab.points3d(dr.x.values, dr.y.values, dr.z.values, scale_factor=0.2)
    
    def get_distrefs(self):
        return self.data[self.data.is_distref==True]
    
    def interp_dists(self, coord_cbl_m=None):
        '''
        Add a cable distance for each polyline point, interpolated from the reference distances        
        '''
        distrefs = self.get_distrefs()
        if len(distrefs) > 1:
            #Create a spline function cbl_dist = spl(euc_dist)
            spl = InterpolatedUnivariateSpline(distrefs.index.values, distrefs.cable_dist.values, k=1)
            self.data.cable_dist = spl(self.data.index.values)
        else:
            #There's only one distance reference point for this section, so we correct offset based on data
            #and slope based on the ratio of the spatial coordinates (meters, if UTM) and cable (usually feet)
            offset = distrefs.cable_dist.values[0] - distrefs.index.values[0] * coord_cbl_m
            self.data.cable_dist = self.data.index.values * coord_cbl_m + offset

def interpolate_temperatures(cable_section, grid_step):
    from pykrige.k3d import Krige3D
    #Next two lines are also in CableSection.plot_w_mayavi -- refactor if it gets used a third time
    dtsnn = cable_section.dts_data[cable_section.dts_data.x.notnull()]
    time_averaged_dts = dtsnn.drop(['x','y','z'],axis='columns').mean(axis='columns').values    
    k = Krige3D(dtsnn.x.values, dtsnn.y.values, dtsnn.z.values, time_averaged_dts)
    grid = []
    for dim in ['x','y','z']:
        grid.append(numpy.arange(dtsnn[dim].min(), dtsnn[dim].max(), grid_step))
    kvals, sigmasq = k.execute('grid',grid[0],grid[1],grid[2])
    return kvals

def find_baths(dts_data, trefs):
    '''
    For each thermistor, correlate against all distances and plot the pearson coefficient values.
    '''
    tref_results={}
    for tref_name, tref in trefs.iteritems():
        results={}
        for dist, temps in dts_data.iterrows():
            results.update({dist:tref.corr(temps)})
        tref_results.update({tref_name:results})
    return tref_results