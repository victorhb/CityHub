import sys
import pickle
import utm
import operator
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Point
import numpy as np
import osmnx as ox
import networkx as nx
import math
import folium
import branca
from sklearn.neighbors import BallTree
from statistics import mean

EARTH_RADIUS = 6371.0

class CityHub:
     """CityHub class for urban data integration.
     
     Attributes:
         city_street_graph (Graph): A networkX undirected graph representing the city street graph.
         city_vert_list (list): List of vertices coordinates tuples of city street graph.
         city_vert_dict (dict): Dictionary of vertices coordinates of city street graph. Keys are coordinates, and values are its city_vert_list indices.
         city_tree (BallTree): BallTree of refined city street graph vertices.
         city_vert_ind_to_nxind_dict (dict): Dictionary to convert from city_vert_list indices to NetworkX indices
         city_vert_nxind_to_ind_dict (dict): Dictionary to convert from NetworkX indices to city_vert_list indices
         refined_city_vert_correspondence_dict (dict): Dictionary of correspondences between original and refined city street graph vertices
         refined_city_vert_list (list): List of refined vertices coordinates tuples of city street graph.
         refined_city_vert_dict (dict): Dictionary of refined vertices coordinates of city street graph. Keys are coordinates, and values are its refined_city_vert_list indices.
         
         RDLayers (list): The list of Regional Domain (RD) Layers.
         RDLayers_vert_list (dict): Dictionary of lists of RD Layers vertices coordinates
         RDLayers_balltree (BallTree): Dictionary of BallTrees of RDLayers
         RDLayers_area_vertices_indices_dict (dict): Dictionary of dictionaries for each RD layer, where the key is an area index and the value is a list with the indices of its vertices.
         RDLayers_area_vertices_coords_dict (dict): Dictionary of dictionaries for each RD layer, where the key is an area index, and the value is a list with the coordinates of its vertices.
         RDLayers_vertices_area_dict (dict): Dictionary of dictionaries for each RD layer, where the key is a vertex index, and the value is its area
         PBLayers (list): The list of Point-based (PB) Layers.
         PBLayers_balltree (dict): Dictionary of BallTrees of PBLayers
         PBLayers_csv_lat_key_column (dict): Dictionary of the key latitude column of each PB Layer.
         PBLayers_csv_lng_key_column (dict): Dictionary of the key longitude column of each PB Layer.   
         PBLayers_corner_projection_dicts (dict): Dictionary of dictionaries, being one for each corner of the city street graph, holding points from the corresponding PBLayer that are projected to that corner.
         PALayers_mesh (list): The list of the meshes of the Polygon-Aggregated (PA) layers.
         PALayers_mesh_keycolumn (dict): Dictionary of the key column of each PA layer mesh.
         PALayers_csv_data (dict): Dictionary of dataframes of each PA layer, holding data from CSV files.  
         PALayers_csv_keycolumn (dict):  Dictionary of the key column of each PA layer CSV data.
         SMLayers_known_measurements (dict): Dictionary with all time series for each measurement point
         SMLayers_temp_agg_funcs (dict): Dictionary with the temporal aggregation function for each measurement variable
         SMLayers_measurements_points_info (DataFrame): Dataframe with information (e.g. lat, lng) of each measurement point      
         SMLayers_estimated_measurements (DataFrame): DataFrame with estimations for all measurement variables for each corner
         
     """

     def __init__(self, filename, EDGE_LENGTH_THRESHOLD=40.0):
        """
        Load a city mesh from a file or city query string. GPICKLE (from nx.read_gpickle()) formats are accepted. Then, data structures are generated.
    
        Parameters:
            filename (string): city mesh filename, including the file extension, which will be used to identify the file format.
            EDGE_LENGTH_THRESHOLD (float):maximum edge length for robust nearest neighbor search. An edge will be subdivided to 
            assure its length is less than the threshold.
        
    """
        
        
        #city streets data
        self.refined_city_vert_correspondence_dict = dict()        
        self.refined_city_vert_list = []
        self.refined_city_vert_correspondence_list = [] 
        
        
        ##layers' data
        
        # Regional-Domain (RD) Layer 
        self.RDLayers = []
        self.RDLayers_balltree = dict()       
        self.RDLayers_vert_list = dict()         
        self.RDLayers_area_vertices_indices_dict = dict()
        self.RDLayers_area_vertices_coords_dict = dict()
        self.RDLayers_vertices_area_dict = dict()
        
        # Point-Based (PB) Layer         
        self.PBLayers = []
        self.PBLayers_balltree = dict()
        self.PBLayers_csv_lat_key_column = dict()
        self.PBLayers_csv_lng_key_column = dict()
        self.PBLayers_corner_projection_dicts = dict()
        
        # Polygon-Aggregated (PA) Layer
        self.PALayers_mesh = []
        self.PALayers_csv_data = dict()
        self.PALayers_mesh_keycolumn = dict()
        self.PALayers_csv_keycolumn = dict()
        self.PALayers_polygon_vertices_dict = dict()
        self.PALayers_aggregated_polygon_vert_list = dict()
        self.PALayers_aggregated_polygon_vert_dict = dict()
        self.PALayers_aggregated_polygon_vertices_indices_dict = dict()
        self.PALayers_aggregated_polygon_vert_setcens_list = dict()
        self.PALayers_refined_aggregated_polygon_vert_list = dict()
        self.PALayers_refined_aggregated_polygon_vert_dict = dict()
        self.PALayers_refined_aggregated_polygon_vert_correspondence_dict = dict()
        self.PALayers_aggregated_polygon_tree = dict()
        
        # Sparse-Measurement (SM) Layer
        self.SMLayers_known_measurements = dict()
        self.SMLayers_temp_agg_funcs = dict()
        self.SMLayers_measurements_points_info = pd.DataFrame()
        self.SMLayers_estimated_measurements = pd.DataFrame()
        
        #optional feature vectors data
        self.feature_vecs = []
        
        #visualization variables
        self.nearest_vertex_index = 0
        self.nearest_subnormal_list = []

        if not self.load_city_street_graph(filename):
            print('Cannot open '+filename)
            return
        self.preprocess_city_mesh(True,EDGE_LENGTH_THRESHOLD)
    
     def save_preprocessed_CityHub(self, filename):
        """
        Save preprocessed data structures of a city street graph to a pickle file .
    
        Parameters
        ----------
        filename : string
            Pickle filename to be written.
        """
        pickle.dump(self, open(filename, 'wb'))

     def load_city_street_graph(self, city_string):
        """
        Load a city street graph from OMSnx library city query string, or a gpickle.
        In case of success, the Graph city_street_graph will be created.
    
        Parameters
        ----------
        city_string (string): city mesh filename in gpickle format (from nx.read_gpickle()), including the file 
        extension, or a string with the city description to download from OSMnx (i.e. 'Sao Paulo, Brazil')
        
        Returns
        -------
            returns True if the city street graph is sucessfully loaded from file.
        """
        
        print('Loading file...')
                                              
        extension_str = city_string[-3:]
        if(extension_str.lower()=='kle'):
            try:
                self.city_street_graph = nx.read_gpickle(city_string)
            except:
                return False
        else:
            try:
                self.city_street_graph = ox.graph_from_place(city_string, network_type='drive')
            except:
                return False
            
        self.city_street_graph=self.city_street_graph.to_undirected()
    
        return True  

    
     def preprocess_city_mesh(self, build_tree=True, EDGE_LENGTH_THRESHOLD=sys.float_info.max):
        """
        Generate data structures to quickly retrieve relevant information. A valid city_street_graph is required.
    
        Parameters
        ----------
        build_tree : bool
            builds a tree for querying
        EDGE_LENGTH_THRESHOLD: float
            maximum edge length for robust nearest neighbor search. An edge will be subdivided to assure its length is less than the threshold.
        
        Returns
        -------
            returns True if the preprocessing succeeds.
        """   
        
        print('Preprocessing city mesh...')
                                              
        l=list(self.city_street_graph.nodes.data())
        self.city_vert_nxind_to_ind_dict = {k[0]: v for v, k in enumerate(l)}
        self.city_vert_ind_to_nxind_dict = {v: k[0] for v, k in enumerate(l)}
        self.city_vert_list = [tuple((node[1]['y'],node[1]['x'])) for node in l]
        self.city_vert_dict=dict(zip(self.city_vert_list,range(len(self.city_vert_list))))
        
        
        if(EDGE_LENGTH_THRESHOLD<sys.float_info.max-1):
            print('Refining mesh...')
            refined_city_vert_coords_correspondence_dict = dict()
            refined_city_vert_set = set()

            for u,v,a in self.city_street_graph.edges(data=True):
                vertA = self.city_vert_list[self.city_vert_nxind_to_ind_dict[u]]
                vertB = self.city_vert_list[self.city_vert_nxind_to_ind_dict[v]]
                great_circle_dist = ox.distance.great_circle_vec(vertA[0],vertA[1],vertB[0],vertB[1])
                latlong_dist = ox.distance.euclidean_dist_vec(vertA[0],vertA[1],vertB[0],vertB[1])
                
                refined_city_vert_set.add(vertA)
                refined_city_vert_coords_correspondence_dict[vertA] = [self.city_vert_nxind_to_ind_dict[u],self.city_vert_nxind_to_ind_dict[v]]
                if(great_circle_dist>EDGE_LENGTH_THRESHOLD):
                    latlong_threshold = EDGE_LENGTH_THRESHOLD * latlong_dist / great_circle_dist
                    line = shapely.geometry.LineString([(vertA[0],vertA[1]), (vertB[0],vertB[1])])
                    num_vert = max(math.ceil(line.length / latlong_threshold), 1)
                
                    for n in range(1,num_vert-1):
                        interpolated_coords = line.interpolate(n / num_vert, normalized=True).coords[0]
                        refined_city_vert_set.add(interpolated_coords)
                        refined_city_vert_coords_correspondence_dict[interpolated_coords] = [self.city_vert_nxind_to_ind_dict[u],self.city_vert_nxind_to_ind_dict[v]]
                refined_city_vert_set.add(vertB)
                refined_city_vert_coords_correspondence_dict[vertB] = [self.city_vert_nxind_to_ind_dict[u],self.city_vert_nxind_to_ind_dict[v]]
                
            self.refined_city_vert_list = list(refined_city_vert_set)
            self.refined_city_vert_dict={k: v for v, k in enumerate(self.refined_city_vert_list)}       
            
            for k in self.refined_city_vert_dict.keys():                
                ind_refined = self.refined_city_vert_dict[k]
                inds_original = refined_city_vert_coords_correspondence_dict[k] 
                self.refined_city_vert_correspondence_dict[ind_refined] = inds_original     
        else:
            self.refined_city_vert_list = self.vert_list.copy()
            self.refined_city_vert_dict = self.vert_dict.copy()
            self.refined_city_vert_correspondence_dict = self.refined_city_vert_dict.copy()
            
        if build_tree:
            self.city_tree = BallTree(np.deg2rad(np.c_[np.array(self.refined_city_vert_list)]), metric='haversine')
        
        return True


     def query_point_in_city_mesh(self, lat, long, return_nearest_index=False):
        """
        Query a point in lat-long format (in degrees), in the city street graph.
        refined_vert_list is used in the tree, but the original vertices will be returned.
    
        Parameters
        ----------
        lat: float
            latitude of the query point
        long: float
            longitude of the query point
        return_nearest_index: bool
            whether to return the nearest point index (according to vert_list indexing) or a lat-long tuple
            
        Returns
        -------
            int (index) or tuple (lat,long) representing the nearest point in the tree.
        """        
        
        try:
            self.city_tree
        except:
            print('BallTree does not exist. building...')
            self.build_tree()
            
        [dist,v] =self.city_tree.query(convert_deg_query(lat,long))
        
        v_cand1 = self.refined_city_vert_correspondence_dict[v[0][0]][0]
        v_cand2 = self.refined_city_vert_correspondence_dict[v[0][0]][1]
        
        v_cand1_coords = self.city_vert_list[v_cand1]
        v_cand2_coords = self.city_vert_list[v_cand2]
        
        great_circle_dist1 = ox.distance.great_circle_vec(lat,long,v_cand1_coords[0],v_cand1_coords[1])
        great_circle_dist2 = ox.distance.great_circle_vec(lat,long,v_cand2_coords[0],v_cand2_coords[1])                                                              
        if great_circle_dist1<great_circle_dist2:
             v = v_cand1                                                                
        else:
             v = v_cand2
                                                                              
        if return_nearest_index:
            return v
        else:
            return self.city_vert_list[v]
        
     def query_point_in_graph_radius_in_city_mesh(self, lat, long, query_radius):
        """
        Query a point in lat-long format (in degrees), in the city street graph. 
        First, the nearest graph vertex is found. Than, vertices which are less than 'radius' distance, using a distance through the graph's edges, are returned.
        
    
        Parameters
        ----------
        lat: float
            latitude of the query point
        long: float
            longitude of the query point
        query_radius: float
            Include all neighbors of distance<=radius from n, using the graph's distance in km
            
        Returns
        -------
            a list of points within maximum graph's distance to the query point
        """        
        
        try:
            self.city_tree
        except:
            print('BallTree does not exist. building...')
            self.build_tree()
            
        [dist,v] =self.city_tree.query(convert_deg_query(lat,long))
        
        v_cand1 = self.refined_city_vert_correspondence_dict[v[0][0]][0]
        v_cand2 = self.refined_city_vert_correspondence_dict[v[0][0]][1]
        
        v_cand1_coords = self.city_vert_list[v_cand1]
        v_cand2_coords = self.city_vert_list[v_cand2]
        
        great_circle_dist1 = ox.distance.great_circle_vec(lat,long,v_cand1_coords[0],v_cand1_coords[1])
        great_circle_dist2 = ox.distance.great_circle_vec(lat,long,v_cand2_coords[0],v_cand2_coords[1])                                                              
        if great_circle_dist1<great_circle_dist2:
             v = v_cand1                                                                
        else:
             v = v_cand2
        ps=nx.single_source_dijkstra_path(self.city_street_graph,self.city_vert_ind_to_nxind_dict[v],cutoff=query_radius*1000.0,weight='length')
        
    
        result_nodes = []
        for k in ps.keys():
            result_nodes.append(self.city_vert_nxind_to_ind_dict[k])    
        return result_nodes  

    
     def load_PALayer_mesh(self, filename, key_column = 'Name', EDGE_LENGTH_THRESHOLD=sys.float_info.max):
        """
        Load a Polygon-Aggregated Layer mesh from a file for polygon-aggregated data purposes, 
        such as census sectors data. KML and SHP formats are accepted.
    
        Parameters
        ----------
        filename : string
            city mesh filename, including the file extension, which will be used to identify the file format.
        key_column : string
            name of the key column of each polygon (or row of the dataframe)
        EDGE_LENGTH_THRESHOLD: float
            maximum edge length for robust nearest neighbor search. An edge will be subdivided to assure its length is less than the threshold.            
        Returns
        -------
            returns the layer index or -1 if it fails.
        """
        
        extension_str = filename[-3:]
        if(extension_str.lower()=='kml'):
            gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
            try:
                self.PALayers_mesh.append(gpd.read_file(filename, driver='KML'))
            except:
                return -1
        elif(extension_str.lower()=='shp'):
            try:
                self.PALayers_mesh.append(aggregated_polygon_mesh=gpd.read_file(filename, driver='SHP'))
            except:
                return -1
        self.PALayers_mesh_keycolumn[len(self.PALayers_mesh)-1]=key_column
        self.preprocess_PALayer_mesh(len(self.PALayers_mesh)-1,True,True,EDGE_LENGTH_THRESHOLD)
        return len(self.PALayers_mesh)-1
        

     def preprocess_PALayer_mesh(self, layer, build_tree=True, swap_coordinates=True, EDGE_LENGTH_THRESHOLD=sys.float_info.max):
        """
        Generate data structures to quickly retrieve relevant information from the polygon meshes of Polygon-Aggregated layers.
        A city street graph is required.
    
        Parameters
        ----------
        layer (int): layer position according to self.PALayers_mesh
        build_tree (bool): builds a BallTree for querying
        swap_coordinates (bool): useful when latitude and longitude coordinates are given in the wrong order.
        EDGE_LENGTH_THRESHOLD (float): maximum edge length for robust nearest neighbor search. An edge
        will be subdivided to assure its length is less than the threshold.
        
        Returns
        -------
            returns True if the preprocessing succeeds.
        """
        
        if len(self.PALayers_mesh)<=layer:
            return False
        
        if self.PALayers_mesh[layer].empty:
            return False
        
        X=0
        Y=1
        if swap_coordinates:
            X=1
            Y=0
        
        
        """ this code generate PALayers_polygon_vertices_dict[layer] from the kml dataframe (self.PALayers_mesh[layer])
        * PALayers_polygon_vertices_dict[layer] is a dictionary where the key is a aggregated-polygon code, and the value is a list with the coordinates of its vertices"""

        polygon_vertices_dict = dict()
        for ap in self.PALayers_mesh[layer].index:
            if(type(self.PALayers_mesh[layer].geometry[ap])==shapely.geometry.multipolygon.MultiPolygon):
                continue
            polygon_vertices_dict[self.PALayers_mesh[layer][self.PALayers_mesh_keycolumn[layer]][ap]] = np.dstack((self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[X],self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[Y])).tolist()[0]
        self.PALayers_polygon_vertices_dict[layer]=polygon_vertices_dict
        
    
        """ this code generates PALayers_aggregated_polygon_vert_list[layer], PALayers_aggregated_polygon_vert_dict[layer] from the original kml dataframe (self.PALayers_mesh[layer])
        * aggregated_polygon_vert_set is the set of vertices tuples
        * self.PALayers_aggregated_polygon_vert_list[layer] is an (indexed) list of tuples of coordinates of the set of (unique) vertices
        * self.PALayers_aggregated_polygon_vert_dict[layer] is a dictionary of vertices where the key is a tuple of lat/long coordinates and the value is the index of vert_list """

        aggregated_polygon_vert_set = set()
        for ap in self.PALayers_mesh[layer].index:
            if(type(self.PALayers_mesh[layer].geometry[ap])==shapely.geometry.multipolygon.MultiPolygon):
                continue
            polygon_vert_list=np.dstack((self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[X],self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[Y])).tolist()[0]
            for vert in polygon_vert_list:
                aggregated_polygon_vert_set.add(tuple(vert))
        self.PALayers_aggregated_polygon_vert_list[layer] = list(aggregated_polygon_vert_set)
        self.PALayers_aggregated_polygon_vert_dict[layer]={k: v for v, k in enumerate(self.PALayers_aggregated_polygon_vert_list[layer])}    

        """ this code generates PALayers_aggregated_polygon_vertices_indices_dict[layer] from vert_dict and the original kml dataframe (self.PALayers_mesh[layer])
        * PALayers_aggregated_polygon_vertices_indices_dict[layer] is a dictionary where the key is an aggregated-polygon code and the value is a list with the indices of its vertices."""

        self.PALayers_aggregated_polygon_vertices_indices_dict[layer] = dict()
        for ap in self.PALayers_mesh[layer].index:
            if(type(self.PALayers_mesh[layer].geometry[ap])==shapely.geometry.multipolygon.MultiPolygon):
              continue
            polygon_vert_list=np.dstack((self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[X],self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[Y])).tolist()[0]
            polygon_vert_ind_list=list()
            for vert in polygon_vert_list:
                polygon_vert_ind_list.append(self.PALayers_aggregated_polygon_vert_dict[layer][tuple(vert)])
            self.PALayers_aggregated_polygon_vertices_indices_dict[layer][self.PALayers_mesh[layer][self.PALayers_mesh_keycolumn[layer]][ap]] = polygon_vert_ind_list     

        """ this code generates PALayers_aggregated_polygon_vert_setcens_list[layer] from PALayers_aggregated_polygon_vert_dict[layer] and the original kml dataframe (self.PALayers_mesh[layer])
        * self.PALayers_aggregated_polygon_vert_setcens_list[layer] is a list comprised of sets of 1-ring polygons (represented with its code) of each vertex index."""

        self.PALayers_aggregated_polygon_vert_setcens_list[layer] = [set() for _ in range(len(self.PALayers_aggregated_polygon_vert_list[layer]))] 
        for ap in self.PALayers_mesh[layer].index:
            if(type(self.PALayers_mesh[layer].geometry[ap])==shapely.geometry.multipolygon.MultiPolygon):
                continue
            polygon_vert_list=np.dstack((self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[X],self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[Y])).tolist()[0]
            for vert in polygon_vert_list:  
                self.PALayers_aggregated_polygon_vert_setcens_list[layer][self.PALayers_aggregated_polygon_vert_dict[layer][tuple(vert)]].add(self.PALayers_mesh[layer][self.PALayers_mesh_keycolumn[layer]][ap])
                 
                    
        ''' this code will refine the aggregated polygon mesh to assure the edge lengths are not larger than EDGE_LENGTH_THRESHOLD. 
        * refined_aggregated_polygon_vert_set holds the new (refined) vert set.
        * self.PALayers_refined_aggregated_polygon_vert_list[layer] holds an (indexed) list of tuples of coordinates of the set of (unique) refined vertices
        * self.PALayers_refined_aggregated_polygon_vert_dict[layer] is a dictionary of refined vertices where the key is a tuple of lat/long coordinates and the value is the index of self.refined_vert_list
        * self.PALayers_refined_aggregated_polygon_vert_correspondence_dict[layer] is a dictionary of refined vertices where the key is an index of a refined vertex, and the value is the index to the nearest vertex of the edge that originated the refined vertex through interpolation.
        '''
        if(EDGE_LENGTH_THRESHOLD<sys.float_info.max-1):
            print('Refining mesh...')
            refined_aggregated_polygon_vert_set = set()
            refined_vert_coords_correspondence_dict = dict()
            
            for sector in self.PALayers_aggregated_polygon_vertices_indices_dict[layer].keys():
                polygon_vert_list=self.PALayers_aggregated_polygon_vertices_indices_dict[layer][sector]
                for i in range(len(polygon_vert_list)-1):
                    vertA = self.PALayers_aggregated_polygon_vert_list[layer][polygon_vert_list[i]]
                    vertB = self.PALayers_aggregated_polygon_vert_list[layer][polygon_vert_list[i+1]]
                    great_circle_dist = ox.distance.great_circle_vec(vertA[0],vertA[1],vertB[0],vertB[1])
                    latlong_dist = ox.distance.euclidean_dist_vec(vertA[0],vertA[1],vertB[0],vertB[1])
                    if(great_circle_dist>EDGE_LENGTH_THRESHOLD):
                        latlong_threshold = EDGE_LENGTH_THRESHOLD * latlong_dist / great_circle_dist
                        line = shapely.geometry.LineString([(vertA[0],vertA[1]), (vertB[0],vertB[1])])
                        num_vert = max(math.ceil(line.length / latlong_threshold), 1)
                        for n in range(num_vert):
                            interpolated_coords = line.interpolate(n / num_vert, normalized=True).coords[0]
                            if( not interpolated_coords in refined_aggregated_polygon_vert_set):
                                refined_aggregated_polygon_vert_set.add(interpolated_coords)
                                if (n<num_vert/2):
                                    refined_vert_coords_correspondence_dict[interpolated_coords] = polygon_vert_list[i]
                                else:
                                    refined_vert_coords_correspondence_dict[interpolated_coords] = polygon_vert_list[i+1]
                    else:
                        refined_aggregated_polygon_vert_set.add(vertA)
                        refined_aggregated_polygon_vert_set.add(vertB)
                        refined_vert_coords_correspondence_dict[vertA] = polygon_vert_list[i]
                        refined_vert_coords_correspondence_dict[vertB] = polygon_vert_list[i+1]
            
            self.PALayers_refined_aggregated_polygon_vert_list[layer] = list(refined_aggregated_polygon_vert_set)
            self.PALayers_refined_aggregated_polygon_vert_dict[layer]={k: v for v, k in enumerate(self.PALayers_refined_aggregated_polygon_vert_list[layer])}       
            self.PALayers_refined_aggregated_polygon_vert_correspondence_dict[layer] = dict()        
            for k in self.PALayers_refined_aggregated_polygon_vert_dict[layer].keys():                
                ind_refined = self.PALayers_refined_aggregated_polygon_vert_dict[layer][k]
                ind_original = refined_vert_coords_correspondence_dict[k] 
                self.PALayers_refined_aggregated_polygon_vert_correspondence_dict[layer][ind_refined] = ind_original
        else:
            self.PALayers_refined_aggregated_polygon_vert_list[layer] = self.PALayers_aggregated_polygon_vert_list[layer].copy()
            self.PALayers_refined_aggregated_polygon_vert_dict[layer] = self.PALayers_aggregated_polygon_vert_dict[layer].copy()
            self.PALayers_refined_aggregated_polygon_vert_correspondence_dict[layer] = self.PALayers_refined_aggregated_polygon_vert_dict[layer]
            
            
        if build_tree:
             self.PALayers_aggregated_polygon_tree[layer] = BallTree(np.deg2rad(np.array(self.PALayers_refined_aggregated_polygon_vert_list[layer])), metric='haversine')
        return True
    
    
     def load_PALayers_csv_data(self, layer, filename, polygon_key_column):
        """
        Load a CSV file made up of polygon-aggregated data and generate datadf dataframe.
    
        Parameters
        ----------
        layer (int): layer position according to self.PALayers_mesh.
        filename (string): CSV filename. 
        polygon_key_column (string): Name of the key column of the CSV file that matches the polygon's key column PALayers_mesh_keycolumn of the corresponding PA Layer.
        
        Returns
        -------
            returns True if the CSV file is succesfully loaded.
        """
        try:
            self.PALayers_csv_data[layer]=pd.read_csv(filename)
        except NameError as e:
             print(e)  
             return False
        
        if polygon_key_column in self.PALayers_csv_data[layer].columns:
            self.PALayers_csv_keycolumn[layer] = polygon_key_column
        else:
            print('invalid polygon_key_column.')
            return False
            
        return True    


     def query_point_in_PALayer_mesh(self, layer, lat, long, return_nearest_index=False):
        """
        Query a point in lat-long format, in degrees, in a PALayer mesh.
    
        Parameters
        ----------
        layer (float): layer position according to self.PALayers_mesh.
        lat: float
            latitude of the query point
        long: float
            longitude of the query point
        return_nearest_index: bool
            whether to return the nearest point index (according to vert_list indexing) or a lat-long tuple
            
        Returns
        -------
            int (index) or tuple (lat,long) representing the nearest point in the tree.
        """        
        
        try:
            self.PALayers_aggregated_polygon_tree[layer]
        except:
            print('BallTree does not exist. building...')
            self.PALayers_aggregated_polygon_tree[layer] = BallTree(np.deg2rad(np.array(self.PALayers_refined_aggregated_polygon_vert_list[layer])), metric='haversine')
            
        [dist,v] =self.PALayers_aggregated_polygon_tree[layer].query(convert_deg_query(lat,long))
        
        v_org = self.PALayers_refined_aggregated_polygon_vert_correspondence_dict[layer][v[0][0]]
        self.nearest_vertex_index=v_org
        if return_nearest_index:
            return v_org
        else:
            return self.PALayers_aggregated_polygon_vert_list[layer][v_org]


     def polygons_from_PALayer_vertex(self, layer, vert_ind):
        """
        Compute the polygons in the star of the aggregated polygon vertex 'vert_ind'.
    
        Parameters
        ----------
        layer (int): layer position according to self.PALayers_mesh
        vert_ind (int): vertex index (according to city_vert_list indexing)
        Returns
        -------
            a list of polygon keys, according to self.PALayers_mesh_keycolumn[layer], or an empty list in case of an issue
        """        
        try:
            self.PALayers_aggregated_polygon_vertices_indices_dict[layer]
        except:
            return []
            
        nearest_subnormal_list = [a for a, b in self.PALayers_aggregated_polygon_vertices_indices_dict[layer].items() if vert_ind in b]      

        return nearest_subnormal_list 
    
    
     def polygon_features_from_PALayer_vertex(self, layer, vert_ind, feature):
        """
        Retrieve a list of features from the the polygons in the star of the aggregated polygon vertex 'vert_ind'.
    
        Parameters
        ----------
        layer (int): layer position according to self.PALayers_mesh.
        vert_ind (int): vertex index (according to city_vert_list indexing)
        feature (string): feature name to be retrieved from each neighbour polygon.
        Returns
        -------
            a list of the features 'feature' from each neighbour polygon keys, or an empty list in case the CSV polygon data is not loaded.
        """
        
        if not feature in self.PALayers_csv_data[layer].columns:
            print('invalid feature name '+feature)
            return []
        
        polygon_list=self.polygons_from_PALayer_vertex(layer,vert_ind)      
        
        feature_list = []
        for polygon in polygon_list:
            dfsearch = self.PALayers_csv_data[layer].loc[self.PALayers_csv_data[layer][self.PALayers_csv_keycolumn[layer]]==int(polygon)]
            if dfsearch.shape[0]>0:
                feature_list.append(dfsearch[feature].to_numpy()[0])
        return feature_list
    
     def polygon_features_from_coords(self, layer, lat, long, feature):
        """
        Retrieve a list of features from the the polygons in the star of the nearest vertex to a query point with coordinates (lat,long), among the aggregated polygon vertices.
    
        Parameters
        ----------
        layer (int): layer position according to self.PALayers_mesh.
        lat (float): latitude of the query point.
        long (float): longitude of the query point.
        feature (string): feature name to be retrieved from each neighbour polygon.
        Returns
        -------
            a list of the features 'feature' from each neighbour polygon keys, or an empty list in case the CSV polygon data is not loaded.
        """        
        
        if not feature in self.PALayers_csv_data[layer].columns:
            print('invalid feature name '+feature)
            return []
        
        vert_ind=self.query_point_in_PALayer_mesh(layer,lat, long, True)
        polygon_list=self.polygons_from_PALayer_vertex(layer,vert_ind)      
        
        feature_list = []
        for polygon in polygon_list:
            feature_val = self.PALayers_csv_data[layer].loc[self.PALayers_csv_data[layer][self.PALayers_csv_keycolumn[layer]]==int(polygon)][feature].to_numpy()[0]
            feature_list.append(feature_val)
        return feature_list

    
     def load_RDLayer(self, filename, preprocess=True, swap_coordinates=True):
        """
        Load a Regional Domain Layer describing a specific area, such as subnormal agglomerates, from a mesh file.
        KML and SHP formats are accepted. 
        In case of success, a GeoDataFrame will be created and appended to the list self.RDLayers.
    
        Parameters
        ----------
        filename (string): mesh filename, including the file extension, which will be used to identify the file format.
        preprocess (bool): calls preprocess_RDLayer if True.
        swap_coordinates (bool): useful when latitude and longitude coordinates are given in the wrong order.
        
        Returns
        -------
            returns an integer representing the layer dataframe position in self.RDLayers if the mesh is
            sucessfully loaded from file, or -1 otherwise.
        """
        extension_str = filename[-3:]
        if(extension_str.lower()=='kml'):
            gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
            try:
                rdlayerdf = gpd.read_file(filename, driver='KML')
            except:
                return -1
        elif(extension_str.lower()=='shp'):
            try:
                rdlayerdf = gpd.read_file(filename, driver='SHP')
            except:
                return -1

        self.RDLayers.append(rdlayerdf)
        if preprocess:
            self.preprocess_RDLayer(len(self.RDLayers)-1, True, swap_coordinates)
        return len(self.RDLayers)-1


     def preprocess_RDLayer(self, layer, build_tree=True, swap_coordinates=True):
        """
        Generate data structures to quickly retrieve relevant information from RDLayer meshes. The layer must be loaded first.
    
        Parameters:
        layer (int): layer position according to self.RDLayers.
        build_tree (bool): builds a BallTree for querying.
        swap_coordinates (bool): useful when latitude and longitude coordinates are given in the wrong order.
        
        Returns:
            bool: returns True if the preprocessing succeeds.
        """
        
        if len(self.RDLayers)<=layer:
            return False
        
        X=0
        Y=1
        if swap_coordinates:
            X=1
            Y=0   
        
        """
        this code generates self.RDLayers_vert_list[layer] from the layer self.RDLayers[layer].
        * vert_set is the set of vertices tuples
        * vert_list is an (indexed) list of tuples of coordinates of the set of (unique) vertices
        * vert_dict is a dictionary of vertices where the key is a tuple of lat/long coordinates and the value is the index of vert_list
        """
        
        vert_set = set()
        for area in self.RDLayers[layer].index:
            if(type(self.RDLayers[layer].geometry[area])==shapely.geometry.multipolygon.MultiPolygon):
              continue    
            polygon_vert_list=np.dstack((self.RDLayers[layer].geometry[area].exterior.coords.xy[X],self.RDLayers[layer].geometry[area].exterior.coords.xy[Y])).tolist()[0]
            for vert in polygon_vert_list:
                vert_set.add(tuple(vert))
        vert_list = list(vert_set)   
        vert_dict={k: v for v, k in enumerate(vert_list)}
        
        self.RDLayers_vert_list[layer] = vert_list
        
        
        """
        this code generates self.RDLayers_area_vertices_indices_dict[layer] from self.RDLayers[layer]
        * self.RDLayers_area_vertices_indices_dict is a dictionary where the key is an area index and the value is a list with the indices of its vertices.
        """
        
        area_vertices_indices_dict = dict()
        for area in self.RDLayers[layer].index:
            if(type(self.RDLayers[layer].geometry[area])==shapely.geometry.multipolygon.MultiPolygon):
              continue
            polygon_vert_list=np.dstack((self.RDLayers[layer].geometry[area].exterior.coords.xy[X],self.RDLayers[layer].geometry[area].exterior.coords.xy[Y])).tolist()[0]
            polygon_vert_ind_list=list()
            for vert in polygon_vert_list:
                polygon_vert_ind_list.append(vert_dict[tuple(vert)])
            area_vertices_indices_dict[layer] = polygon_vert_ind_list

        self.RDLayers_area_vertices_indices_dict[layer]=area_vertices_indices_dict  
        
        """ this code generates self.RDLayers_area_vertices_coords_dict[layer] from self.RDLayers[layer]
        * self.RDLayers_area_vertices_coords_dict[layer] is a dictionary where the key is an area index, and the value is a list with the coordinates of its vertices
        """
            
        area_vertices_coords_dict = dict()
        for area in self.RDLayers[layer].index:
            if(type(self.RDLayers[layer].geometry[area])==shapely.geometry.multipolygon.MultiPolygon):
              continue
            area_vertices_coords_dict[area] = np.dstack((self.RDLayers[layer].geometry[area].exterior.coords.xy[X],self.RDLayers[layer].geometry[area].exterior.coords.xy[Y])).tolist()[0]
        
        self.RDLayers_area_vertices_coords_dict[layer]=area_vertices_coords_dict  
        
        
        """ this code generates self.layers_vertices_area_dict[layer]
        * self.RDLayers_vertices_area_dict[layer] is a dict where the key is a vertex index, and the value is its area
        """
            
        self.RDLayers_vertices_area_dict[layer] = dict()
        for area in self.RDLayers[layer].index:
            if(type(self.RDLayers[layer].geometry[area])==shapely.geometry.multipolygon.MultiPolygon):
              continue
            area_vert_list=np.dstack((self.RDLayers[layer].geometry[area].exterior.coords.xy[X],self.RDLayers[layer].geometry[area].exterior.coords.xy[Y])).tolist()[0]
            for vert in area_vert_list:
                vert_ind=vert_dict[tuple(vert)]
                self.RDLayers_vertices_area_dict[layer][vert_ind] = area
        
            
        if build_tree:
             self.RDLayers_balltree[layer] = BallTree(np.deg2rad(np.c_[np.array(self.RDLayers_vert_list[layer])]), metric='haversine')
        return True


     def query_point_in_RDLayer(self, layer, lat, long, radius, unique=False):
        """
        Query all points in a RDLayer which are within a Euclidian distance 'radius' from a query point in 
        lat-long format.
    
        Parameters
        ----------
        layer (int): layer position to be searched, according to self.RDLayers
        lat (float): latitude of the query point
        long (float): longitude of the query point       
        radius (float): maximum Euclidian distance to search points, in kilometers
        unique (bool): whether to return a single point per area (polygon), within the distance radius. The nearest point of each area is returned.
            
        Returns
        -------
            a list of points within maximum distance to the query point, sorted by distance.
        """
        
        try:
            self.RDLayers_balltree[layer]
        except:
            print('BallTree does not exist. building...')
            self.RDLayers_balltree[layer] = BallTree(np.deg2rad(np.c_[np.array(self.RDLayers_vert_list[layer])]), metric='haversine')

        near_point_list=self.RDLayers_balltree[layer].query_radius(convert_deg_query(lat,long),
                                                                   radius/EARTH_RADIUS,return_distance=True,sort_results=True)[0][0]   
            
        if not unique:
            return near_point_list.tolist()
        
        unique_list = []
        unique_areas = set()
        
        for point in near_point_list:
            area = self.RDLayers_vertices_area_dict[layer][point]
            if area not in unique_areas:
                unique_areas.add(area)
                unique_list.append(point)
        
        return unique_list

    
     def load_PBLayer(self, filename, zone_number=23, zone_letter='K', build_tree=True, error_x=0.0, error_y=0.0):
        """
        Load a SHP or KML file made up of places that are georreferenced as single points (positions) in UTM format. Points are converted to lat-long format and stored in self.PBLayers[layer]['latlong'] as tuples.
    
        Parameters
        ----------
        filename : string
            SHP filename. 
        zone_number: int
            Zone number. For SÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â£o Paulo, use 23.
        zone_letter: string
            Zone letter. For SÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â£o Paulo, use "K".
        build_tree: bool
            Whether to build a BallTree of the layer points.
        error_x: float
            Translate the first UTM coordinates if the source is wrong.
        error_y: float
            Translate the second UTM coordinates if the source is wrong.
        
        Returns
        -------
            returns the layer index or -1 if it fails.
        """
        extension_str = filename[-3:]
        if(extension_str.lower()=='kml'):
            gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
            try:
                pointdf = gpd.read_file(filename, driver='KML')
                pointdf['latlong'] = [(p.y,p.x) for p in pointdf.geometry]
            except:
                return -1
        elif(extension_str.lower()=='shp'):
            try:
                pointdf = gpd.read_file(filename, driver='SHP')
                pointdf['latlong'] = [utm_to_latlon(tuple(map(operator.add, geom.coords[0], tuple((error_x,error_y)))), zone_number, zone_letter) for geom in pointdf.geometry]
            except:
                return -1
        
        self.PBLayers.append(pointdf)
        layer=len(self.PBLayers)-1
        if build_tree:
             self.PBLayers_balltree[layer] = BallTree(np.deg2rad(np.c_[np.array(pointdf['latlong'].to_list())]), metric='haversine')
        return layer


     def load_PBLayer_csv_data(self, filename, lat_key_column, lng_key_column, build_tree=True, project_to_corners=False):
        """
        Load a CSV file made up of coordinates data and appends it to a PBLayer.
    
        Parameters
        ----------
        filename (string): CSV filename. 
        lat_key_column (string): Name of the key column of the CSV file that matches the latitude.
        lng_key_column (string): Name of the key column of the CSV file that matches the longitude.
        BUILD_TREE
        PROJECT_TO_CORNERS
        
        Returns
        -------
            returns True if the CSV file is succesfully loaded.
        """
        try:
            pointdf = pd.read_csv(filename)
        except NameError as e:
             print(e)  
             return False
        
        if not lat_key_column in pointdf.columns:
            print('invalid lat_key_column.')
            return False

        if not lng_key_column in pointdf.columns:
            print('invalid lng_key_column.')
            return False


        pointdf = pointdf[pointdf[lat_key_column].notna()]
        
        if(isinstance(pointdf.iloc[0][lat_key_column], str)):
            pointdf['latlong'] = [(float(pointdf.iloc[p][lat_key_column].replace(',','.')),float(pointdf.iloc[p][lng_key_column].replace(',','.'))) for p in range(pointdf.shape[0])]

        else:
            pointdf['latlong'] = [(pointdf.iloc[p][lat_key_column],pointdf.iloc[p][lng_key_column]) for p in range(pointdf.shape[0])]

        #pointdf = gpd.GeoDataFrame(pointdf, geometry=gpd.points_from_xy(pointdf['latlong'].apply(lambda x: x[1]), pointdf['latlong'].apply(lambda x: x[0])))

        self.PBLayers.append(pointdf)
        layer=len(self.PBLayers)-1

        self.PBLayers_csv_lat_key_column[layer] = lat_key_column
        self.PBLayers_csv_lng_key_column[layer] = lng_key_column

        if project_to_corners:
            self.project_PBLayer_to_corners(layer)

        if build_tree:
             self.PBLayers_balltree[layer] = BallTree(np.deg2rad(np.c_[np.array(pointdf['latlong'].to_list())]), metric='haversine')
        return layer
     
     def project_PBLayer_to_corners(self, layer):
        """
        Projects PBLayer points to the nearest city street graph corner. The layer must be loaded first. 
    
        Parameters
        ----------
        layer (int): layer position according to self.PBLayers
        build_tree (bool): builds a BallTree for querying
        
        Returns
        -------
            returns True if the preprocessing succeeds.
        """

        if len(self.PBLayers)<=layer:
            return False
        
        #pointdf["vert"] = pointdf.apply(lambda row: self.query_point_in_city_mesh(row['latlong'][0],row['latlong'][1], True), axis = 1)

        vertices_points = [[]] * len(self.city_vert_dict.values())
        vert_column = []
        for point in self.PBLayers[layer].index:
            vrt = self.query_point_in_city_mesh(self.PBLayers[layer].loc[point]['latlong'][0],self.PBLayers[layer].loc[point]['latlong'][1], True)
            vertices_points[vrt] = vertices_points[vrt] + [point]
            vert_column.append(vrt)
        
        self.PBLayers_corner_projection_dicts[layer] = vertices_points     
        self.PBLayers[layer]["vert"] = vert_column

        return True
    
    
     def query_point_in_PBLayer(self, layer, lat, long, radius, return_nearest_indices=False):
        """
        Query all points in a PBLayer which are within a Euclidean distance 'radius' from a query point in 
        lat-long format.
        
        Parameters
        ----------
        layer (int): layer to be searched, according to its position in self.PBLayers
        lat (float): latitude of the query point
        long (float): longitude of the query point       
        radius (float): maximum Euclidian distance to search points, in kilometers
        return_nearest_indices: bool
            whether to return the nearest point indices (according to PBLayers[layer] indexing) or a list of lat-long tuples
        Returns
        -------PALayers_aggregated_polygon_vertices_indices_dict
            a list of points within maximum distance to the query point, sorted by Euclidean distance.
        """
        
        try:
            self.PBLayers_balltree[layer]
        except:
            print('BallTree does not exist. building...')         
            self.PBLayers_balltree[layer] = BallTree(np.deg2rad(np.c_[np.array(self.PBLayers[layer])]), metric='haversine')
        nearest_inds_list=self.PBLayers_balltree[layer].query_radius(convert_deg_query(lat,long),
                                                                     radius/EARTH_RADIUS,return_distance=True,sort_results=True)[0][0].tolist()
        
        if return_nearest_indices:
            return nearest_inds_list
        else:
            return self.PBLayers[layer]['latlong'][nearest_inds_list].tolist()        
        

     def load_SMLayers_known_measurements(self, measu_point_name, measu_filename):
        """
        Load a dataframe from csv file corresponding to time series of
        measurements taken on a single measurement point tagged 'measu_point_name'.
    
        Parameters
        ----------
        measu_point_name : str
            tag for the name of the measurement measu_point_name.
        
        measu_filename : str
            filename in csv format containing temporal informations about
            measurements from a single measurement point.
        
        Returns
        -------
            returns True if the file is sucessfully loaded.
        """
                                              
        try:
            self.SMLayers_known_measurements[measu_point_name] = pd.read_csv(measu_filename, sep=';')
        except:
            return False
        
        self.SMLayers_known_measurements[measu_point_name] = self.SMLayers_known_measurements[measu_point_name].set_index(self.SMLayers_known_measurements[measu_point_name].columns[0])
        self.SMLayers_known_measurements[measu_point_name].index = pd.to_datetime(self.SMLayers_known_measurements[measu_point_name].index,dayfirst=True)
        return True
    
     def load_SMLayers_temp_agg_funcs(self, info_filename):
       """
       Load a dataframe from csv file corresponding to informations about all
       measurement temporal aggregation functions (e.g. sum, mean, max, min).
   
       Parameters
       ----------        
       info_filename : str
           filename in csv format containing informations about all
           measurement temporal aggregation functions.
       
       Returns
       -------
           returns True if the file is sucessfully loaded.
       """
                                             
       try:
           self.SMLayers_temp_agg_funcs = pd.read_csv(info_filename, sep=';')
       except:
           return False
       
       self.SMLayers_temp_agg_funcs = self.SMLayers_temp_agg_funcs.set_index(self.SMLayers_temp_agg_funcs.columns[0])
       self.SMLayers_temp_agg_funcs = self.SMLayers_temp_agg_funcs.iloc[:,0].to_dict()
       return True 
    
     def load_SMLayers_measurements_points_info(self, info_filename):
        """
        Load a dataframe from csv file corresponding to informations about all
        measurement points (e.g. lat, lng, initial measurement date, final
        measurement date).
    
        Parameters
        ----------        
        info_filename : str
            filename in csv format containing informations about all
            measurement points.
        
        Returns
        -------
            returns True if the file is sucessfully loaded.
        """
                                              
        try:
            self.SMLayers_measurements_points_info = pd.read_csv(info_filename, sep=';')
        except:
            return False
        
        self.SMLayers_measurements_points_info = self.SMLayers_measurements_points_info.set_index(self.SMLayers_measurements_points_info.columns[0])
        return True

     def kernel_gaussian(self, X, Xm, sigma):
        """
        Apply gaussian kernel centered at measurement points Xm with
        standard deviation sigma and assess kernel coefficients at X.
        
        Parameters
        ----------        
        X : 2-array
            Nx2 array indicating lat-lng location of each point in which
            estimation will be made (N = number of points of interest).
            
        Xm : 2-array
            Mx2 array indicating lat-lng location of each measurement point on
            which estimaton will be based (M = number of points of known
            measurements).
            
        sigma : float
            Standard deviation of the gaussian kernel employed.
        
        Returns
        -------
            returns a NxM array indicating the kernel coefficient for each
            point of interest versus each measurement point.
        """
        
        # Compute all distances between X and Xm
        M = Xm.shape[0]
        D = [ ox.distance.great_circle_vec(X[:,0],X[:,1],Xm[m,0],Xm[m,1]) for m in range(M) ]
        D = np.array(D).T
        
        # Compute kernel coefficientes
        K = np.exp(-(D/sigma)**2/2)
        return K

     def kernel_estimation(self, X, Xm, Ym, kernel_function='auto'):
        """
        Auxiliar function for sparse_estimation(). Uses kernel coefficients to
        estimate outputs from X based on known Xm to Ym mapping and
        kernel_function.
    
        Parameters
        ----------        
        X : 2-array
            Nx2 array indicating lat-lng location of each point in which
            estimation will be made (N = number of points of interest).
            
        Xm : 2-array
            Mx2 array indicating lat-lng location of each measurement point on
            which estimaton will be based (M = number of points of known
            measurements).
            
        Ym : 2-array
            TxM array indicating time series measured in each measurement
            point (T = number of time slots).
        
        kernel_function : function or 'auto'
            function expecting as input X (Nx2) and Xm (Mx2) as above and
            returning a NxM array of kernel coefficients for each point of
            interest versus each measurement point. By default ('auto'),
            computes coefficients using a gaussian kernel with standard
            deviation equals half the average distance between two different
            measurement point.
        
        Returns
        -------
            returns a NxT array indicating estimation for each point of
            interest and time slot.
        """
        
        M = Xm.shape[0]
        N = X.shape[0]
        if M==1:
            Y = np.tile( Ym.T , (N,1) )
            return Y
        
        if kernel_function=='auto':
            sigma = .5*np.mean([ ox.distance.great_circle_vec(Xm[i,0],Xm[i,1],Xm[j,0],Xm[j,1]) for i in range(M) for j in range(1+i,M)])
            K = self.kernel_gaussian(X,Xm,sigma)
        else:
            K = kernel_function(X,Xm)
        #Knorm = K.div(K.sum(axis=1), axis=0) #for pd dataframe
        Knorm = K / K.sum(axis=1).reshape((N,1)) #for np array
        Y = Knorm @ Ym.T
        return Y
        
         
     def sparse_ts_estimation(self, measu_name, initial_date, final_date, query_points='all', lat_column='lat', lng_column='lng', kernel_function='auto'):
        """
        Estimates time series of measurement measu_name from initial_date to
        final_date in points query_points based on time series measured at the
        measurement points.
    
        Parameters
        ----------        
        measu_name : str
            column name from self.SMLayers_known_measurements which estimation
            will be made.
            
        initial_date : str
            date (in format dd-mm-yyyy) from which estimation will be made.
            
        final_date : str
            date (in format dd-mm-yyyy) until which estimation will be made.
            
        query_points : pandas dataframe or 'all'
            points in which estimation will be made. If a pandas dataframe is
            passed, columns lat_column and lng_column are expected to show
            points latitude and longitude respectively. If 'all' is passed
            (default), estimation will be made in all graph vertices.
        
        lat_column : str
            name of the column with latitude information.
            
        lng_column : str
            name of the column with longitude information.
        
        kernel_function : function or 'auto'
            function expecting as input X (Nx2) and Xm (Mx2) as above and
            returning a NxM array of kernel coefficients for each point of
            interest versus each measurement point. By default ('auto'),
            computes coefficients using a gaussian kernel with standard
            deviation equals half the average distance between two different
            measurement point.
        
        Returns
        -------
            returns a dataframe of shape (num_query_points, num_time_steps),
            or False if parameters passed were not possible to make
            estimations. If a dataframe is passed in query_points, dataframe
            returned is input dataframe with added columns.
        """
        
        initial_date = pd.to_datetime(initial_date,dayfirst=True)
        final_date = pd.to_datetime(final_date,dayfirst=True)
        
        if final_date <= initial_date:
            print('Initial date must be older than final date')
            return False
        
        if type(query_points)==str:
            if query_points=='all':
                query_points = [ [self.city_vert_nxind_to_ind_dict[n],
                                  self.city_street_graph.nodes[n]['y'],
                                  self.city_street_graph.nodes[n]['x']] for n in self.city_street_graph.nodes() ]
                query_points = pd.DataFrame(query_points, columns=['corner',lat_column,lng_column])
                query_points = query_points.set_index('corner')
            else:
                print('query_points must be dataframe or str \'all\'')
        
        slice_measu = pd.DataFrame(columns=self.SMLayers_measurements_points_info.index)
        for measurement_point in self.SMLayers_measurements_points_info.index:
            slice_measu[measurement_point] = self.SMLayers_known_measurements[measurement_point][measu_name]
        slice_measu = slice_measu.loc[(pd.to_datetime(slice_measu.index)>=initial_date) & (pd.to_datetime(slice_measu.index)<final_date)]
        
        measu_estimation = self.kernel_estimation(np.array([query_points[lat_column],query_points[lng_column]]).T,
                                                  np.array([self.SMLayers_measurements_points_info[lat_column],self.SMLayers_measurements_points_info[lng_column]]).T,
                                                  np.array(slice_measu),
                                                  kernel_function)
        
        measu_estimation = pd.DataFrame(measu_estimation,columns=slice_measu.index)
        measu_estimation[query_points.index.name] = query_points.index
        measu_estimation = measu_estimation.set_index(query_points.index.name)
        
        return measu_estimation
    
     def sparse_agg_estimation(self, measu_name, temp_agg_func, initial_date, final_date, query_points='all', lat_column='lat', lng_column='lng', kernel_function='auto'):
        """
        Aggregate measurements measu_name from initial_date to final_date
        using temp_agg_func to obtain a single number for each measurement
        point. Then estimates the measurement in query_points based on these
        numbers.
    
        Parameters
        ----------        
        measu_names : str
            column name from self.SMLayers_known_measurements which estimation
            will be made.
        
        temp_agg_func : str
            string with function to be used in the temporal aggregation of
            measuarements for each measurement point (can be sum, mean, min or
            max).
        
        initial_date : str
            date (in format dd-mm-yyyy) from which estimation will be made.
            
        final_date : str
            date (in format dd-mm-yyyy) until which estimation will be made.
            
        query_points : pandas dataframe or 'all'
            points in which estimation will be made. If a pandas dataframe is
            passed, columns lat_column and lng_column are expected to show
            points latitude and longitude respectively. If 'all' is passed
            (default), estimation will be made in all graph vertices.
        
        lat_column : str
            name of the column with latitude information.
            
        lng_column : str
            name of the column with longitude information.
        
        kernel_function : function or 'auto'
            function expecting as input X (Nx2) and Xm (Mx2) as above and
            returning a NxM array of kernel coefficients for each point of
            interest versus each measurement point. By default ('auto'),
            computes coefficients using a gaussian kernel with standard
            deviation equals half the average distance between two different
            measurement point.
        
        Returns
        -------
            returns a dataframe of shape (num_query_points, num_time_steps),
            or False if parameters passed were not possible to make estimations.
            If a dataframe is passed in query_points, dataframe returned is input
            dataframe with added columns.
        """
        
        if temp_agg_func not in ['sum', 'mean', 'min', 'max']:
            print('Possible temp_agg_func strings: \'sum\', \'mean\', \'min\' or \'max\'.')
            return False
        
        initial_date = pd.to_datetime(initial_date,dayfirst=True)
        final_date = pd.to_datetime(final_date,dayfirst=True)
        
        if final_date <= initial_date:
            print('Initial date must be older than final date')
            return False
        
        if type(query_points)==str:
            if query_points=='all':
                query_points = [ [self.city_vert_nxind_to_ind_dict[n],
                                  self.city_street_graph.nodes[n]['y'],
                                  self.city_street_graph.nodes[n]['x']] for n in self.city_street_graph.nodes() ]
                query_points = pd.DataFrame(query_points, columns=['corner',lat_column,lng_column])
                query_points = query_points.set_index('corner')
            else:
                print('query_points must be dataframe or str \'all\'')
            
        available_measurement_points = self.SMLayers_measurements_points_info

        agg_measu = pd.DataFrame(columns=[measu_name])
        for measurement_point in self.SMLayers_measurements_points_info.index:
            ts = self.SMLayers_known_measurements[measurement_point][measu_name].loc[
                                       (pd.to_datetime(self.SMLayers_known_measurements[measurement_point].index)>=initial_date)&
                                       (pd.to_datetime(self.SMLayers_known_measurements[measurement_point].index)< final_date)]
            if temp_agg_func=='sum':
                agg_measu.loc[measurement_point,:] = ts.sum()
            elif temp_agg_func=='mean':
                agg_measu.loc[measurement_point,:] = ts.mean()
            elif temp_agg_func=='min':
                agg_measu.loc[measurement_point,:] = ts.min()
            else:
                agg_measu.loc[measurement_point,:] = ts.max()
        
        measu_estimation = self.kernel_estimation(np.array([query_points[lat_column],query_points[lng_column]]).T,
                                                  np.array([available_measurement_points[lat_column],available_measurement_points[lng_column]]).T,
                                                  np.array(agg_measu).T,
                                                  kernel_function)
        measu_estimation = pd.DataFrame(measu_estimation,columns=[measu_name])
        measu_estimation[query_points.index.name] = query_points.index
        measu_estimation = measu_estimation.set_index(query_points.index.name)
        
        return measu_estimation
 

     def compute_layer_sparse(self, initial_date, final_date, kernel_function='auto'):
       """
       Aggregate measurements from all sparse variable from initial_date to
       final_date using temp_agg_funcs given by self.SMLayers_temp_agg_funcs
       to obtain a single number for each measurement point. Then estimates
       the measurement in all graph vertices based on these numbers.
   
       Parameters
       ----------                      
       initial_date : str
           date (in format dd-mm-yyyy) from which estimation will be made.
           
       final_date : str
           date (in format dd-mm-yyyy) until which estimation will be made.
       
       kernel_function : function or 'auto'
           function expecting as input X (Nx2) and Xm (Mx2) as above and
           returning a NxM array of kernel coefficients for each point of
           interest versus each measurement point. By default ('auto'),
           computes coefficients using a gaussian kernel with standard
           deviation equals half the average distance between two different
           measurement point.
       
       Returns
       -------
           returns True if the computation if sucessfull (this method modifies
           the property self.SMLayers_estimated_measurements).
       """
       
       measurement_point = self.SMLayers_measurements_points_info.index[0]
       measu_names = list(self.SMLayers_known_measurements[measurement_point].columns)
       for measurement_point in self.SMLayers_measurements_points_info.index:
           if list(self.SMLayers_known_measurements[measurement_point].columns) != measu_names:
               print('All measurement variables must have the same name.')
               return False
       
       temp_agg_funcs = self.SMLayers_temp_agg_funcs.values()
       for temp_agg_func in temp_agg_funcs:
           if temp_agg_func not in ['sum', 'mean', 'min', 'max']:
               print('Possible temp_agg_func strings: \'sum\', \'mean\', \'min\' or \'max\'.')
               return False
           
       if len(measu_names) != len(temp_agg_funcs):
           print('The amount of aggregation functions must match the number of measurements variables.')
           return False
       
       initial_date = pd.to_datetime(initial_date,dayfirst=True)
       final_date = pd.to_datetime(final_date,dayfirst=True)
       
       query_points = [ [self.city_vert_nxind_to_ind_dict[n],
                         self.city_street_graph.nodes[n]['y'],
                         self.city_street_graph.nodes[n]['x']] for n in self.city_street_graph.nodes() ]
       query_points = pd.DataFrame(query_points, columns=['corner','lat','lng'])
       query_points = query_points.set_index('corner')
    
       measu_estimation = pd.DataFrame(columns=measu_names)
       for measu_name in measu_names:
           measu_estimation[measu_name] = self.sparse_agg_estimation(measu_name, self.SMLayers_temp_agg_funcs[measu_name], initial_date, final_date, query_points, kernel_function=kernel_function)
       
       self.SMLayers_estimated_measurements = measu_estimation
       return True

     def compute_feature_vector(self, vert_ind):
        """
        Create your own method to compute a feature vector given a vertex index.
        
        Parameters
        ----------
        vert_ind (int): vertex index from city mesh
        Returns
        -------
            A list with a feature vector from vert_ind 
        """     

        return []

     def compute_all_feature_vectors(self):
        """
        Compute the feature vectors of all vertices of the city street graph, and assign them to self.feature_vecs
        
        Parameters
        ----------
        
        Returns
        -------
            
        """
        self.feature_vecs = [self.compute_feature_vector(vert_ind) for vert_ind in range(len(self.city_vert_list))]


def convert_deg_query(lat, long):
    return np.array([np.deg2rad(lat),np.deg2rad(long)]).reshape(1,-1)

def load_preprocessed_CityHub(filename):
    """
        Load a preprocessed CityHub file
        
        Parameters
        ----------
        filename (string): pickle file
        Returns
        -------
        Return a CityHub instance
        """ 

    with open(filename, 'rb') as f:
        return pickle.load(f)
    
def utm_to_latlon(coords, zone_number, zone_letter):
    easting = coords[0]
    northing = coords[1]
    return utm.to_latlon(easting, northing, zone_number, zone_letter)  
