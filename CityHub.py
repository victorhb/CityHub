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
         city_street_graph (Graph): A networkX undirected graph.
         
         RDLayers (list): The list of Regional Domain (RD) Layers.
         PBLayers (list): The list of Point-based (PB) Layers.
         PALayers_mesh (list): The list of the meshes of the Polygon-Aggregated (PA) layers.
         PALayers_keycolumns (dict): Dictionary of the key column of each PA layer.
         
         
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
        self.refined_city_vert_set = set()
        
        
        ##layers' data
        
        # Regional Domain (RD) Layer 
        self.RDLayers = []
        self.RDLayers_balltree = dict()       
        self.RDLayers_vert_list = dict() 
        
        # Point-Based (PB) Layer         
        self.PBLayers = []
        self.PBLayers_balltree = dict()
        
        # Polygon-Aggregated (PA) Layer
        self.PALayers_mesh = []
        self.PALayers_keycolumns = dict()
        
        self.PALayers_polygon_vertices_dict = dict()
        self.PALayers_aggregated_polygon_vert_list = dict()
        self.PALayers_aggregated_polygon_vert_dict = dict()
        self.PALayers_aggregated_polygon_vertices_indices_dict = dict()
        self.PALayers_aggregated_polygon_vert_setcens_list = dict()
        self.PALayers_refined_aggregated_polygon_vert_list = dict()
        self.PALayers_refined_aggregated_polygon_vert_dict = dict()
        self.PALayers_refined_aggregated_polygon_vert_correspondence_dict = dict()
        self.PALayers_aggregated_polygon_tree = dict()
        
        
        
        
        self.layers_area_vertices_indices_dict = dict()
        self.layers_area_vertices_coords_dict = dict()
        self.layers_vertices_area_dict = dict()

        
        #feature vectors data
        self.feature_vecs = []
        
        #visualization variables
        self.nearest_vertex_index = 0
        self.nearest_subnormal_list = []
     
        #variables from sparse measurements
        self.station_measu = dict()
        self.measu_temp_agg_funcs = dict()
        self.station_info = pd.DataFrame()

        if not self.load_city_street_graph(filename):
            print('Cannot open '+filename)
            return
        self.preprocess_city_mesh(True,EDGE_LENGTH_THRESHOLD)
            
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
    
    
     def load_PALayer_mesh(self, filename, key_column = 'Name', EDGE_LENGTH_THRESHOLD=sys.float_info.max):
        """
        Load a Polygon-Aggregated Layer mesh from a file for polygon-aggregated data purposes, 
        such as census sectors data. KML and SHP formats are accepted.
        In case of success, the GeoDataFrame city_mesh_df will be created.
    
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
            returns True if the PA Layer mesh is sucessfully loaded from file.
        """
        
        extension_str = filename[-3:]
        if(extension_str.lower()=='kml'):
            gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
            try:
                self.PALayers_mesh.append(gpd.read_file(filename, driver='KML'))
            except:
                return False
        elif(extension_str.lower()=='shp'):
            try:
                self.PALayers_mesh.append(aggregated_polygon_mesh=gpd.read_file(filename, driver='SHP'))
            except:
                return False
        self.PALayers_keycolumns[len(self.PALayers_mesh)-1]=key_column
        self.preprocess_PALayer_mesh(len(self.PALayers_mesh)-1,True,True,EDGE_LENGTH_THRESHOLD)
        return True

    
     def save_preprocessed_CityHub(self, filename):
        """
        Save preprocessed data structures of a city street graph to a pickle file .
    
        Parameters
        ----------
        filename : string
            Pickle filename to be written.
        """
        pickle.dump(self, open(filename, 'wb'))
        
        
     def preprocess_city_mesh(self, build_tree=True, EDGE_LENGTH_THRESHOLD=sys.float_info.max):
        """
        Generate data structures to quickly retrieve relevant information. A valid city_mesh_df is required.
    
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

            for u,v,a in self.city_street_graph.edges(data=True):
                vertA = self.city_vert_list[self.city_vert_nxind_to_ind_dict[u]]
                vertB = self.city_vert_list[self.city_vert_nxind_to_ind_dict[v]]
                great_circle_dist = ox.distance.great_circle_vec(vertA[0],vertA[1],vertB[0],vertB[1])
                latlong_dist = ox.distance.euclidean_dist_vec(vertA[0],vertA[1],vertB[0],vertB[1])
          #      self.city_street_graph[u][v]['weight'] = great_circle_dist
                
                self.refined_city_vert_set.add(vertA)
                refined_city_vert_coords_correspondence_dict[vertA] = [self.city_vert_nxind_to_ind_dict[u],self.city_vert_nxind_to_ind_dict[v]]
                if(great_circle_dist>EDGE_LENGTH_THRESHOLD):
                    latlong_threshold = EDGE_LENGTH_THRESHOLD * latlong_dist / great_circle_dist
                    line = shapely.geometry.LineString([(vertA[0],vertA[1]), (vertB[0],vertB[1])])
                    num_vert = max(math.ceil(line.length / latlong_threshold), 1)
                
                    for n in range(1,num_vert-1):
                        interpolated_coords = line.interpolate(n / num_vert, normalized=True).coords[0]
                        self.refined_city_vert_set.add(interpolated_coords)
                        refined_city_vert_coords_correspondence_dict[interpolated_coords] = [self.city_vert_nxind_to_ind_dict[u],self.city_vert_nxind_to_ind_dict[v]]
                self.refined_city_vert_set.add(vertB)
                refined_city_vert_coords_correspondence_dict[vertB] = [self.city_vert_nxind_to_ind_dict[u],self.city_vert_nxind_to_ind_dict[v]]
                
            self.refined_city_vert_list = list(self.refined_city_vert_set)
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
            self.build_city_tree()
        
        return True

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
        for ap in self.self.PALayers_mesh[layer].index:
            if(type(self.self.PALayers_mesh[layer].geometry[ap])==shapely.geometry.multipolygon.MultiPolygon):
                continue
            polygon_vertices_dict[self.PALayers_mesh[layer][self.key_column][ap]] = np.dstack((self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[X],self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[Y])).tolist()[0]
        PALayers_polygon_vertices_dict[layer]=polygon_vertices_dict
        
    
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
            self.PALayers_aggregated_polygon_vertices_indices_dict[layer][self.PALayers_mesh[layer][self.key_column][ap]] = polygon_vert_ind_list     

        """ this code generates PALayers_aggregated_polygon_vert_setcens_list[layer] from PALayers_aggregated_polygon_vert_dict[layer] and the original kml dataframe (self.PALayers_mesh[layer])
        * self.PALayers_aggregated_polygon_vert_setcens_list[layer] is a list comprised of sets of 1-ring polygons (represented with its code) of each vertex index."""

        self.PALayers_aggregated_polygon_vert_setcens_list[layer] = [set() for _ in range(len(self.PALayers_aggregated_polygon_vert_list[layer]))] 
        for ap in self.PALayers_mesh[layer].index:
            if(type(self.PALayers_mesh[layer].geometry[ap])==shapely.geometry.multipolygon.MultiPolygon):
                continue
            polygon_vert_list=np.dstack((self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[X],self.PALayers_mesh[layer].geometry[ap].exterior.coords.xy[Y])).tolist()[0]
            for vert in polygon_vert_list:  
                self.PALayers_aggregated_polygon_vert_setcens_list[layer][self.PALayers_aggregated_polygon_vert_dict[layer][tuple(vert)]].add(self.PALayers_mesh[layer][self.key_column][ap])
                 
                    
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
            for k in self.PALayers_refined_aggregated_polygon_vert_list[layer].keys():                
                ind_refined = self.PALayers_refined_aggregated_polygon_vert_list[layer][k]
                ind_original = refined_vert_coords_correspondence_dict[k] 
                self.PALayers_refined_aggregated_polygon_vert_correspondence_dict[layer][ind_refined] = ind_original
        else:
            self.PALayers_refined_aggregated_polygon_vert_list[layer] = self.PALayers_aggregated_polygon_vert_list[layer].copy()
            self.PALayers_refined_aggregated_polygon_vert_dict[layer] = self.PALayers_aggregated_polygon_vert_dict[layer].copy()
            self.PALayers_refined_aggregated_polygon_vert_correspondence_dict[layer] = self.PALayers_refined_aggregated_polygon_vert_list[layer]
            
            
        if build_tree:
             self.PALayers_aggregated_polygon_tree[layer] = BallTree(np.deg2rad(np.array(self.PALayers_refined_aggregated_polygon_vert_list[layer])), metric='haversine')
        return True
    
    
     def load_polygon_csv_data(self, filename, polygon_key_column):
        """
        Load a CSV file made up of polygon-aggregated data and generate datadf dataframe.
    
        Parameters
        ----------
        filename : string
            CSV filename. 
        polygon_key_column: string
            Name of the key column of the CSV file that matches the polygon's key column of the city mesh dataframe (city_mesh_df)
        
        Returns
        -------
            returns True if the CSV file is succesfully loaded.
        """
        try:
            self.datadf=pd.read_csv(filename)
            self.csv_polygon_key_column = polygon_key_column
        except NameError as e:
             print(e)  
             return False
        
        if polygon_key_column in self.datadf.columns:
            self.csv_polygon_key_column = polygon_key_column
        else:
            print('invalid polygon_key_column.')
            
        return True    
    
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
            returns an integer representing the layer dataframe position in self.RDLayers if the mesh is sucessfully 
            loaded from file, or -1 otherwise.
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
    
     def load_layer_points(self, filename, zone_number=23, zone_letter='K', build_tree=True, error_x=0.0, error_y=0.0):
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
            returns True if the SHP file is succesfully loaded.
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
        return True


     def load_point_csv_data(self, filename, lat_key_column, lng_key_column, build_tree=True):
        """
        Load a CSV file made up of coordinates data and appends it to points layer.
    
        Parameters
        ----------
        filename : string
            CSV filename. 
        lat_key_column: string
            Name of the key column of the CSV file that matches the latitude
        
        Returns
        -------
            returns True if the CSV file is succesfully loaded.
        """
        try:
            pointdf = pd.read_csv(filename)
            self.csv_lat_key_column = lat_key_column
            self.csv_lng_key_column = lng_key_column
        except NameError as e:
             print(e)  
             return False
        
        if lat_key_column in pointdf.columns:
            self.csv_lat_key_column = lat_key_column
        else:
            print('invalid lat_key_column.')
            return False

        if lng_key_column in pointdf.columns:
            self.csv_lng_key_column = lng_key_column
        else:
            print('invalid lng_key_column.')
            return False

        pointdf = pointdf[pointdf[lat_key_column].notna()]
        
        if(isinstance(pointdf.iloc[0][lat_key_column], str)):
            pointdf['latlong'] = [(float(pointdf.iloc[p][lat_key_column].replace(',','.')),float(pointdf.iloc[p][lng_key_column].replace(',','.'))) for p in range(pointdf.shape[0])]
        else:
            pointdf['latlong'] = [(pointdf.iloc[p][lat_key_column],pointdf.iloc[p][lng_key_column]) for p in range(pointdf.shape[0])]


        self.PBLayers.append(pointdf)
        layer=len(self.PBLayers)-1
        if build_tree:
             self.PBLayers_balltree[layer] = BallTree(np.deg2rad(np.c_[np.array(pointdf['latlong'].to_list())]), metric='haversine')
        return True
     
     def preprocess_PBLayer(self, layer, swap_coordinates = False):
        """
        Generate data structures to quickly retrieve relevant information from PBLayer. The layer must be loaded first.
    
        Parameters
        ----------
        layer (int): layer position according to self.PBLayers
        build_tree (bool): builds a BallTree for querying
        swap_coordinates (bool): useful when latitude and longitude coordinates are given in the wrong order.
        
        Returns
        -------
            returns True if the preprocessing succeeds.
        """

        if len(self.PBLayers)<=layer:
            return False
        
        X=0
        Y=1
        if swap_coordinates:
            X=1
            Y=0

        vertices_points = [[]] * len(self.city_vert_dict.values())
        for point in self.PBLayers[layer].index:
            vrt = self.query_point_in_city_mesh(self.PBLayers[layer].loc[point]['latlong'][X],self.PBLayers[layer].loc[point]['latlong'][Y], True)
            vertices_points[vrt] = vertices_points[vrt] + [point]

        if hasattr(self, 'layer_points_vertices'):
            self.layer_points_vertices.update({layer: vertices_points})
        else:
            self.layer_points_vertices = {layer: vertices_points}

        return True


     def load_crimes_layer(self, filename, lat_key_column, lng_key_column):
        """
        Load a CSV file made up of coordinates data of crimes.
    
        Parameters
        ----------
        filename : string
            CSV filename. 
        lat_key_column: string
            Name of the key column of the CSV file that matches the latitude
        
        Returns
        -------
            returns True if the CSV file is succesfully loaded.
        """
        try:
            pointdf = pd.read_csv(filename, index_col = 0)
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

        pointdf["vert"] = pointdf.apply(lambda row: self.query_point_in_city_mesh(row['latlong'][0],row['latlong'][1], True), axis = 1)

        pointdf = pointdf.drop_duplicates()

        self.layer_crimes = pointdf

        #if build_tree:
        #     self.layer_crimes_tree = BallTree(np.deg2rad(np.c_[np.array(pointdf['latlong'].to_list())]), metric='haversine')
        return True  
     
     def compute_layer_crimes(self, initial_date = None, final_date = None):
        """
        Computes crimes per day period for each vertice

        Parameters
        ----------
        build_tree : bool
            builds a BallTree for querying
        swap_coordinates: bool
            useful when latitude and longitude coordinates are given in the wrong order.

        Returns
        -------
            returns True if the preprocessing succeeds.
        """

        tmp_layer_crimes = self.layer_crimes.copy()
        if((not initial_date is None) or (not final_date is None)):
            tmp_layer_crimes["DATAOCORRENCIA"] = pd.to_datetime(tmp_layer_crimes["DATAOCORRENCIA"], format = "%d/%m/%Y", errors = 'coerce')
            tmp_layer_crimes = tmp_layer_crimes[tmp_layer_crimes['DATAOCORRENCIA'].notna()]
            if(not initial_date is None):
                tmp_layer_crimes = tmp_layer_crimes[tmp_layer_crimes["DATAOCORRENCIA"] >= initial_date]
            if(not final_date is None):
                tmp_layer_crimes = tmp_layer_crimes[tmp_layer_crimes["DATAOCORRENCIA"] <= final_date]

        aux = pd.DataFrame(tmp_layer_crimes[["vert","OBJETO"]].value_counts().keys().tolist(), columns = ["vert","OBJETO"])

        aux["n"] = tmp_layer_crimes[["vert","OBJETO"]].value_counts().values
        aux = aux.pivot(index='vert', columns='OBJETO')
        aux.columns = [ x[1] for x in aux.columns ]

        aux["TODOS"] = aux.sum(axis=1)

        crimes_columns = list(tmp_layer_crimes['OBJETO'].unique())
        crimes_count = [[0]*(len(crimes_columns)+1) for _ in range(len(self.city_vert_dict.values()))]

        crimes_count = pd.DataFrame(crimes_count,columns = aux.columns) + aux
        crimes_count = crimes_count.fillna(0)
        
           
        crimes_count.columns = ['CRIME_' + name.replace(' ','_') for name in list(crimes_count.columns)]

        #self.layer_crimes_vertices = vertices_points
        self.layer_crimes_count = crimes_count
        self.layer_crimes_count.rename(columns={"CRIME_CELULAR": "crime_mobile", "CRIME_VEÍCULOS": "crime_vehicle", "CRIME_TODOS": "crime_all"},inplace=True)

        return True

     def load_layer_trips(self, filename):
        """
        Load occurences trip data in CSV format
    
        Parameters
        ----------
        filename : string
            CSV filename. 
        
        Returns
        -------
            returns True if the CSV file is succesfully loaded.
        """
        print("hola")
        try:
            self.tripdf = pd.read_csv(filename, encoding='latin-1')
        except NameError as e:
            print(e)
            return False
        
        self.tripdf = self.tripdf.loc[self.tripdf['fora_capital_start'] == 0]
        self.tripdf = self.tripdf.loc[self.tripdf['fora_capital_dest'] == 0]
        
        return True
    
     def compute_trips_density(self, radius, DENSE_THRESHOLD = 10):
        self.trip_vert_occurrences_starting = [0] * len(self.city_vert_list)
        self.trip_vert_occurrences_dest = [0] * len(self.city_vert_list)
        start_dict=self.tripdf['closest_node_starting'].to_dict()
        end_dict=self.tripdf['closest_node_dest'].to_dict()
        for value in start_dict.values():
            self.trip_vert_occurrences_starting[value] += 1
        for value in end_dict.values():
            self.trip_vert_occurrences_dest[value] += 1       
            
        self.trip_vert_occurrences_starting_density = [0] * len(self.city_vert_list)
        self.trip_vert_occurrences_dest_density = [0] * len(self.city_vert_list)
        
        self.city_tree_coarse = BallTree(np.deg2rad(np.c_[np.array(self.city_vert_list)]), metric='haversine')
        
        
        for i in range(len(self.city_vert_list)):
            coords = self.city_vert_list[i]
            near_point_list=self.city_tree_coarse.query_radius(convert_deg_query(coords[0],coords[1]),radius/EARTH_RADIUS)[0]
            for vert_ind in near_point_list:
                self.trip_vert_occurrences_starting_density[i]+=self.trip_vert_occurrences_starting[vert_ind]
                self.trip_vert_occurrences_dest_density[i]+=self.trip_vert_occurrences_dest[vert_ind]    
        
        self.dense_label = [(self.trip_vert_occurrences_starting_density[i]>DENSE_THRESHOLD)*1 for i in range(len(self.trip_vert_occurrences_starting_density))]

    
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
        this code will ignore areas with multipolygons and generate self.layers_area_vertices_indices_dict[layer] from self.RDLayers[layer]
        * self.layers_vertices_indices_dict[layer] is a dictionary where the key is an area index and the value is a list with the indices of its vertices.
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

        self.layers_area_vertices_indices_dict[layer]=area_vertices_indices_dict  
        
        """ this code will ignore areas with multipolygons and generate self.layers_area_vertices_coords_dict[layer] from self.RDLayers[layer]
        * self.layers_areas_vertices_dict[layer] is a dictionary where the key is an area index, and the value is a list with the coordinates of its vertices
        """
            
        area_vertices_coords_dict = dict()
        for area in self.RDLayers[layer].index:
            if(type(self.RDLayers[layer].geometry[area])==shapely.geometry.multipolygon.MultiPolygon):
              continue
            area_vertices_coords_dict[area] = np.dstack((self.RDLayers[layer].geometry[area].exterior.coords.xy[X],self.RDLayers[layer].geometry[area].exterior.coords.xy[Y])).tolist()[0]
        
        self.layers_area_vertices_coords_dict[layer]=area_vertices_coords_dict  
        
        
        """ this code will ignore areas with multipolygons and generate self.layers_vertices_area_dict[layer]
        * self.layers_vertices_area_dict[layer] is a dict where the key is a vertex index, and the value is its area
        """
            
        self.layers_vertices_area_dict[layer] = dict()
        for area in self.RDLayers[layer].index:
            if(type(self.RDLayers[layer].geometry[area])==shapely.geometry.multipolygon.MultiPolygon):
              continue
            area_vert_list=np.dstack((self.RDLayers[layer].geometry[area].exterior.coords.xy[X],self.RDLayers[layer].geometry[area].exterior.coords.xy[Y])).tolist()[0]
            for vert in area_vert_list:
                vert_ind=vert_dict[tuple(vert)]
                self.layers_vertices_area_dict[layer][vert_ind] = area
        
            
        if build_tree:
             self.RDLayers_balltree[layer] = BallTree(np.deg2rad(np.c_[np.array(self.RDLayers_vert_list[layer])]), metric='haversine')
        return True
    
    
     def build_city_tree(self):
        """
        Build a BallTree with the vertices of the city mesh. Requires preprocessing the city mesh (preprocess_city_mesh). 
        refined_city_vert_list is used to build the tree, but the original vertices will be returned by queries.
    
        Parameters
        ----------
        None
        
        Returns
        -------
            returns True if the preprocessing succeeds.
        """
        #Building the BallTree data structure
        
        print('Building tree...')
        try:
        #    self.tree = KDTree(np.array(self.vert_list))
             self.city_tree = BallTree(np.deg2rad(np.c_[np.array(self.refined_city_vert_list)]), metric='haversine')
        except NameError as e:
        #    print(e)
            return False
        return True
    
    
     def query_point_in_city_mesh(self, lat, long, return_nearest_index=False):
        """
        Query a point in lat-long format, in degrees. Requires a BallTree, which will be built if it does not already exist.
        refined_vert_list is used to build the tree, but the original vertices will be returned by queries.
    
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
        Query a point in lat-long format, in degrees. First, the nearest graph vertex is found. Than, vertices which are less than 'radius' distance, through the graph's edges, are returned.

    
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
                
    #    H = nx.ego_graph(self.city_street_graph, n=self.city_vert_ind_to_nxind_dict[v], radius=query_radius*1000, center=True, undirected=True, distance='length')
    
        ps=nx.single_source_dijkstra_path(self.city_street_graph,self.city_vert_ind_to_nxind_dict[v],cutoff=query_radius*1000.0,weight='length')
        
    
        result_nodes = []
        for k in ps.keys():
            result_nodes.append(self.city_vert_nxind_to_ind_dict[k])    
        return result_nodes  
  
        
     def query_point_in_polygon_mesh(self, layer, lat, long, return_nearest_index=False):
        """
        Query a point in lat-long format, in degrees, in a PALayer. Requires a BallTree, which will be built if it does not already exist.
        refined_aggregated_polygon_vert_list is used to build the tree, but the original vertices will be returned by queries.
    
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
            self.aggregated_polygon_tree
        except:
            print('BallTree does not exist. building...')
            self.aggregated_polygon_tree = BallTree(np.deg2rad(np.array(self.refined_aggregated_polygon_vert_list)), metric='haversine')
            
        [dist,v] =self.aggregated_polygon_tree.query(convert_deg_query(lat,long))
        
        v_org = self.refined_aggregated_polygon_vert_correspondence_dict[v[0][0]]
        self.nearest_vertex_index=v_org
        if return_nearest_index:
            return v_org
        else:
            return self.PALayers_aggregated_polygon_vert_list[layer][v_org]
        
     def query_point_in_RDLayer(self, lat, long, layer, radius, unique=False):
        """
        Query all points in a RDLayer which are within a Euclidian distance 'radius' from a query point in 
        lat-long format. Requires a layer BallTree, which will be built if it does not already exist.
    
        Parameters
        ----------
        lat (float): latitude of the query point
        long (float): longitude of the query point       
        layer (int): layer position to be searched, according to self.RDLayers
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
            area = self.layers_vertices_area_dict[layer][point]
            if area not in unique_areas:
                unique_areas.add(area)
                unique_list.append(point)
        
        return unique_list
    
    
     def query_point_in_PBLayer(self, lat, long, layer, radius, return_nearest_indices=False):
        """
        Query all points in a PBLayer which are within a Euclidean distance 'radius' from a query point in 
        lat-long format. Requires a layer BallTree, which will be built if it does not already exist.
    
        Parameters
        ----------
        lat (float): latitude of the query point
        long (float): longitude of the query point       
        layer (int): layer to be searched, according to its position in self.PBLayers
        radius (float): maximum Euclidian distance to search points, in kilometers
        return_nearest_indices: bool
            whether to return the nearest point indices (according to vert_list indexing) or a list of lat-long tuples
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
        
        
     def polygons_from_vertex(self, vert_ind):
        """
        Compute the polygons in the star of the aggregated polygon vertex 'vert_ind'.
    
        Parameters
        ----------
        vert_ind: int
            vertex index (according to vert_list indexing)
        Returns
        -------
            a list of polygon keys, according to self.key_column
        """        
        self.nearest_subnormal_list = [a for a, b in self.sector_vertices_indices_dict.items() if vert_ind in b]      

        return self.nearest_subnormal_list 
    
    
     def polygon_features_from_vertex(self, vert_ind, feature):
        """
        Retrieve a list of features from the the polygons in the star of the aggregated polygon vertex 'vert_ind'.
    
        Parameters
        ----------
        vert_ind: int
            vertex index (according to vert_list indexing)
        feature: string
            feature name to be retrieved from each neighbour polygon.
        Returns
        -------
            a list of the features 'feature' from each neighbour polygon keys, or an empty list in case the CSV polygon data is not loaded.
        """
        
        if not feature in self.datadf.columns:
            print('invalid feature name '+feature)
            return []
        
        polygon_list=self.polygons_from_vertex(vert_ind)      
        
        feature_list = []
        for polygon in polygon_list:
            dfsearch = self.datadf.loc[self.datadf[self.csv_polygon_key_column]==int(polygon)]
            if dfsearch.shape[0]>0:
                feature_list.append(dfsearch[feature].to_numpy()[0])
        return feature_list
    
     def polygon_features_from_coords(self, layer, lat, long, feature):
        """
        Retrieve a list of features from the the polygons in the star of the nearest vertex to a query point with coordinates (lat,long), among the aggregated polygon vertices.
    
        Parameters
        ----------
        layer (float): layer position according to self.PALayers_mesh.
        lat: float
            latitude of the query point
        long: float
            longitude of the query point
        feature: string
            feature name to be retrieved from each neighbour polygon.
        Returns
        -------
            a list of the features 'feature' from each neighbour polygon keys, or an empty list in case the CSV polygon data is not loaded.
        """        
        
        if not feature in self.datadf.columns:
            print('invalid feature name '+feature)
            return []
        
        vert_ind=self.query_point_in_polygon_mesh(layer,lat, long, True)
        polygon_list=self.polygons_from_vertex(vert_ind)      
        
        feature_list = []
        for polygon in polygon_list:
            feature_val = self.datadf.loc[self.datadf[self.csv_polygon_key_column]==int(polygon)][feature].to_numpy()[0]
            feature_list.append(feature_val)
        return feature_list


     def visualization(self, layer, lat, long, feature, polygon_features):
       """
           Visualization of census sectors close to the closest vertex of a given coordinate.    
        
       Parameters
       ----------
       lat: float
           latitude of the query point
       long: float
           longitude of the query point
       feature: string
           feature name to be retrieved from each neighbour polygon.
       polygon_features: list
        list of generated values from polygon_features_from_coords
       Returns
       -------
           returns a visualization and saves it in HTML
       """        
       print('Generating visualization...')
       
       zip_iterator = zip(self.nearest_subnormal_list,polygon_features)
       sector_feature_dic = dict(zip_iterator)

       m = folium.Map(location=[lat,long], zoom_start=15, tiles='CartoDB positron')
       j=0
       #Real colormap
       colorscale = branca.colormap.linear.YlOrRd_09.scale(self.datadf[feature].min(),self.datadf[feature].max())
       #Colormap considering only the values ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¹ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â¹of the sectors involved
       #colorscale = branca.colormap.linear.YlOrRd_09.scale(sector_feature_dic[min(sector_feature_dic, key=sector_feature_dic.get)],sector_feature_dic[max(sector_feature_dic, key=sector_feature_dic.get)])

       def getcolor(feature):
         return colorscale(sector_feature_dic[self.aggregated_polygon_mesh.Name.loc[ int(feature['id']) ]])
       for i in self.nearest_subnormal_list:
         folium.GeoJson(self.aggregated_polygon_mesh[self.aggregated_polygon_mesh.Name == i].geometry, style_function= lambda feature: {
                 'fillColor': getcolor(feature),
                 'weight': 0.7,
                 'fillOpacity': 0.75,
         }).add_child(folium.Popup(i+' : '+str(polygon_features[j]))).add_to(m)
         polyverts = self.PALayers_polygon_vertices_dict[layer][i]
         j=j+1
         for vert in polyverts:
             folium.CircleMarker([vert[0],vert[1]], radius=3, color='orange').add_to(m)

       folium.CircleMarker([lat,long], radius=7, color='red').add_to(m)                                         #Query point
       folium.CircleMarker(self.PALayers_aggregated_polygon_vert_list[layer][self.nearest_vertex_index], radius=7, color='blue').add_to(m)         #Nearest corner to point

       folium.LayerControl().add_to(m)
       folium.LatLngPopup().add_to(m)
       m
       m.save("map.html")
       return m

     def compute_feature_vector(self, vert_ind):
        
#        feature_vec = pd.Series([vert_ind], index = ['vert_ind'])
        feature_vec = pd.Series([])
#        layer_crimes = 0
        layer_pontodeonibus = 0
        layer_estacaodemetro = 1
        layer_estacaodetrem = 2
        layer_terminaldeonibus = 3
        
        query_coords=self.city_vert_list[vert_ind]


#       crimes
        feature_vec = feature_vec.append(self.layer_crimes_count.loc[vert_ind])

#       ibge
        ibge_list = pd.Series([0.0]*7)
        ibge_list.index = ['Renda_media_por_domicilio', 'Renda_media_responsaveis', 'Responsaveis_sem_renda_taxa', 'Alfabetizados_de_7_a_15_anos', 'menores_de_18_anos_taxa', '18_a_65_anos_taxa', 'maiores_de_65_anos_taxa']
        
        for i in range(len(ibge_list.index)):
            feats = self.polygon_features_from_vertex(vert_ind,  ibge_list.index[i])
            if len(feats)==0:
                return pd.Series([])
            for f in feats:
                if math.isnan(f):
                    return pd.Series([])
            ibge_list[i]=mean(feats)
        feature_vec = feature_vec.append(ibge_list)    
        
#       transporte
        transporte_indices = ['Pontos_de_onibus','Estacao_de_metro', 'Estacao_de_trem', 'Terminal_de_onibus']
        transporte_list = pd.Series([0.0]*len(transporte_indices), index = transporte_indices)
        transporte_list['Pontos_de_onibus'] = len(self.query_point_in_PBLayer(query_coords[0],query_coords[1],layer_pontodeonibus,0.2))
        transporte_list['Estacao_de_metro'] = len(self.query_point_in_PBLayer(query_coords[0],query_coords[1],layer_estacaodemetro,0.2))
        transporte_list['Estacao_de_trem'] = len(self.query_point_in_PBLayer(query_coords[0],query_coords[1],layer_estacaodetrem,0.2))
        transporte_list['Terminal_de_onibus'] = len(self.query_point_in_PBLayer(query_coords[0],query_coords[1],layer_terminaldeonibus,0.2))
        
        feature_vec = feature_vec.append(transporte_list)  
        
#       favelas
        favela_indices = ['Favela_proxima']
        favela_list = pd.Series([0.0]*len(favela_indices), index = favela_indices)
        favela_list['Favela_proxima'] = (len(self.query_point_in_mesh_layer(query_coords[0],query_coords[1],0,0.5,True))>0)*1
        
        feature_vec = feature_vec.append(favela_list)    

#       weather        
        feature_vec = feature_vec.append(self.measu_estimation.loc[vert_ind,:])
        
        return feature_vec

     def compute_all_feature_vectors(self):
        self.feature_vecs = [self.compute_feature_vector(vert_ind) for vert_ind in range(len(self.city_vert_list))]
     
     def load_station_measurements(self, station, measu_filename):
        """
        Load a dataframe from csv file corresponding to time series of
        measurements taken a single place tagged 'station'.
    
        Parameters
        ----------
        station : str
            tag for the name of the measurement station.
        
        measu_filename : str
            filename in csv format containing temporal informations about
            measurements from a single place.
        
        Returns
        -------
            returns True if the file is sucessfully loaded.
        """
                                              
        try:
            self.station_measu[station] = pd.read_csv(measu_filename, sep=';')
        except:
            return False
        
        self.station_measu[station] = self.station_measu[station].set_index('Data')
        self.station_measu[station].index = pd.to_datetime(self.station_measu[station].index,dayfirst=True)
        return True
    
     def load_measu_temp_agg_funcs(self, info_filename):
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
           self.measu_temp_agg_funcs = pd.read_csv(info_filename, sep=';')
       except:
           return False
       
       self.measu_temp_agg_funcs = self.measu_temp_agg_funcs.set_index('Measu_Name')
       self.measu_temp_agg_funcs = self.measu_temp_agg_funcs["Temp_Agg_Func"].to_dict()
       return True 
    
     def load_station_info(self, info_filename):
        """
        Load a dataframe from csv file corresponding to informations about all
        measurement stations (e.g. lat, lng, initial measurement date, final
        measurement date).
    
        Parameters
        ----------        
        info_filename : str
            filename in csv format containing informations about all
            measurement stations.
        
        Returns
        -------
            returns True if the file is sucessfully loaded.
        """
                                              
        try:
            self.station_info = pd.read_csv(info_filename, sep=';')
        except:
            return False
        
        self.station_info['INITIAL DATE'] = pd.to_datetime(self.station_info['INITIAL DATE'],dayfirst=True)
        self.station_info['FINAL DATE'] = pd.to_datetime(self.station_info['FINAL DATE'],dayfirst=True)
        self.station_info = self.station_info.set_index('STATION')
        return True
     
     def kernel_estimation(self, X, Xm, Ym):
        """
        Auxiliar function for sparse_estimation(). Uses kernel coefficients to
        estimate outputs from X based on known Xm to Ym mapping.
    
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
        
        Returns
        -------
            returns a NxT array indicating estimation for each point of
            interest and time slot.
        """
        
        # Compute all distances between elements of X and Xm
        M = Xm.shape[0]
        N = X.shape[0]
        if M==1:
            Y = np.tile( Ym.T , (N,1) )
            return Y
        
        D = [ ox.distance.great_circle_vec(X[:,0],X[:,1],Xm[m,0],Xm[m,1]) for m in range(M) ]
        D = np.array(D).T
        
        # Get kernel coefficientes and make estimations
        sig = .5*np.mean([ ox.distance.great_circle_vec(Xm[i,0],Xm[i,1],Xm[j,0],Xm[j,1]) for i in range(M) for j in range(1+i,M)])
        K = np.exp(-(D/sig)**2/2)
        #Knorm = K.div(K.sum(axis=1), axis=0) #for pd dataframe
        Knorm = K / K.sum(axis=1).reshape((N,1)) #for np array
        Y = Knorm @ Ym.T
        return Y
        
         
     def sparse_ts_estimation(self, measu_name, initial_date, final_date, query_points = 'all', lat_column = 'lat', lng_column = 'lng'):
        """
        Estimates time series of measurement measu_name from initial_date to
        final_date in points query_points based on time series measured in
        stations. Manage initial and final dates to use stations correctly.
    
        Parameters
        ----------        
        measu_name : str
            column name from self.station_measurements which estimation will
            be made.
            
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
        
        Returns
        -------
            returns a dataframe 2-array of shape (num_query_points, num_time_steps),
            or False if parameters passed were not possible to make estimations.
            If a dataframe is passed in query_points, dataframe returned is input
            dataframe with added columns.
        """
        
        initial_date = pd.to_datetime(initial_date,dayfirst=True)
        final_date = pd.to_datetime(final_date,dayfirst=True)
        
        if final_date <= initial_date:
            print('Initial date must be older than final date')
            return False
        
        older_station_initial_date = self.station_info['INITIAL DATE'].min()
        if initial_date < older_station_initial_date:
            print('Older initial date possible: {}'.format(older_station_initial_date))
            return False
        
        sooner_station_final_date = self.station_info['FINAL DATE'].max()
        if final_date > sooner_station_final_date:
            print('Sooner final date possible: {}'.format(sooner_station_final_date))
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
            
        critical_dates = pd.unique(pd.concat( (self.station_info['INITIAL DATE'],self.station_info['FINAL DATE']) ))
        critical_dates = critical_dates[(critical_dates>initial_date)&(critical_dates<final_date)]
        critical_dates = list(critical_dates) + [final_date]
        critical_dates.sort()
        
        measu_estimation = None
        slice_initial_date = initial_date
        for critical_date in critical_dates:
            slice_final_date = critical_date
            available_stations = self.station_info.loc[(self.station_info['INITIAL DATE'] <= slice_initial_date) &
                                                          (self.station_info['FINAL DATE'] >= slice_final_date)]
            
            available_measu = pd.DataFrame(columns=available_stations.index)
            for station in available_stations.index:
                available_measu[station] = self.station_measu[station][measu_name].loc[
                                           (pd.to_datetime(self.station_measu[station].index)>=slice_initial_date)&
                                           (pd.to_datetime(self.station_measu[station].index)< slice_final_date)]
                    
            slice_estimation = self.kernel_estimation(np.array([query_points[lat_column],query_points[lng_column]]).T,
                                                      np.array([available_stations['LAT'],available_stations['LNG']]).T,
                                                      np.array(available_measu))
            
            slice_estimation = pd.DataFrame(slice_estimation,columns=available_measu.index)
            slice_estimation[query_points.index.name] = query_points.index
            slice_estimation = slice_estimation.set_index(query_points.index.name)
            if measu_estimation is None:
                measu_estimation = slice_estimation
                #first_column = measu_estimation.columns[0]
                #measu_estimation = measu_estimation.rename(columns={first_column: measu_name+' '+first_column})
            else:
                measu_estimation = pd.concat( (measu_estimation,slice_estimation) , axis=1)
            slice_initial_date = slice_final_date
        
        return measu_estimation
    
     def sparse_agg_estimation(self, measu_name, temp_agg_func, initial_date, final_date, query_points = 'all', lat_column = 'lat', lng_column = 'lng'):
        """
        Aggregate measurements in measu_name from initial_date to
        final_date using temp_agg_func to obtain a single number for each
        station. Then estimates the measurement in points query_points based
        on these numbers. Only stations with measurements in all time slots
        between initial_date and final_date are used.
    
        Parameters
        ----------        
        measu_names : str
            column name from self.station_measurements which estimation will
            be made.
        
        temp_agg_func : str
            string with function to be used in the temporal aggregation of
            measuarements for each station (can be sum, mean, min or max).
        
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
        
        Returns
        -------
            returns a dataframe 2-array of shape (num_query_points, num_time_steps),
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
        
        older_station_initial_date = self.station_info['INITIAL DATE'].min()
        if initial_date < older_station_initial_date:
            print('Older initial date possible: {}'.format(older_station_initial_date))
            return False
        
        sooner_station_final_date = self.station_info['FINAL DATE'].max()
        if final_date > sooner_station_final_date:
            print('Sooner final date possible: {}'.format(sooner_station_final_date))
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
            
        available_stations = self.station_info.loc[(self.station_info['INITIAL DATE'] <= initial_date) &
                                                          (self.station_info['FINAL DATE'] >= final_date)]

        agg_measu = pd.DataFrame(columns=[measu_name])
        for station in available_stations.index:
            ts = self.station_measu[station][measu_name].loc[
                                       (pd.to_datetime(self.station_measu[station].index)>=initial_date)&
                                       (pd.to_datetime(self.station_measu[station].index)< final_date)]
            if temp_agg_func=='sum':
                agg_measu.loc[station,:] = ts.sum()
            elif temp_agg_func=='mean':
                agg_measu.loc[station,:] = ts.mean()
            elif temp_agg_func=='min':
                agg_measu.loc[station,:] = ts.min()
            else:
                agg_measu.loc[station,:] = ts.max()
        
        measu_estimation = self.kernel_estimation(np.array([query_points[lat_column],query_points[lng_column]]).T,
                                                  np.array([available_stations['LAT'],available_stations['LNG']]).T,
                                                  np.array(agg_measu).T)
        measu_estimation = pd.DataFrame(measu_estimation,columns=[measu_name])
        measu_estimation[query_points.index.name] = query_points.index
        measu_estimation = measu_estimation.set_index(query_points.index.name)
        
        return measu_estimation
    
     def compute_layer_weather(self, initial_date, final_date):
       """
       Aggregate measurements from all sparse variable from initial_date to
       final_date using temp_agg_funcs given by self.temp_agg_funcs to obtain
       a single number for each station. Then estimates the measurement in all
       graph vertices based on these numbers. Only stations with measurements
       in all time slots between initial_date and final_date are used.
   
       Parameters
       ----------                      
       initial_date : str
           date (in format dd-mm-yyyy) from which estimation will be made.
           
       final_date : str
           date (in format dd-mm-yyyy) until which estimation will be made.
       
       Returns
       -------
           returns True if the computation if sucessfull (this method modifies
           the property self.measu_estimation).
       """
       
       station = self.station_info.index[0]
       measu_names = list(self.station_measu[station].columns)
       for station in self.station_info.index:
           if list(self.station_measu[station].columns) != measu_names:
               print('All station measurements must have the same name.')
               return False
       
       temp_agg_funcs = self.measu_temp_agg_funcs.values()
       for temp_agg_func in temp_agg_funcs:
           if temp_agg_func not in ['sum', 'mean', 'min', 'max']:
               print('Possible temp_agg_func strings: \'sum\', \'mean\', \'min\' or \'max\'.')
               return False
           
       if len(measu_names) != len(temp_agg_funcs):
           print('The amount of aggregation functions must match the number of measurements variables.')
           return False
       
       initial_date = pd.to_datetime(initial_date,dayfirst=True)
       final_date = pd.to_datetime(final_date,dayfirst=True)
       
       older_station_initial_date = self.station_info['INITIAL DATE'].min()
       if initial_date < older_station_initial_date:
           print('Older initial date possible: {}'.format(older_station_initial_date))
           return False
       
       sooner_station_final_date = self.station_info['FINAL DATE'].max()
       if final_date > sooner_station_final_date:
           print('Sooner final date possible: {}'.format(sooner_station_final_date))
           return False
       
       query_points = [ [self.city_vert_nxind_to_ind_dict[n],
                         self.city_street_graph.nodes[n]['y'],
                         self.city_street_graph.nodes[n]['x']] for n in self.city_street_graph.nodes() ]
       query_points = pd.DataFrame(query_points, columns=['corner','lat','lng'])
       query_points = query_points.set_index('corner')
    
       measu_estimation = pd.DataFrame(columns=measu_names)
       for measu_name in measu_names:
           measu_estimation[measu_name] = self.sparse_agg_estimation(measu_name, self.measu_temp_agg_funcs[measu_name], initial_date, final_date, query_points)
       
       self.measu_estimation = measu_estimation
       return True

def convert_deg_query(lat, long):
    return np.array([np.deg2rad(lat),np.deg2rad(long)]).reshape(1,-1)

def load_preprocessed_CityHub(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)
    
def utm_to_latlon(coords, zone_number, zone_letter):
    easting = coords[0]
    northing = coords[1]
    return utm.to_latlon(easting, northing, zone_number, zone_letter)  
