# CityHub library
CityHub is a library to handle multiple urban datasets. It was conceived to represent and organize multiple sources of urban data in the form of data layers that considers four different types of layers. Using edge/node projection algorithms, the integrated data is preprocessed to a common spatial domain: the nodes of a city street graph. The resulting data structures may be easily used to export feature vectors associated to the nodes of a street graph, enabling a multitude of analytical procedures. Each layer may hold multiple elements. For more details, see our paper [CityHub: A Library for Urban Data Integration](http://sibgrapi.sid.inpe.br/col/sid.inpe.br/sibgrapi/2022/09.14.21.46/doc/sibgrapi2022_cityhub-preprint.pdf)

In addition to city street graph data, four types of data layers are currently available:

1. Point-based (PB) data layers: to represent specific types of places that can be described as a single point location, such as bus stations, restaurants, bars, etc. 
2. Regional domain (RD) data layers: may represent areas with specific characteristics, such as subnormal agglomerates.
3. Polygon-aggregated (PA) data layers: useful for data that is aggregated in polygonal regions, such as census sectors.
4. Sparse data (SD) layers: used to represent data that is sparse in space and requires interpolation, such as temperature and rainfall precipitation.


The file CityHub.py contains the source code of the CityHub library, with the following functionalities:

1. load a city street graph from a Gpickle (networkX graph in pickle format); or from a city string (i.e., 'Sao Paulo, Brazil')
2. load a polygon mesh for polygon-agreggated data purposes from a KML or SHP file
3. generate data structures to quickly retrieve relevant information, including BallTrees for fast nearest neighbour search using the haversine metric.
4. load PA layers from a CSV file made up of polygon-aggregated data and its corresponding polygon mesh from KML or SHP files.
5. load RD layers from KML or SHP files.
6. load PB layers from KML, SHP or CSV files.
7. load SD layers from CSV files and interpolate through Gaussian Kernels or custom kernels.
8. query points in the city street graph using (lat,lon) pairs within a search radius by Euclidian distance or graph distance
9. query points in PA, RD and PB layers using (lat,lon) pairs
10. save a preprocessed CityHub as a pickle file
11. load a preprocessed CityHub as a pickle file


The following code first loads São Paulo's street mesh ('saopaulo_street_graph.gpickle'), which is refined to assure maximum edge length of 40 meters (for edge projection purposes). Then, a PA layer is loaded in two steps: a polygon mesh (saopaulo_malha_census_units-2010.kml) is loaded and also refined to assure maximum edge length of 40 meters; and the associated CSV polygon data file (saopaulo_census_sectors_data-2010.csv) is loaded. Two PB layers, representing bus stops and train stations, and a RD layer representing subnormal agglomerates are also loaded. Finally, a SM layer is loaded with information from 3 different measurement stations. The preprocessed CityHub object is saved to the file 'saopaulo.bin' for faster loading times. 

```
import CityHub
ch = CityHub.CityHub('saopaulo_street_graph.gpickle',40.0)
ch.load_PALayer_mesh('saopaulo_malha_census_units-2010.kml','Name',40.0)
ch.load_PALayers_csv_data(0,'saopaulo_census_sectors_data-2010.csv','Cod_setor')
ch.load_PBLayer('saopaulo_bus_stops.shp',23,'K',True,-45.0,-45.0)
ch.load_PBLayer('saopaulo_train_stations.shp',23,'K',True,-45.0,-45.0)
ch.load_RDLayer('saopaulo_subnormal_agglomerates.shp')
ch.load_SMLayers_known_measurements('Barueri','Barueri.CSV')
ch.load_SMLayers_known_measurements('Interlagos','Interlagos.CSV')
ch.load_SMLayers_known_measurements('Mirante','Mirante.CSV')
ch.load_SMLayers_temp_agg_funcs('Sparse_temp_agg_funcs.CSV')
ch.load_SMLayers_measurements_points_info('Stations_info.CSV')
ch.compute_layer_sparse('01-01-2020','31-12-2020')
ch.save_preprocessed_CityHub('saopaulo.bin')
```

In the following, a preprocessed CityHub is loaded, and a specific lat-long point is given as input to retrieve a list of specific features ('Renda_media_por_domicilio') from polygon aggregated data, in its neighborhood. Note that the given lat-long is first projected to the nearest vertex in the polygon mesh vertex, before retrieving the polygons' data. Note that such features could then be assigned to the city mesh vertex that is nearest to lat-long. 
 ```
import CityHub
ch = CityHub.load_preprocessed_CityHub('data/saopaulo.bin')
polygon_features = ch.polygon_features_from_coords(-23.5368789998025, -46.453812,  'Renda_media_por_domicilio')
```

Following the same code, a specific (lat,lon) pair is queried in the subnormal agglomerates layer, using a 1km search radius. The last input parameter (True) specifies that a unique vertex per area is returned (in case you only need to identify any vertex of nearby subnormal agglomerates).  Then, another (lat,lon) pair is queried in the second PB Layer (train stations), using a 500m search radius, and returning the indices of the nearby vertices (the last input parameter is True).

```
ch.query_point_in_RDLayer(0,-23.5368789998025, -46.453812, 1.0, True)
ch.query_point_in_PBLayer(1,-23.5988, -46.63542,0.5,True)
```

Please cite our paper in your work:

@inproceedings{CityHub2022,
author={Karelia Salinas and Thales Gonçalves and Victor Barella and Thales Vieira and Luis Gustavo Nonato},
title={CityHub: A Library for Urban Data Integration},
booktitle={2022 35nd SIBGRAPI Conference on Graphics, Patterns and Images (SIBGRAPI)},
year={2022},
organization={IEEE}
}
