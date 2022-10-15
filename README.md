# CityHub library
CityHub is a library to handle multiple urban datasets. It was conceived to represent and organize multiple sources of urban data in the form of data layers that considers four different types of layers. Using edge/node projection algorithms, the integrated data is preprocessed to a common spatial domain: the nodes of a city street graph. The resulting data structures may be easily used to export feature vectors associated to the nodes of a street graph, enabling a multitude of analytical procedures. Each layer may hold multiple elements.

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


The following code first loads São Paulo's street mesh ('data/sp_network.gpickle'), which is refined to assure maximum edge length of 40 meters. Then, a polygon mesh (malha-sp-2010.kml) is loaded and also refined to assure maximum edge length of 40 meters. The associated CSV polygon data file (sectorsindices.csv) is then loaded. Then, the preprocessed CityHub objectred is saved to the file 'data/saopaulo.bin' for faster loading times. 
```
import CityHub
ch = CityHub.CityHub('data/sp_network.gpickle',40.0)
ch.load_polygon_mesh('data/malha-sp-setores-2010.kml','Name',40.0)
ch.load_polygon_csv_data('data/sectorsindices.csv','Cod_setor')
ch.save_preprocessed_CityHub('data/saopaulo.bin')
polygon_features = ch.polygon_features_from_coords(-23.5368789998025, -46.453812,  'Renda_media_por_domicilio')
ch.visualization(-23.5368789998025, -46.453812,  'Renda_media_por_domicilio',polygon_features)
```

In the following, a preprocessed CityHub is loaded, and a specific lat-long point is given as input to retrieve a list of specific features ('Renda_media_por_domicilio') from polygon aggregated data, in its neighborhood. Note that the given lat-long is first projected to the nearest vertex in the polygon mesh vertex, before retrieving the polygons' data. Note that such features could then be assigned to the city mesh vertex that is nearest to lat-long. Finally, the result is visualized in the map.
 
 ```
import CityHub
ch = CityHub.load_preprocessed_CityHub('data/saopaulo.bin')
polygon_features = ch.polygon_features_from_coords(-23.5368789998025, -46.453812,  'Renda_media_por_domicilio')
ch.visualization(-23.5368789998025, -46.453812,  'Renda_media_por_domicilio',polygon_features)
```

Next, in the following code a layer mesh with several distinct polygonal areas (or agglomerates) is loaded, and a query for the layer mesh neighbor vertices of a given lat-long, within a given radius vertices, is performed. The last input parameter (True) specifies that a unique vertex per area is returned.

```
ch.load_layer_mesh('AGSN_2019.shp')
ch.query_point_in_mesh_layer(-23.5368789998025, -46.453812, 0, 1.0, True)
```

Finally, in the following example a preprocessed CityHub file is loaded; four point-based layers holding bus stops positions, subway stations, train stations, and bus terminals from São Paulo are loaded from the shapefiles 'SAD69-96_SHP_pontoonibus.shp', 'SAD69-96_SHP_estacaometro_point.shp', 'SAD69-96_SHP_estacaotrem_point.shp' and 'sad6996_terminal_onibus.shp'; and all bus stops that lie inside a circle of 1km radius centered at query_coords.
Note: the last parameters -45.0, -45.0 are corrections on the input of the georreferenced data files.

```
ch = CityHub.load_preprocessed_CityHub('data/saopaulo.bin')
ch.load_layer_points('SAD69-96_SHP_pontoonibus.shp',23,'K',True,-45.0,-45.0))
ch.load_layer_points('SAD69-96_SHP_estacaometro_point.shp',23,'K',True,-45.0,-45.0))
ch.load_layer_points('SAD69-96_SHP_estacaotrem_point.shp',23,'K',True,-45.0,-45.0))
ch.load_layer_points('sad6996_terminal_onibus.shp',23,'K',True,-45.0,-45.0))
query_coords = (-23.5988, -46.63542)
ch.query_point_in_points_layer(query_coords[0],query_coords[1],0,1)
```
