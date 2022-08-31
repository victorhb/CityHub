# CityHub class
The CityHub class is aimed to provide data structures and algorithms to work with different types of data on a city streets' mesh. In addition to the city mesh, three kinds of data are currently available:

1. Polygon-aggregated city mesh and data: useful for census data, where data is aggregated by census sectors, which are generally defined as polygons. 
2. Poligonal areas layers: useful to describe areas with specific characteristics, such as subnormal agglomerates. Many layers may be loaded separately.
3. Point-of-interest layers: may represent specific types of places that can be described as a single point location, such as bus stations, restaurants, bars, etc. Many layers may be loaded separately.


The file CityHub.py contains the source code of the CityHub class, with the following functionalities:

1. load a city mesh from a Gpickle (networkX graph in pickle format); or from a city string (i.e., 'Sao Paulo, Brazil')
2. load a polygon mesh for polygon-agreggated data purposes from a KML or SHP file
3. generate data structures to quickly retrieve relevant information
4. build BallTrees for fast nearest neighbour search using the haversine metric
8. refine a city mesh or a polygon mesh to improve search accuracy
9. load a CSV file made up of polygon-aggregated data (this functionality requires an associated polygon mesh)
10. load a layer mesh describing poligonal areas with specific characteristics
11. load a point-based layer describing points of interest
12. save a preprocessed CityGraph as a pickle file
13. load a preprocessed CityGraph as a pickle file
14. query points in the city mesh, in the polygon mesh, in layer meshes, and in point-based layers. Query may be performed by radius.
15. retrieve polygons in the star of a polygon mesh vertex
16. retrieve features of polygons in the star of a vertex

Information retrieval may be performed by giving a specific city mesh vertex index, or by providing lat-long coordinates, which will be projected to the nearest city mesh vertex.

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
