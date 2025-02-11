# Description of shapefiles

These shapefiles were generated from CT-wide data from 2006 and 2015.  The original data is available from the [UConn CLEAR website](https://clear.uconn.edu/data/index.htm).  They provide [2006 Forest Fragmentation & 1985-2006 Change data](http://clear.uconn.edu/projects/landscape/v2/forestfrag/download2.htm), as well as a [Landscape Fragmentation Tool 2](http://clear.uconn.edu/tools/lft/lft2/index.htm) that we used used to convert [2015 Land Cover data](http://clear.uconn.edu/projects/landscape/download.htm) into a similar format.

## 2015 data

### ct2015_core_wgs84.shp
These are all of the core forest polygons in the state of CT from the 2015 CLEAR dataset in the WGS84 projection.

### ct2015_FEN_core_nad83.shp
These are the core forest polygons in the NAD83 projection for each of the 32 FEN sites from the 2015 CLEAR dataset. They include fragment size and perimeter calculated in NAD83.

### ct2015_FEN_core_wgs84.shp
These are the core forest polygons in the WGS84 projection for each of the 32 FEN sites from the 2015 CLEAR dataset. They include fragment size and perimeter calculated in NAD83.

### ct2015_FEN_1km_nad83.shp
This contains all of the forest within 1km of the borders of each site's core forest polygon from the 2015 CLEAR dataset. It includes area calculated in NAD83 for each of the forest polygons.

### ct2015_FEN_conn_nad83.shp
This contains all of the forest directly connected to each site's core forest polygon from the 2015 CLEAR dataset. It includes area calculated in NAD83 for each of the forest polygons.

## 2006 data

### ct2006_FEN_core.shp
These are the core forest polygons for each of the 32 FEN sites from the 2006 CLEAR dataset. They include area and perimeter calculated in WGS84.

### ct2006_FEN_connected.shp
These are all the forest polygons connected to the FEN core forest areas from the 2006 CLEAR dataset. They include area calculated in WGS84. In some cases, the forests are much larger than the immediate forest around the core fragment. 

### ct2006_FEN_mainforest.shp
These are the forest polygons immediately surrounding the FEN core forest area from the 2006 CLEAR dataset. Areas of > ~10 forest grid cells that are connected by < 3 edge forest cells were removed. For the remaining forest, area is calculated in WGS84

### ct2006_FEN_1km.shp
This includes all forest within 1 km of each FEN site core forest. For many sites, these forest areas overlap. Area for each polygon of forest is included.  This is intended to allow for a characterization of the overall forest within a km of the core.


# State Level Data

These shapefileswere taken from the [Connecticut GIS Data](http://magic.lib.uconn.edu/connecticut_data.html) website, last updated in 2010

### ct2010_state_wgs84.shp
1:100,000 scale state map

### ct2010_counties_wgs84.shp
1:100,000 scale CT counties map
