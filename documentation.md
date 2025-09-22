# GIS workflow documentation

This documentation highlights the technical pipeline of PV identification of rooftops and leads generation. 

## Singapore

In QGIS, use QUICKOSM to extract all building polygons in Singapore (123 763) as of Aug 2025. Building layer is by OpenStreetMaps, the world's largest community-contributed geospatial tool, functioning as a "mapping equivalent of Wikipedia"

Filter out buildings types, keeping only Commercial & Industrial buildings. Use other additional layers provided my data.gov to filter out layers. *Note* OSM's building tagging or building polygons may not be accurate/up to date, like all community-driven projects

Use Canvas Extent tool in QGIS to draw boundary around Singapore in size 19x19 tiles 

Run Roofpedia's modified ResNet model's inference on created Singapore dataset

Results is a shp file of classified Photovoltaic building polygons 

Execute Overlap Analysis with targeted buildings and PV model

Select features using expression and export only overlap/ non-overlap

Run reverse geo-coding using OneMap API (onemap_rg.py)

Run company matching on results joining on postal code (postal.ipynb)

Check for duplicate companies, OneMap API may tag 2 separate buildings as the same postal code (also postal.ipynb)

Upload list on companies on Apollo.io

Select Location, Job Title & Management level filters

Extract people's contacts information 

Export list & join back company with reverse geo-coded list 

Seb pivots tables to generate aggregated summary page 

Done

## Japan