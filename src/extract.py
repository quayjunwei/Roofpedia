import os
import argparse

from tqdm import tqdm
from PIL import Image
import geopandas as gp
import numpy as np
from shapely.geometry import shape, mapping

from src.tiles import tiles_from_slippy_map
from src.features.building import Roof_features

# flatten 3D geometries (if any) to 2D
def flatten_geometry(g):
    return shape(mapping(g))  # removes Z if present

def mask_to_feature(mask_dir):
    handler = Roof_features()
    tiles = list(tiles_from_slippy_map(mask_dir))

    for tile, path in tqdm(tiles, ascii=True, unit="mask"):
        image = np.array(Image.open(path).convert("P"), dtype=np.uint8)
        mask = (image == 1).astype(np.uint8)
        handler.apply(tile, mask)

    feature = handler.jsonify()
    return feature

def intersection(target_type, city_name, mask_dir):
    print("\n Converting Prediction Masks to GeoJSON Features")
    features = mask_to_feature(mask_dir)
    prediction = gp.GeoDataFrame.from_features(features, crs="EPSG:4326")

    # load building polygons
    city_path = os.path.join('results', '01City', city_name + '.geojson')
    city = gp.read_file(city_path)[['geometry']]

    # clean geometry + update CRS
    city['geometry'] = city['geometry'].apply(flatten_geometry)
    city = city[city.is_valid]
    city['geometry'] = city['geometry'].to_crs('EPSG:4326')
    city['area'] = city['geometry'].area

    # spatial join between prediction & buildings
    intersections = gp.sjoin(city, prediction, how="inner", predicate='intersects')
    intersections = intersections.drop_duplicates(subset=['geometry'])

    # save results
    output_path = os.path.join('results', '04Results', f"{city_name}_{target_type}.geojson")
    intersections.to_file(output_path, driver='GeoJSON')

    print(f"\n Process complete. Footprints with {target_type} roofs saved at:\n{output_path}")
    return intersections

def intersection_from_file(prediction_path, target_type, city_name, mask_dir):
    print("\n Loading Prediction GeoJSON Features")
    prediction = gp.read_file(prediction_path)[['geometry']]

    city_path = os.path.join('results', '01City', city_name + '.geojson')
    city = gp.read_file(city_path)[['geometry']]

    city['geometry'] = city['geometry'].apply(flatten_geometry)
    city = city[city.is_valid]
    city['geometry'] = city['geometry'].to_crs('EPSG:4326')
    city['area'] = city['geometry'].area

    intersections = gp.sjoin(city, prediction, how="inner", predicate='intersects')
    intersections = intersections.drop_duplicates(subset=['geometry'])

    output_path = os.path.join('results', '04Results', f"{city_name}_{target_type}.geojson")
    intersections.to_file(output_path, driver='GeoJSON')

    print(f"\n Process complete. Footprints with {target_type} roofs saved at:\n{output_path}")
    return intersections

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("city", help="City to be predicted, must match dataset name")
    parser.add_argument("type", help="Roof typology: 'Green' or 'Solar'")
    args = parser.parse_args()

    city_name = args.city
    target_type = args.type
    mask_dir = os.path.join("results", "03Masks", target_type, city_name)

    intersection(target_type, city_name, mask_dir)