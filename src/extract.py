import os
import argparse
from collections import defaultdict

from tqdm import tqdm
from PIL import Image
import geopandas as gp
import pandas as pd
import numpy as np

from shapely.ops import unary_union
from shapely.geometry import box
from pyproj import Geod
import mercantile

from src.tiles import tiles_from_slippy_map
from src.features.building import Roof_features

# Global geodesic calculator
GEOD = Geod(ellps="WGS84")


def geodesic_area_m2(geom):
    """
    Compute the geodesic area of a geometry in square metres.

    Uses the WGS84 ellipsoid for accurate area calculation anywhere in the
    world, avoiding the distortion introduced by flat projections.

    Args:
        geom (shapely.geometry.Polygon or MultiPolygon): Input geometry in
              WGS84 (EPSG:4326) coordinates.

    Returns:
        float: Area in square metres (m²). Returns 0.0 for empty or
               unsupported geometry types.
    """
    if geom.is_empty:
        return 0.0

    if geom.geom_type == "Polygon":
        lon, lat = geom.exterior.coords.xy
        area, _ = GEOD.polygon_area_perimeter(lon, lat)
        return abs(area)

    elif geom.geom_type == "MultiPolygon":
        return sum(geodesic_area_m2(p) for p in geom.geoms)

    return 0.0


def mask_to_feature(mask_dir):
    """
    Convert prediction mask tiles to GeoJSON features.

    Reads mask images from a slippy map directory structure, applies each
    mask to the tile handler, and returns the resulting feature collection.

    Args:
        mask_dir (str): Path to the directory containing prediction mask tiles
                        in slippy map format (z/x/y.png).

    Returns:
        tuple:
            - features (list): GeoJSON feature collection of predicted polygons.
            - invalid_summary (dict or None): Summary of invalid/simplified tiles
              from the handler, or None if not supported.
    """
    handler = Roof_features()
    tiles = list(tiles_from_slippy_map(mask_dir))

    for tile, path in tqdm(tiles, ascii=True, unit="mask"):
        image = np.array(Image.open(path).convert("P"), dtype=np.uint8)
        mask = (image == 1).astype(np.uint8)
        handler.apply(tile, mask)

    feature = handler.jsonify()

    # Get invalid tiles summary if available
    invalid_summary = None
    if hasattr(handler, "get_invalid_tiles_summary"):
        invalid_summary = handler.get_invalid_tiles_summary()

    return feature, invalid_summary


def intersection(target_type, city_name, mask_dir):
    """
    Route intersection logic based on roof type.

    Args:
        target_type (str): Roof typology — 'solar' or 'green' (case-insensitive).
        city_name (str): Name of the city, used to locate input files and
                         name output files.
        mask_dir (str): Path to the directory containing prediction mask tiles.

    Returns:
        GeoDataFrame: Result from the appropriate intersection function.

    Raises:
        ValueError: If target_type is not 'solar' or 'green'.
    """
    if target_type.lower() == "solar":
        return intersection_solar(city_name, mask_dir)
    elif target_type.lower() == "green":
        return intersection_green(city_name, mask_dir)
    else:
        raise ValueError(f"Unknown target_type: {target_type}. Use 'Solar' or 'Green'")


def intersection_green(city_name, mask_dir):
    """
    **Original simple logic for green roofs - UNCHANGED
    Process green roof predictions and match with building footprints.

    Converts prediction masks to GeoJSON features, performs a spatial join
    with building polygons, and saves matched buildings as a GeoJSON file.

    Args:
        city_name (str): Name of the city, used to locate input files and
                         name output files.
        mask_dir (str): Path to the directory containing prediction mask tiles.

    Returns:
        GeoDataFrame: Buildings that intersect with at least one green roof
                      prediction, with an additional 'area' column (m²).
    """
    print()
    print("Converting Prediction Masks to GeoJson Features")

    features, _ = mask_to_feature(mask_dir)  # Unpack but ignore invalid_summary
    prediction = gp.GeoDataFrame.from_features(features, crs=4326)

    # loading building polygons
    city_path = f"results/01City/{city_name}.geojson"
    city = gp.read_file(city_path)[["geometry"]]
    city["area"] = city.geometry.map(geodesic_area_m2)

    intersections = gp.sjoin(city, prediction, how="inner", op="intersects")
    intersections = intersections.drop_duplicates(subset=["geometry"])

    output_path = f"results/04Results/{city_name}_Green.geojson"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    intersections.to_file(output_path, driver="GeoJSON")

    print()
    print(f"Process complete, footprints with Green roofs saved at {output_path}")

    return intersections


def intersection_solar(city_name, mask_dir):
    """
    Updated logic for solar roofs with detailed metrics and Excel output.
    Process solar roof predictions and match with building footprints.

    Converts prediction masks to GeoJSON features, performs centroid-based
    spatial matching with building polygons, computes per-building solar
    coverage metrics, and saves results as both GeoJSON and Excel.

    Buildings are classified as solar if predicted PV coverage is >= 20%
    of the roof area. All buildings with any solar detection are included
    in the Excel output regardless of the threshold.

    Args:
        city_name (str): Name of the city, used to locate input files and
                         name output files.
        mask_dir (str): Path to the directory containing prediction mask tiles.

    Returns:
        GeoDataFrame: Buildings with >= 20% solar coverage, with columns:
                      'area', 'solar_area_m2', 'solar_coverage_pct',
                      'available_roof_m2', 'num_predictions'.
                      Returns empty GeoDataFrame if no predictions found
                      or no buildings meet the threshold.
    """
    # predicted features
    features, invalid_summary = mask_to_feature(mask_dir)
    prediction = gp.GeoDataFrame.from_features(features, crs=4326)

    if len(prediction) == 0:
        print("No predictions found")
        return gp.GeoDataFrame()

    # loading building polygons - keep ALL columns
    city_path = f"results/01City/{city_name}.geojson"
    city = gp.read_file(city_path, engine="pyogrio")

    # avoid sjoin column conflicts
    for col in ["index_left", "index_right"]:
        if col in city.columns:
            city = city.rename(columns={col: f"{col}_orig"})

    # calculate building area using ellipsoidal method (worldwide accurate)
    city["area"] = city.geometry.map(geodesic_area_m2)

    # centroid-based matching - only match if prediction center is inside building
    prediction_projected = prediction.to_crs("EPSG:3857")
    prediction["centroid"] = prediction_projected.geometry.centroid.to_crs(
        prediction.crs
    )

    # spatial join using centroids
    pred_centroids = prediction.set_geometry("centroid")
    intersections = gp.sjoin(city, pred_centroids, how="inner", predicate="contains")

    if len(intersections) == 0:
        print("No intersections found")
        return gp.GeoDataFrame()

    # group predictions by building
    building_to_preds = defaultdict(list)
    for idx, row in intersections.iterrows():
        building_to_preds[idx].append(row["index_right"])

    # buildings classified as solar (>=20% coverage)
    solar_rows = []
    all_buildings_with_solar = []
    COVERAGE_THRESHOLD = 20

    # process each building once, combining ALL its predictions
    for building_idx, pred_indices in building_to_preds.items():
        building = city.loc[building_idx].copy()
        building["_building_idx"] = building_idx

        building_geom = building.geometry
        building_area = building["area"]

        # collect and union all solar geometries for this building
        solar_geoms = prediction.loc[pred_indices].geometry.tolist()
        combined_solar = (
            solar_geoms[0] if len(solar_geoms) == 1 else unary_union(solar_geoms)
        )

        # clip to building boundary
        clipped = building_geom.intersection(combined_solar)

        if clipped.is_empty or clipped.geom_type not in ["Polygon", "MultiPolygon"]:
            continue

        # calculate metrics
        solar_area = geodesic_area_m2(clipped)
        coverage_pct = (solar_area / building_area) * 100 if building_area > 0 else 0

        building["solar_area_m2"] = round(solar_area, 2)
        building["solar_coverage_pct"] = round(coverage_pct, 2)
        building["available_roof_m2"] = round(building_area - solar_area, 2)
        building["num_predictions"] = len(pred_indices)

        all_buildings_with_solar.append(building)

        if coverage_pct >= COVERAGE_THRESHOLD:
            solar_rows.append(building)

    # Save GeoJSON (only buildings with >=20% solar coverage)
    if solar_rows:
        result = gp.GeoDataFrame(solar_rows, crs=city.crs)
        output_geojson = f"results/04Results/{city_name}_Solar.geojson"
        os.makedirs(os.path.dirname(output_geojson), exist_ok=True)
        result.to_file(output_geojson, driver="GeoJSON")
    else:
        result = gp.GeoDataFrame()

    # EXCEL EXPORT WITH TWO SHEETS

    # Sheet 1: ALL buildings with solar metrics
    full = city.copy()
    full["solar_area_m2"] = 0.0
    full["solar_coverage_pct"] = 0.0
    full["available_roof_m2"] = full["area"]
    full["num_predictions"] = 0
    full["has_solar"] = False

    for b in all_buildings_with_solar:
        idx = b["_building_idx"]
        full.loc[idx, "solar_area_m2"] = b["solar_area_m2"]
        full.loc[idx, "solar_coverage_pct"] = b["solar_coverage_pct"]
        full.loc[idx, "available_roof_m2"] = b["available_roof_m2"]
        full.loc[idx, "num_predictions"] = b["num_predictions"]
        full.loc[idx, "has_solar"] = b["solar_coverage_pct"] >= COVERAGE_THRESHOLD

    # Sheet 2: Buildings affected by invalid tiles
    invalid_features_data = []

    if invalid_summary and (
        invalid_summary["total_invalid"] > 0
        or invalid_summary["total_simplified_invalid"] > 0
    ):

        all_invalid_tiles = invalid_summary.get(
            "invalid_tiles", []
        ) + invalid_summary.get("simplified_invalid_tiles", [])

        print(f"\nProcessing {len(all_invalid_tiles)} invalid tiles...")

        for invalid_tile_info in all_invalid_tiles:
            tile_path = invalid_tile_info["tile"]
            reason = invalid_tile_info["reason"]

            parts = tile_path.split("/")
            if len(parts) == 3:
                z, x, y = int(parts[0]), int(parts[1]), int(parts[2])

                tile = mercantile.Tile(x=x, y=y, z=z)
                bounds = mercantile.bounds(tile)
                tile_bbox = box(bounds.west, bounds.south, bounds.east, bounds.north)

                for idx, building in city.iterrows():
                    if building.geometry.intersects(tile_bbox):
                        intersection_geom = building.geometry.intersection(tile_bbox)
                        overlap_pct = (
                            intersection_geom.area / building.geometry.area
                        ) * 100

                        building_data = {
                            "building_index": idx,
                            "osm_id": building.get("osm_id", "N/A"),
                            "name": building.get("name", ""),
                            "fclass": building.get("fclass", ""),
                            "type": building.get("type", ""),
                            "invalid_tile": tile_path,
                            "failure_reason": reason,
                            "tile_overlap_pct": round(overlap_pct, 2),
                            "building_area_m2": round(building["area"], 2),
                            "has_solar_detected": full.loc[idx, "has_solar"],
                            "solar_coverage_pct": full.loc[idx, "solar_coverage_pct"],
                            "num_predictions": full.loc[idx, "num_predictions"],
                        }

                        invalid_features_data.append(building_data)

    # Create Excel with two sheets
    excel_path = f"results/04Results/{city_name}_Solar_buildings.xlsx"

    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        # Sheet 1: Main results
        output_df = pd.DataFrame(full.drop(columns="geometry"))
        output_df.to_excel(writer, sheet_name="output", index=False)

        # Sheet 2: Invalid features
        if invalid_features_data:
            invalid_df = pd.DataFrame(invalid_features_data)
            invalid_df_unique = invalid_df.drop_duplicates(
                subset=["building_index", "invalid_tile"]
            )
            invalid_df_unique.to_excel(
                writer, sheet_name="invalid_features", index=False
            )

            print(f"\nInvalid features sheet: {len(invalid_df_unique)} records")
            print(
                f"  Unique buildings affected: {invalid_df_unique['building_index'].nunique()}"
            )
        else:
            empty_df = pd.DataFrame(
                columns=[
                    "building_index",
                    "osm_id",
                    "name",
                    "fclass",
                    "type",
                    "invalid_tile",
                    "failure_reason",
                    "tile_overlap_pct",
                    "building_area_m2",
                    "has_solar_detected",
                    "solar_coverage_pct",
                    "num_predictions",
                ]
            )
            empty_df.to_excel(writer, sheet_name="invalid_features", index=False)

    print(
        f"\nGeoJSON saved: {output_geojson if solar_rows else 'None (no buildings >=20%)'}"
    )
    print(f"Excel saved: {excel_path}")
    print(f"  - Sheet 1 'output': {len(full)} buildings")
    print(f"  - Sheet 2 'invalid_features': {len(invalid_features_data)} records")

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("city", help="City name")
    parser.add_argument("type", help="Roof Typology: Green or Solar")
    args = parser.parse_args()

    mask_dir = os.path.join("results", "03Masks", args.type, args.city)

    if not os.path.exists(mask_dir):
        print(f"ERROR: Mask directory not found: {mask_dir}")
        print("Please run prediction first!")
        exit(1)

    intersection(args.type, args.city, mask_dir)
