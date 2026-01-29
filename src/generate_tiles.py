#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download only the slippy map tiles (zoom 19) that intersect building footprints.

Input:
- Vector layer of building polygons (GeoJSON / GPKG / SHP), any CRS

Output:
- Raster tiles stored as {OUT_ROOT}/{z}/{x}/{y}.png
- CSV index of downloaded tiles for reproducibility and ML pipelines

Notes:
- Geometries are buffered in meters (EPSG:3857) before tile enumeration
  to ensure complete spatial coverage near tile edges.
- Buffer size is dynamically scaled based on footprint dimensions.
"""

import os
import sys
import math
import csv
from io import BytesIO
from typing import Iterator, Optional, Set, Tuple

import geopandas as gpd
from shapely.geometry import Polygon, box as shapely_box
from shapely.ops import unary_union
from pyproj import CRS

import requests
from PIL import Image
from tqdm import tqdm

# Try mercantile, else use fallback
try:
    import mercantile

    MERCANTILE_AVAILABLE = True
except Exception:
    MERCANTILE_AVAILABLE = False


# User configuration

# Path to buildings layer (any CRS)
INPUT_PATH = r"PATH\TO\YOUR\BUILDINGS.gpkg"
LAYER_NAME = None  # if your .gpkg has multiple layers; otherwise None

# Tile source URL)
TILE_URL = "https://api.mapbox.com/styles/v1/mapbox/satellite-v9/tiles/256/{z}/{x}/{y}?access_token=YOUR_MAPBOX_TOKEN"

# Slippy zoom
ZOOM = 19

# Output directory for z/x/y.png
OUT_ROOT = r"PATH\TO\OUTPUT\TILES"

# Dynamic buffer parameters (in meters, EPSG:3857)
MIN_BUFFER_M = 10.0  # minimum buffer (small buildings)
SIZE_FRACTION = 0.10  # fraction of max(width,height) to add as buffer
MAX_BUFFER_M = 60.0  # cap to avoid excessive expansion (tune for your dataset)

# HTTP settings
TIMEOUT = 12
USER_AGENT = "Mozilla/5.0 (TilesDownloader/1.0)"
MAX_WORKERS = 8
RETRY_COUNT = 4

# If many buildings are adjacent, you can union them to reduce duplicates (optional)
USE_UNARY_UNION = False


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def out_path_png(root: str, z: int, x: int, y: int) -> str:
    folder = os.path.join(root, str(z), str(x))
    ensure_dir(folder)
    return os.path.join(folder, f"{y}.png")


def format_tile_url(tmpl: str, z: int, x: int, y: int) -> str:
    return tmpl.replace("&amp;", "&").format(z=z, x=x, y=y)


def read_buildings(path: str, layer: Optional[str] = None) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(path, layer=layer) if layer else gpd.read_file(path)

    if gdf.crs is None:
        raise ValueError("Input data has no CRS. Please assign the correct CRS first.")

    # Fix invalid geometries
    gdf["geometry"] = gdf.geometry.buffer(0)
    return gdf


def to_epsg3857(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    return gdf.to_crs(epsg=3857)


def to_epsg4326(geom) -> Polygon:
    return gpd.GeoSeries([geom], crs=3857).to_crs(4326).iloc[0]


def dynamic_buffer_m(geom_3857: Polygon) -> float:
    minx, miny, maxx, maxy = geom_3857.bounds
    max_dim = max(maxx - minx, maxy - miny)
    buf = MIN_BUFFER_M + SIZE_FRACTION * max_dim
    return max(MIN_BUFFER_M, min(buf, MAX_BUFFER_M))


def lonlat_to_tile_xy(lon: float, lat: float, zoom: int) -> Tuple[int, int]:
    n = 2**zoom
    x = int((lon + 180.0) / 360.0 * n)
    lat_rad = math.radians(lat)
    y = int((1 - math.log(math.tan(lat_rad) + 1 / math.cos(lat_rad)) / math.pi) / 2 * n)
    return x, y


def tile_bounds_4326(x: int, y: int, zoom: int):
    n = 2**zoom
    lon_w = x / n * 360 - 180
    lon_e = (x + 1) / n * 360 - 180

    def tile_y_to_lat(yy):
        s = math.pi * (1 - 2 * yy / n)
        return math.degrees(math.atan(math.sinh(s)))

    lat_n = tile_y_to_lat(y)
    lat_s = tile_y_to_lat(y + 1)
    return lon_w, lat_s, lon_e, lat_n


def tile_poly_4326_from_xy(x: int, y: int, zoom: int) -> Polygon:
    w, s, e, n = tile_bounds_4326(x, y, zoom)
    return shapely_box(w, s, e, n)


def enumerate_tiles_fallback(geom_4326: Polygon, zoom: int):
    minx, miny, maxx, maxy = geom_4326.bounds
    x0, y1 = lonlat_to_tile_xy(minx, miny, zoom)
    x1, y0 = lonlat_to_tile_xy(maxx, maxy, zoom)

    for x in range(min(x0, x1), max(x0, x1) + 1):
        for y in range(min(y0, y1), max(y0, y1) + 1):
            if tile_poly_4326_from_xy(x, y, zoom).intersects(geom_4326):
                yield zoom, x, y


def enumerate_tiles_merc(geom_4326: Polygon, zoom: int):
    minx, miny, maxx, maxy = geom_4326.bounds
    for t in mercantile.tiles(minx, miny, maxx, maxy, zooms=[zoom]):
        b = mercantile.bounds(t)
        tile_poly = shapely_box(b.west, b.south, b.east, b.north)
        if tile_poly.intersects(geom_4326):
            yield zoom, t.x, t.y


def collect_tiles(
    buildings_3857: gpd.GeoDataFrame, zoom: int
) -> Set[Tuple[int, int, int]]:
    tiles = set()
    enumerator = (
        enumerate_tiles_merc if MERCANTILE_AVAILABLE else enumerate_tiles_fallback
    )

    for geom in tqdm(buildings_3857.geometry, desc="Enumerating tiles", unit="bldg"):
        if geom is None or geom.is_empty:
            continue

        buffered = geom.buffer(dynamic_buffer_m(geom))
        geom_4326 = to_epsg4326(buffered)

        for z, x, y in enumerator(geom_4326, zoom):
            tiles.add((z, x, y))

    return tiles


def save_as_png(img_bytes: bytes, path: str):
    try:
        if img_bytes[:8] == b"\x89PNG\r\n\x1a\n":
            with open(path, "wb") as f:
                f.write(img_bytes)
        else:
            img = Image.open(BytesIO(img_bytes)).convert("RGB")
            img.save(path, format="PNG")
    except Exception:
        img = Image.open(BytesIO(img_bytes)).convert("RGB")
        img.save(path, format="PNG")


def download_tile(session, tmpl, z, x, y, timeout, retries):
    url = format_tile_url(tmpl, z, x, y)
    for _ in range(retries + 1):
        try:
            r = session.get(url, timeout=timeout)
            r.raise_for_status()
            return r.content
        except requests.RequestException:
            pass
    raise RuntimeError(f"Failed to download tile {z}/{x}/{y}")


def main():
    buildings = read_buildings(INPUT_PATH, LAYER_NAME)
    buildings_3857 = to_epsg3857(buildings)

    tiles = collect_tiles(buildings_3857, ZOOM)

    if not tiles:
        print("No tiles detected. Check geometries and buffer settings.")
        return

    ensure_dir(OUT_ROOT)

    tiles_list = sorted(list(tiles))
    csv_path = os.path.join(OUT_ROOT, "tile_index.csv")

    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    from concurrent.futures import ThreadPoolExecutor, as_completed

    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["z", "x", "y", "image_path"])

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
            futures = {}

            for z, x, y in tiles_list:
                out_path = out_path_png(OUT_ROOT, z, x, y)
                rel_path = f"{z}/{x}/{y}.png"

                writer.writerow([z, x, y, rel_path])

                if os.path.exists(out_path):
                    continue

                futures[
                    pool.submit(
                        download_tile, session, TILE_URL, z, x, y, TIMEOUT, RETRY_COUNT
                    )
                ] = out_path

            for fut in tqdm(
                as_completed(futures), total=len(futures), desc="Downloading tiles"
            ):
                try:
                    save_as_png(fut.result(), futures[fut])
                except Exception as e:
                    sys.stderr.write(f"[WARN] {e}\n")

    print(f"Finished. Tiles saved under: {OUT_ROOT}")
    print(f"Tile index written to: {csv_path}")
    print(f"Total unique tiles: {len(tiles_list)}")


if __name__ == "__main__":
    main()
