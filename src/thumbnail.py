import os
import math
import requests
from io import BytesIO
from PIL import Image, ImageDraw
from openpyxl import load_workbook
from openpyxl.drawing.image import Image as ExcelImage
from qgis.core import (
    QgsProject,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform
)

# === CONFIGURATION ===
layer_name = "test17"
zoom = 17
tile_size = 256
grid_size = 3
stitched_size = tile_size * grid_size
tile_url_template = "https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}"  # Hybrid layer

output_base = "C:/Users/junwei.quay/Documents/qgis_japan/kanto_buildings/kanto_10000/test17/img"
excel_path = "C:/Users/junwei.quay/Documents/qgis_japan/kanto_buildings/kanto_10000/test17/test17.xlsx"
output_excel = "C:/Users/junwei.quay/Documents/qgis_japan/kanto_buildings/kanto_10000/test17/test17_tn.xlsx"
thumbnail_column_letter = "H"
id_column_index = 1  # Column index for osm_id

os.makedirs(output_base, exist_ok=True)

# === HELPER FUNCTIONS ===
def latlon_to_tile(lat, lon, zoom):
    lat_rad = math.radians(lat)
    n = 2.0 ** zoom
    x_tile = (lon + 180.0) / 360.0 * n
    y_tile = (1.0 - math.log(math.tan(lat_rad) + 1.0 / math.cos(lat_rad)) / math.pi) / 2.0 * n
    return x_tile, y_tile

def latlon_to_pixel_offset(lat, lon, zoom, center_tile_x, center_tile_y):
    x_tile_float, y_tile_float = latlon_to_tile(lat, lon, zoom)
    dx = (x_tile_float - (center_tile_x - 1)) * tile_size
    dy = (y_tile_float - (center_tile_y - 1)) * tile_size
    return int(dx), int(dy)

# === LOAD QGIS LAYER ===
layer = QgsProject.instance().mapLayersByName(layer_name)[0]
crs_src = layer.crs()
crs_dest = QgsCoordinateReferenceSystem("EPSG:4326")
xform = QgsCoordinateTransform(crs_src, crs_dest, QgsProject.instance())

# === LOAD EXCEL ===
wb = load_workbook(excel_path)
ws = wb.active

# === BUILD FEATURE LOOKUP BY osm_id ===
feature_dict = {}
for feature in layer.getFeatures():
    osm_id = str(feature["osm_id"])
    feature_dict[osm_id] = feature

# === LOOP THROUGH EXCEL ROWS ===
headers = {"User-Agent": "Mozilla/5.0"}

for row_idx in range(2, ws.max_row + 1):
    osm_id = str(ws.cell(row=row_idx, column=id_column_index).value)
    if osm_id not in feature_dict:
        print(f"⚠️ No matching feature for osm_id {osm_id}")
        continue

    feature = feature_dict[osm_id]
    geom = feature.geometry().centroid().asPoint()
    wgs_point = xform.transform(geom)
    lat, lon = wgs_point.y(), wgs_point.x()

    center_tile_x, center_tile_y = map(int, latlon_to_tile(lat, lon, zoom))
    stitched_img = Image.new("RGB", (stitched_size, stitched_size))

    success = True
    for dx in range(-1, 2):
        for dy in range(-1, 2):
            tile_x = center_tile_x + dx
            tile_y = center_tile_y + dy
            tile_url = tile_url_template.format(x=tile_x, y=tile_y, z=zoom)
            try:
                response = requests.get(tile_url, headers=headers)
                if response.status_code == 200:
                    tile_img = Image.open(BytesIO(response.content)).convert("RGB")
                    stitched_img.paste(tile_img, ((dx + 1) * tile_size, (dy + 1) * tile_size))
                else:
                    print(f"❌ Failed to download tile {tile_x}, {tile_y}")
                    success = False
            except Exception as e:
                print(f"❌ Error downloading tile {tile_x}, {tile_y}: {e}")
                success = False

    if not success:
        ws.cell(row=row_idx, column=8).value = "Tile download failed"
        continue

    # Calculate pixel offset of centroid within stitched image
    pixel_x, pixel_y = latlon_to_pixel_offset(lat, lon, zoom, center_tile_x, center_tile_y)

    # Draw red cross at centroid pixel
    draw = ImageDraw.Draw(stitched_img)
    cross_size = 8
    draw.line((pixel_x - cross_size, pixel_y, pixel_x + cross_size, pixel_y), fill="red", width=3)
    draw.line((pixel_x, pixel_y - cross_size, pixel_x, pixel_y + cross_size), fill="red", width=3)

    # Save image
    output_path = os.path.join(output_base, f"{osm_id}.png")
    stitched_img.save(output_path)
    print(f"✅ Saved and annotated image for {osm_id}")

    # Insert into Excel
    try:
        excel_img = ExcelImage(output_path)
        excel_img.width = 200
        excel_img.height = 200
        ws.add_image(excel_img, f"{thumbnail_column_letter}{row_idx}")
        print(f"✅ Inserted thumbnail for {osm_id}")
    except Exception as e:
        print(f"❌ Failed to insert image into Excel for {osm_id}: {e}")
        ws.cell(row=row_idx, column=8).value = "Image not found"

# === SAVE EXCEL ===
wb.save(output_excel)
print("🎉 All tiles stitched, annotated, and inserted into Excel.")