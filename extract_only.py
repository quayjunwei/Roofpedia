import os
import argparse
from src.extract import intersection

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