import numpy as np
import json

def dict_to_molecule(geometry):
    geometry_dict = None
    with open("molecule.json", "r") as f:
        geometry_dict = json.loads(geometry)

    print(geometry_dict)
