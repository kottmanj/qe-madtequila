import numpy as np
import json

def dict_to_molecule(molgeometry):
    geometry = None
    with open(molgeometry) as f:
        geometry = json.load(f)

    geometry_str = ""
    for atom in geometry["sites"]:
        geometry_str += "{} {} {} {}\n".format(
            atom["species"], atom["x"], atom["y"], atom["z"]
        )
    geometry_str["schema"] = "string"
    with open('geometry_str.json', 'w') as f:
        f.write(json.dumps(geometry_str))
