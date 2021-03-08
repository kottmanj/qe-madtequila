import numpy as np
import json

def dict_to_molecule(geometry):
    geometry_str = f"{0} {1}\n"
    for atom in geometry["sites"]:
        geometry_str += "{} {} {} {}\n".format(
            atom["species"], atom["x"], atom["y"], atom["z"]
        )

    geometry_str += "\nunits angstrom\n"
    print(geometry_str)
