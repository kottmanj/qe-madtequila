import numpy as np
import json

def dict_to_molecule(molgeometry):
    geometry = None
    with open(molgeometry) as f:
        geometry = json.load(f)
    print("here is the fucking geometery: ",geometry)
