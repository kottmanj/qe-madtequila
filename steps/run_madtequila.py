import json
from qemadtequila import make_hamiltonian
from zquantum.core.utils import SCHEMA_VERSION

def run_madtequila(geometry, n_pno, transformation):
    result = make_hamiltonian(geometry=geometry, n_pno=n_pno, transformation=transformation)
    result["schema"] = SCHEMA_VERSION + "-orbitals_calc"

    with open("orbitalscalc-results.json", "w") as f:
        f.write(json.dumps(result, indent=2))
