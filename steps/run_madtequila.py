import json
import qemadtequila as qemadtq
SCHEMA_VERSION="schema" 

def run_madness(geometry, n_pno, **kwargs):

    mol = qemadtq.run_madness(geometry=geometry, n_pno=n_pno, **kwargs)
    results_dict = {}
    results_dict["schema"] = SCHEMA_VERSION + "-madresults"
    results_dict["kwargs"] = kwargs
    results_dict["geometry"] = geometry
    results_dict["n_pno"] = n_pno
    json_string = qemadtq.mol_to_json(mol)
    results_dict["mol"]=json_string
    with open("madmolecule.json", "w") as f:
        f.write(json.dumps(results_dict, indent=2))

def make_qubit_operator(madmolecule, transformation="JordanWigner", **kwargs):
    # madmolecule is the result of run_madness
    # re-initialize tq molecule
    mol = qemadtq.mol_from_json(madmolecule, transformation="JordanWigner", **kwargs)
    H = mol.make_hamiltonian()

def compute_pno_upccd(madmolcule, **kwargs):
    # madmolecule is the result of run_madness
    # re-initialize tq molecule
    mol = qemadtq.mol_from_json(madmolecule, transformation="JordanWigner", **kwargs)
    H = mol.make_hamiltonian()
    U = mol.make_pno_upccgsd_ansatz()
    E = tq.ExpectationValue(H=H, U=U)
    result = tq.minimize(E, **kwargs)

    energy = result.energy
    with open("final_energy", "w") as f:
        f.write(json_dumps(energy, indent=2))
