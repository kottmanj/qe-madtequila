import json
import qemadtequila as qemadtq
SCHEMA_VERSION="schema"

def run_madness(geometry, n_pno):
    geometry_str = None
    with open(geometry) as f:
        geometry_str = json.load(f)
    kwargs = {}
    mol = qemadtq.run_madness(geometry=geometry_str, n_pno=n_pno)
    results_dict = {}
    results_dict["schema"] = SCHEMA_VERSION + "-madresults"
    results_dict["info"] = """
    Result dictionary contains data from a tequila molecule obtained with madness as backend\n
    tequila infostring={}\n\nContent of this dictionary:\n
    one_body_integrals:the one body integrals\n
    two_body_integrals:the two body integrals in openfermion convention\n
    orbitals:information about the used orbitals\n
    nuclear_repulsion:the nuclear repulsion of the current molecule\n
    kwargs:additional keyword arguments passed to tq.Molecule""".format(str(mol))
    results_dict["one_body_integrals"] = mol.compute_one_body_integrals().tolist()
    results_dict["two_body_integrals"] = mol.compute_two_body_integrals().tolist()
    results_dict["orbitals"] = str(mol.orbitals) # workaround ... not ideal
    results_dict["nuclear_repulsion"] = mol.molecule.nuclear_repulsion

    with open("madresults.json", "w") as f:
        f.write(json.dumps(results_dict, indent=2))

def make_qubit_operator(geometry, n_pno):
    from zquantum.core.openfermion import save_interaction_operator # has import errors, probably issues with numpy>=1.20
    mol = run_madness(geometry=geometry, n_pno=n_pno)
    H = mol.make_hamiltonian()
    qubit_operator = H.to_openfermion()
    save_interaction_operator(hamiltonian, "hamiltonian.json")
