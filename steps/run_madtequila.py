import json
import qemadtequila as madtq
from os import PathLike
from typing import Union, Dict
import numpy
SCHEMA_VERSION=madtq.SCHEMA_VERSION

AnyPath = Union[str, bytes, PathLike]

def run_madness(geometry, n_pno, **kwargs):
    """
    wrapper function that takes care of json i/o business
    geometry: Either a json filename or a geometry string
              json should encode a dictionary of the form
              {"sited": {[{"species":"H", "x":0.0, "y":0.0, "z":0.0}, {"species":"H", "x":0.0, "y":0.0, "z":0.7}, ...]} }
              the geometry string is in xyz format without the preamble e.g.
              "H 0.0 0.0 0.0\nH 0.0 0.0 0.7"
    n_pno: number of pnos that shall be used for the qubit hamiltonian
    kwargs: other keyword arguments for the tequila-madness interface (see also madtq.run_madness)
    """
    
    molgeometry=None
    geometry_str = ""
    if "sites" in geometry:
        molgeometry = geometry
    elif ".json" in geometry:
        with open(geometry) as f:
            molgeometry = json.load(f)
    else:
        geometry_str = geometry
    
    if molgeometry is not None:
        for atom in molgeometry["sites"]:
            geometry_str += "{} {} {} {}\n".format(
                atom["species"], atom["x"], atom["y"], atom["z"]
            )

    kwargs = {}
    mol = madtq.run_madness(geometry=geometry_str, n_pno=n_pno, **kwargs)
    results_dict = {}
    results_dict["schema"] = SCHEMA_VERSION + "-madresults"
    results_dict["kwargs"] = kwargs
    results_dict["geometry"] = geometry
    results_dict["n_pno"] = n_pno
    json_string = madtq.mol_to_json(mol)
    results_dict["mol"]=json_string
    with open("madmolecule.json", "w") as f:
        f.write(json.dumps(results_dict, indent=2))

def compute_pno_upccd(madmolecule, **kwargs):
    """
    Small example of a VQE in tequila using the previously computed madness molecule
    madmolecule: result of run_madness given as json filename
    kwargs: further arguments for the tq.minimize function
    """
    # madmolecule is the result of run_madness
    # re-initialize tq molecule
    mol = madtq.mol_from_json(madmolecule, transformation="JordanWigner", **kwargs)
    H = mol.make_hamiltonian()
    U = mol.make_pno_upccgsd_ansatz()
    E = madtq.tq.ExpectationValue(H=H, U=U)
    result = madtq.tq.minimize(E, **kwargs)

    energy = result.energy
    with open("final_energy.json", "w") as f:
        f.write(json.dumps(energy, indent=2))

def make_interaction_operator(madmolecule, **kwargs):
    #from zquantum.core.openfermion import save_interaction_operator # import problems in combination with custom image, use the function only with standard runtime
    mol = madtq.mol_from_json(madmolecule, transformation="JordanWigner", **kwargs)
    hamiltonian = mol.make_molecular_hamiltonian()
    # leaving this here since it might be useful to know
    # h = mol.compute_one_body_integrals() # actually get function in this case
    # g = mol.compute_two_body_integrals() # same

    # there is an issue with importing ans using the save_interaction_operator function
    # from zquantum.core.openfermion, due to some depencdencies so we copied the functions
    # here to work around that
    madtq.save_interaction_operator(hamiltonian, "hamiltonian.json")

def optimize_measurements(circuit:str, hamiltonian:str):
    """
    Use the tequila implementation of the grouping Algorithm from T.C Yen et.al. 
    circuit: open-qasm-2 string
    hamiltonian: openfermion::QubitOperator or tequila hamiltonian as string

    Rerurn: A json dictionary with "measurement_count":total number of optimized measurments
    "groupings":a list of dictionary containing 'circuit' (qasm str) and 'hamiltonian':openfermion_string
    """

    try: 
        H = madtq.tq.QubitHamiltonian.from_string(hamiltonian, openfermion_format=True)
    except:
        H = madtq.tq.QubitHamiltonian.from_string(hamiltonian)
    U = madtq.tq.import_open_qasm(circuit)

    # optimize_measurements will decompose the expectation value
    # into a sum of expectation values with transformed circuits
    # and Hamiltonians that are build from Pauli-Z only
    # in tequila this objective can be used like any other
    # in the following we will extract it's components in order to 
    # potentially use them outside of tequila
    E = madtq.tq.ExpectationValue(H=H, U=U, optimize_measurements=True)

    # now pull the circuits out and give them back as qasm lists
    # note that the measurement optimization will change the circuits (adding basis rotations)
    result = {"schema":"schema", "measurement_count":E.count_expectationvalues()}
    groups = []
    for expv in E.get_expectationvalues():
        # summation doesn't do much here
        # note that this hamiltonian is an all-z hamiltonian (one measurement)
        h = sum(expv.H, 0.0)
        h = h.to_openfermion()
        u = madtq.tq.export_open_qasm(expv.U)
        groups.append({"circuit":u,"hamiltonian":str(h)})
    result["groups"]=groups
    with open("groupings.json", "w") as f:
        f.write(json.dumps(result, indent=2))


    





if __name__ == "__main__":
    #run_madness("he 0.0 0.0 0.0", 1)
    #compute_pno_upccd(madmolecule="madmolecule.json")
    U = madtq.tq.gates.Ry(angle="a", target=0) + madtq.tq.gates.CNOT(0,1)
    qasm = madtq.tq.export_open_qasm(U, variables={"a":1.0})
    optimize_measurements(circuit=qasm, hamiltonian="1.0*X(0)+2.0*X(0)Y(1)")
