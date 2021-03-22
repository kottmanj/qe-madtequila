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

def make_qubit_hamiltonian(madmolecule, transformation="jordanwigner", **kwargs):
    # take a madmolecule and create a qubit-operator
    # export as string
    mol = madtq.mol_from_json(madmolecule, transformation=transformation, **kwargs)
    hamiltonian = mol.make_hamiltonian()
    result = {"schema":"schema"}
    # if you want openfermion strings: str(hamiltonian.to_openfermion())
    result["qubit_hamiltonian"]=str(hamiltonian)
    with open("qubit_hamiltonian.json", "w") as f:
        f.write(json.dumps(result, indent=2))

def load_json(x:str):
    if ".json" in x:
        with open(x) as f:
            y = json.load(f)
        return y
    elif isinstance(x, str):
        return json.loads(x)
    else:
        return x

def optimize_measurements(qubit_hamiltonian:str, circuit:str=None):
    """
    Use the tequila implementation of the grouping Algorithm from T.C Yen et.al. 
    circuit: open-qasm-2 string; if no circuit is given we will create an empty one
    hamiltonian: openfermion::QubitOperator or tequila hamiltonian as string

    Rerurn: A json dictionary with "measurement_count":total number of optimized measurments
    "groupings":a list of dictionary containing 'circuit' (qasm str) and 'hamiltonian':openfermion_string
    if no circuit was given, the circuits in the grouping list contain only the basis changes
    """
    
    # get strings from json
    hamiltonian = load_json(qubit_hamiltonian)["qubit_hamiltonian"]
     
    # convert to tq objects
    try: 
        H = madtq.tq.QubitHamiltonian.from_string(hamiltonian, openfermion_format=True)
    except:
        H = madtq.tq.QubitHamiltonian.from_string(hamiltonian)
    
    if circuit is None:
        U = madtq.tq.QCircuit()
    else:
        circuit = load_json(circuit)["circuit"]
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
    result = {"schema":"schema", "measurement_count":E.count_measurements()}
    groups = []
    for expv in E.get_expectationvalues():
        # summation doesn't do much here
        # note that this hamiltonian is an all-z hamiltonian (one measurement)
        h = sum(expv.H, 0.0)
        h = h.to_openfermion()
        u = madtq.tq.export_open_qasm(expv.U)
        groups.append({"circuit":u,"hamiltonian":str(h)})
    result["groups"]=groups
    with open("measurement_groups.json", "w") as f:
        f.write(json.dumps(result, indent=2))

    return result

def compute_pyscf_energy(madmolecule, method="fci", **kwargs):
    # Step needs pyscf installed (will add that in next update of the image, currently it needs to be added to requirements)
    # not 100% sure this does what it's supposed to ... maybe be carful
    mol = madtq.mol_from_json(madmolecule, **kwargs)
    energy = madtq.compute_pyscf_energy(mol, method=method, **kwargs)
    result = {"SCHEMA":"schema",
            "info":"{} - {}/MRA-PNO({},{})".format(mol.parameters.name, "method", mol.n_electrons, 2*mol.n_orbitals),
            "energy":energy}
    with open("energy.json", "w") as f:
        f.write(json.dumps(result, indent=2))
    return energy


if __name__ == "__main__":
    run_madness("Li 0.0 0.0 0.0\n H 0.0 0.0 1.6", 1, name="lih")
    X=compute_pyscf_energy(madmolecule="madmolecule.json", method="fci")
    print(X)
    Y=compute_pyscf_energy(madmolecule="madmolecule.json", method="ccsd")
    print(Y)
    Z=compute_pyscf_energy(madmolecule="madmolecule.json", method="ccsd(t)")
    print(Z)
    A=compute_pyscf_energy(madmolecule="madmolecule.json", method="all")
    print(A)
    
    mol = madtq.mol_from_json("madmolecule.json")
    v, vv = numpy.linalg.eigh(mol.make_hamiltonian().to_matrix())
    # GS is 3-electron state
    for i in range(5):
        print(v[i])
        print(madtq.tq.QubitWaveFunction(vv[:,i]))
    #compute_pno_upccd(madmolecule="madmolecule.json")
    #U = madtq.tq.gates.Ry(angle="a", target=0) + madtq.tq.gates.CNOT(0,1)
    #qasm = madtq.tq.export_open_qasm(U, variables={"a":1.0})
    #optimize_measurements(circuit={"circuit":qasm}, qubit_hamiltonian={"qubit_hamiltonian":"1.0*X(0)+2.0*X(0)Y(1)"})
    #asd=optimize_measurements(qubit_hamiltonian={"qubit_hamiltonian":"1.0*X(0)+2.0*X(0)Y(1)"})
    #print(asd["measurement_count"])
