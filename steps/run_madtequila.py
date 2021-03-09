import json
import qemadtequila as qemadtq
from openfermion import (
    InteractionOperator,
    QubitOperator,
    IsingOperator,
    SymbolicOperator,
    InteractionRDM,
)
from os import PathLike
from typing import Union, Dict
import numpy
SCHEMA_VERSION="schema"

AnyPath = Union[str, bytes, PathLike]

def run_madness(geometry, n_pno):
    geometry_str = None
    with open(geometry) as f:
        geometry_str = json.load(f)
    kwargs = {}
    mol = qemadtq.run_madness(geometry=geometry_str, n_pno=n_pno)
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

def compute_pno_upccd(madmolecule, **kwargs):
    # madmolecule is the result of run_madness
    # re-initialize tq molecule
    mol = qemadtq.mol_from_json(madmolecule, transformation="JordanWigner", **kwargs)
    H = mol.make_hamiltonian()
    U = mol.make_pno_upccgsd_ansatz()
    E = qemadtq.tq.ExpectationValue(H=H, U=U)
    result = qemadtq.tq.minimize(E, **kwargs)

    energy = result.energy
    with open("final_energy.json", "w") as f:
        f.write(json.dumps(energy, indent=2))

def save_interaction_operator(interaction_operator: InteractionOperator, filename: AnyPath) -> None:
    """Save an interaction operator to file.
    Args:
        interaction_operator (InteractionOperator): the operator to be saved
        filename (str): the name of the file
    """

    with open(filename, "w") as f:
        f.write(
            json.dumps(convert_interaction_op_to_dict(interaction_operator), indent=2)
        )

def convert_interaction_op_to_dict(op: InteractionOperator) -> dict:
    """Convert an InteractionOperator to a dictionary.
    Args:
        op (openfermion.ops.InteractionOperator): the operator
    Returns:
        dictionary (dict): the dictionary representation
    """

    dictionary = {"schema": SCHEMA_VERSION + "-interaction_op"}
    dictionary["constant"] = convert_array_to_dict(numpy.array(op.constant))
    dictionary["one_body_tensor"] = convert_array_to_dict(numpy.array(op.one_body_tensor))
    dictionary["two_body_tensor"] = convert_array_to_dict(numpy.array(op.two_body_tensor))

    return dictionary

def convert_array_to_dict(array: numpy.ndarray) -> dict:
    """Convert a numpy array to a dictionary.
    Args:
        array (numpy.array): a numpy array
    Returns:
        dictionary (dict): the dict containing the data
    """

    dictionary = {}
    if numpy.iscomplexobj(array):
        dictionary["real"] = array.real.tolist()
        dictionary["imag"] = array.imag.tolist()
    else:
        dictionary["real"] = array.tolist()

    return dictionary

def make_qubit_operator(madmolecule, **kwargs):
    #from zquantum.core.openfermion import save_interaction_operator # import problems in combination with custom image, use the function only with standard runtime
    mol = qemadtq.mol_from_json(madmolecule, transformation="JordanWigner", **kwargs)
    H = mol.make_hamiltonian()
    # leaving this here since it might be useful to know
    # h = mol.compute_one_body_integrals() # actually get function in this case
    # g = mol.compute_two_body_integrals() # same
    qubit_operator = H.to_openfermion()
    save_interaction_operator(hamiltonian, "hamiltonian.json")

if __name__ == "__main__":
    run_madness("he 0.0 0.0 0.0", 1)
    compute_pno_upccd(madmolecule="madmolecule.json")
