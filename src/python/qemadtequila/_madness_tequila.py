import sys
# activate the path to tequila on the madness-tequila docker image
sys.path.append("/usr/local/lib/python3.7/dist-packages/")
import tequila as tq
import os
import json
import numpy
SCHEMA_VERSION="schema"

def run_madness(geometry, n_pno, mra_threshold=1.e-4, localize="boys", orthogonalization_method="cholesky", diagonal=True, maxrank=None, **kwargs):
    """
    geometry: a string that either points to an xyz file (e.g. geometry="my_file.xyz") or carries the molecular structure
              in xyz syntax (without the preamble, e.g. geometry="H 0.0 0.0 0.0\nH 0.0 0.0 0.7")
    n_pno: number of pnos to compute (number of qubits is then 2*n_pno + n_electrons)
    mra_threshold: MRA accuracy threshold (1.e-4 is usually fine)
    localize: localization method ("boys", "pm", "canon")
    orthogonalization_method: Global orthogonalization of the PNOs ("cholesky", "symmetric")
    diagonal: Apply a diagonal approximation to the MP2 surrogate (only compute diagonal pairs)
    maxrank: maximum number of PNOs computed for each MP2 pair (not the number picked in the end; in some cases it's good to be able to control that)
    """
    dft={"econv":mra_threshold, "localize":localize}
    if "dft" in kwargs:
        dft = {**dft, **kwargs["dft"]}

    pnoint={"orthog":orthogonalization_method}
    if "pnoint" in kwargs:
        pnoint = {**pnoint, **kwargs["pnoint"]}

    pno={"thresh":mra_threshold, "diagonal":diagonal}
    if "pno" in kwargs:
        pno = {**pno, **kwargs["pno"]}

    if maxrank is not None:
        pno["maxrank"]=maxrank

    kwargs["pnoint"]=pnoint
    kwargs["dft"]=dft
    kwargs["pno"]=pno


    #exe = tq.quantumchemistry.QuantumChemistryMadness.find_executable("/app/madroot/")
    #mol = tq.Molecule(geometry=geometry, n_pno=n_pno, executable=exe, **kwargs)
    mol = tq.Molecule(geometry=geometry, n_pno=n_pno, **kwargs)
    return mol

def fetch_integrals(mol, two_body_ordering="mulliken", *args, **kwargs):
    # small helper function
    # will be integrated in tq for next versions
    if mol.active_space is not None and len(mol.active_space.frozen_reference_orbitals)>0:
        c, h1, h2 = mol.molecule.get_active_space_integrals(active_indices=mol.active_space.active_orbitals,
                                                        occupied_indices=mol.active_space.frozen_reference_orbitals)
    else:
        c=0.0
        h1 = mol.compute_one_body_integrals()
        h2 = mol.compute_two_body_integrals()
    eri = tq.quantumchemistry.NBodyTensor(h2, ordering="openfermion")
    eri = eri.reorder(to=two_body_ordering).elems

    return c, h1, eri

def compute_fci(mol, *args, **kwargs):
    from pyscf import fci
    c, h1, h2 = fetch_integrals(mol)
    norb = mol.n_orbitals
    nelec = mol.n_electrons
    e, fcivec = fci.direct_spin1.kernel(h1, h2, norb, nelec, **kwargs)
    return e + c + mol.molecule.nuclear_repulsion

def compute_ccsd(mol, **kwargs):
    import pyscf
    c, h1, h2 = fetch_integrals(mol)
    norb = mol.n_orbitals
    nelec = mol.n_electrons
    
    mo_coeff = numpy.eye(norb)
    mo_occ = numpy.zeros(norb)
    mo_occ[:nelec//2] = 2
    
    # pray this works
    hf_dummy = pyscf.scf.RHF(pyscf.gto.M())
    hf_dummy._eri = h2
    hf_dummy.get_hcore = lambda *args: h1
    from pyscf import cc
    ccsd = cc.ccsd.CCSD(hf_dummy, mo_coeff=mo_coeff, mo_occ=mo_occ, **kwargs)
    ccsd.diis_start_cycle=0 # following recommendation from pyscf doc.
    ccsd.diis_start_energy_diff = 1e2
    ccsd.kernel()
    
    return ccsd.e_tot


###
#From here on: Only JSON stuff. Use mol_to_json and mol_from_json to serialize molecules
#Works only for the madness interface
###


class TqMadnessMoleculeEncoder(json.JSONEncoder):
    def default(self, mol):
        one_body_integrals = mol.compute_one_body_integrals()
        one_body_integrals = {"shape":list(one_body_integrals.shape), "data":[float(x) for x in one_body_integrals.flatten()]}
        two_body_integrals = mol.compute_two_body_integrals()
        two_body_integrals = {"shape":list(two_body_integrals.shape), "data":[float(x) for x in two_body_integrals.flatten()]}
        nuc_rep = float(mol.molecule.nuclear_repulsion)
        orbital_data = self.encode_pnoinfo(mol.orbitals)
        parameters = mol.parameters.__dict__
        parameters["basis_set"]="madness"
        return {"one_body_integrals":one_body_integrals, "two_body_integrals":two_body_integrals, "nuclear_repulsion":nuc_rep, "orbital_data":orbital_data, "name":mol.parameters.name, "parameters":parameters}

    def encode_pnoinfo(self, orbitals):
        pairinfo=""
        occinfo=""
        for orbital in orbitals:
           p = orbital.pno_pair
           if len(p)==1:
               pairinfo+="{},".format(*p)
           elif len(p)==2:
               pairinfo+="{}.{},".format(*p)
           occinfo += "{},".format(orbital.occ)

        return {"pairinfo":pairinfo.rstrip(","), "occinfo":occinfo.rstrip(",")}

def mol_to_json(mol):
    return json.dumps(mol, indent=2, cls=TqMadnessMoleculeEncoder)

def mol_from_json(json_data:str, name=None, **kwargs):

    if hasattr(json_data, "lower") and ".json" in json_data.lower():
        with open(json_data, "r") as f:
            json_dict = json.load(f)
            print(type(json_dict))
            print(json_dict)
    elif hasattr(json_data, "lower()"):
        json_dict = json.loads(json_data)
    else:
        json_dict = json_data

    if "mol" in json_dict:
        json_dict=json.loads(json_dict["mol"])

    parameters = json_dict["parameters"]
    if name is None:
        name=parameters["name"]
    else:
        parameters["name"]=name
    del parameters["multiplicity"]


    obi_data=json_dict["one_body_integrals"]
    one_body_integrals=numpy.asarray(obi_data["data"], dtype=float).reshape(obi_data["shape"])
    tbi_data=json_dict["two_body_integrals"]
    two_body_integrals=numpy.asarray(tbi_data["data"], dtype=float).reshape(tbi_data["shape"])
    numpy.save("{}_htensor.npy".format(name), arr=one_body_integrals)
    numpy.save("{}_gtensor.npy".format(name), arr=two_body_integrals)
    orbital_data = json_dict["orbital_data"]
    pairinfo = orbital_data["pairinfo"]
    occinfo = orbital_data["occinfo"]
    with open("{}_pnoinfo.txt".format(name), "w") as f:
        print("MADNESS MRA-PNO INFORMATION", file=f)
        print("pairinfo={}".format(pairinfo), file=f)
        print("nuclear_repulsion={}".format(json_dict["nuclear_repulsion"]), file=f)
        print("occinfo={}".format(occinfo), file=f)
    mol = tq.Molecule(n_pno=None, **parameters, **kwargs)
    return mol


def save_interaction_operator(interaction_operator, filename) -> None:
    """Save an interaction operator to file.
    Args:
        interaction_operator (InteractionOperator): the operator to be saved
        filename (str): the name of the file
    """

    with open(filename, "w") as f:
        f.write(
            json.dumps(convert_interaction_op_to_dict(interaction_operator), indent=2)
        )

def convert_interaction_op_to_dict(op) -> dict:
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
