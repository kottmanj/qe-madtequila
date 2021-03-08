import sys
# activate the path to tequila on the madness-tequila docker image
sys.path.append("/usr/local/lib/python3.7/dist-packages/")
import tequila as tq
import os 
import json
import numpy

def run_madness(geometry, n_pno, mra_threshold=1.e-4, localize="boys", orthogonalization_method="cholesky", **kwargs):

    dft={"econv":mra_threshold, "localize":localize}
    if "dft" in kwargs:
        dft = {**dft, **kwargs["dft"]}

    pnoint={"orthog":orthogonalization_method}
    if "pnoint" in kwargs:
        pnoint = {**pnoint, **kwargs["pnoint"]}

    kwargs["pnoint"]=pnoint
    kwargs["dft"]=dft

    #exe = tq.quantumchemistry.QuantumChemistryMadness.find_executable("/app/madroot/")
    #mol = tq.Molecule(geometry=geometry, n_pno=n_pno, executable=exe, **kwargs)
    mol = tq.Molecule(geometry=geometry, n_pno=n_pno, **kwargs)
    return mol


###
#From here on: Only JSON stuff. Use mol_to_json and mol_from_json to serialize molecules
#Works only for the madness interface
###


class TqMadnessMoleculeEncoder(json.JSONEncoder):
    def default(self, mol):
        print("whats up ", type(mol))
        one_body_integrals = mol.compute_one_body_integrals()
        one_body_integrals = {"shape":list(one_body_integrals.shape), "data":[float(x) for x in one_body_integrals.flatten()]}
        two_body_integrals = mol.compute_one_body_integrals()
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

def mol_from_json(json_data:dict, name=None, **kwargs):
    print(type(json_data))
    print(json_data)
    if json_data is not dict:
        json_dict = json.loads(json_data)
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
    numpy.save("one_body_integrals.npy", arr=one_body_integrals)
    numpy.save("two_body_integrals.npy", arr=two_body_integrals)
    orbital_data = json_dict["orbital_data"]
    pairinfo = orbital_data["pairinfo"] 
    occinfo = orbital_data["occinfo"] 
    with open("{}_pnoinfo.txt".format(name), "w") as f:
        print("MADNESS MRA-PNO INFORMATION", file=f)
        print("pairinfo={}".format(pairinfo), file=f)
        print("nuclear_repulsion={}".format(json_dict["nuclear_repulsion"]), file=f)
        print("occinfo={}".format(occinfo), file=f)
    
    return tq.Molecule(n_pno=None, **parameters, **kwargs)
