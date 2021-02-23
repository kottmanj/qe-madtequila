import tequila as tq

def make_hamiltonian(geometry, n_pno, **kwargs):

    mol = tq.Molecule(geometry=geometry, n_pno=n_pno, **kwargs)

    results_dict["info"] = "Result dictionary contains data from a tequila molecule obtained with madness as backend\n
    tequila infostring={}\n\nContent of this dictionary:\n
    hamiltonian:openfermion::QubitOperator\n
    one_body_integrals:the one body integrals\n
    two_body_integrals:the two body integrals in openfermion convention\n
    orbitals:information about the used orbitals\n
    nuclear_repulsion:the nuclear repulsion of the current molecule\n
    kwargs:additional keyword arguments passed to tq.Molecule".format(str(mol))
    results_dict["hamiltonian"] = H.to_openfermion()
    results_dict["one_body_integrals"] = mol.compute_one_body_integrals()
    results_dict["two_body_integrals"] = mol.compute_two_body_integrals()
    results_dict["orbitals"] = mol.orbitals
    results_dict["nuclear_repulsion"] = mol.molecule.nuclear_repulsion
    results_dict["kwargs"] = kwargs

    return results_dict
