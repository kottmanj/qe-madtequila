import sys
sys.path.append("/usr/local/lib/python3.7/dist-packages/")
import tequila as tq
import os 

def run_madness(geometry, n_pno, **kwargs):
    exe = tq.quantumchemistry.QuantumChemistryMadness.find_executable()
    print("madness exe:", exe)
    madness_root_dir = str(os.environ.get("MAD_ROOT_DIR"))
    print("MAD_ROOT_DIR=", madness_root_dir)
    mol = tq.Molecule(geometry=geometry, n_pno=n_pno, **kwargs)
    return mol
