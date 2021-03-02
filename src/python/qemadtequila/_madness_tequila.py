import sys
sys.path.append("/usr/local/lib/python3.7/dist-packages/")
import tequila as tq

def run_madness(geometry, n_pno, **kwargs):
    exe = tq.quantumchemistry.QuantumChemistryMadness.find_executable()
    print("madness exe:", exe)
    mol = tq.Molecule(geometry=geometry, n_pno=n_pno, **kwargs)
    return mol
