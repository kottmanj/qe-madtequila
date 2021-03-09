import setuptools
import os

## failsafe (no better idea currently)
#try:
#    import tequila as tq
#    mad_exe = tq.quantumchemistry.QuantumChemistryMadness.find_executable()
#    if mad_exe is None:
#        raise Exception("pno_integrals executable not found")
#except Exception as E:
#    raise Exception("{}\n...something is wrong with the environment, make sure to use the madness-tequila docker image.".format(E))

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="qe-madtequila",
    version="1.0",
    author="Jakob S. Kottmann",
    author_email="jakob.kottmann@utoronto.ca",
    description="Madness PNO Integration for Quantum Engine",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="tba",
    packages=["qemadtequila"],
    package_dir={"": "src/python"},
    classifiers=(
        "Programming Language :: Python :: 3",
        "Operating System :: Hopefully OS Independent",
    ),
    install_requires=["ruamel.yaml"],

)
