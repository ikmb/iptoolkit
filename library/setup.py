import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="IPTKL", 
    version="0.4.3",
    author="Hesham ElAbd",
    author_email="h.elabd@ikmb.uni-kiel.de",
    description="IPTK is a Pythonic library specialized in the analysis of HLA-peptidomes identified through an Immunopeptiomics pipeline.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ikmb/iptoolkit",
    packages=setuptools.find_packages(),
    install_requires=['numpy','pandas','biopython','seaborn' ,'matplotlib', 'plotly' ,'mhcnames', 'pyteomics', 'h5py', 'logomaker', 'colour', 'lxml', 'nglview', 'sklearn', 'scipy','statannot'], 
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
