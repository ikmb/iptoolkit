import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="iptkl", 
    version="0.6.7",
    author="Hesham ElAbd",
    author_email="h.elabd@ikmb.uni-kiel.de",
    description="IPTK is a library specialized in the analysis of HLA-peptidomes identified through an immunopeptidomics pipeline.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        'Documentation': 'https://iptk.readthedocs.io/en/latest/',
        'Funding': 'https://www.genes-environment-inflammation.de', 
        'Tutorials': 'https://github.com/ikmb/iptoolkit/tree/master/Tutorials',
        'Tracker': 'https://github.com/ikmb/iptoolkit/issues',
    },
    url='https://github.com/ikmb/iptoolkit',
    packages=setuptools.find_packages(),
    install_requires=['pandas','biopython','seaborn' ,'matplotlib', 'plotly' ,'mhcnames', 'pyteomics', 'h5py', 'logomaker', 'colour', 'lxml',
             'nglview', 'sklearn', 'scipy','statannot', 'numba', 'pyopenms','bokeh','holoviews','tqdm','numba','goatools','dash'], 
    classifiers=[
        'Development Status :: 3 - Alpha',
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Intended Audience :: Science/Research",
        "Typing :: Typed"
    ],
    keywords='HLA immunopeptidomics antigen-processing antigen-presentation computational-immunology interactive-data-analysis MHC',
    python_requires='>=3.7',
)
