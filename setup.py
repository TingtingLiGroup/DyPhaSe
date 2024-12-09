from setuptools import setup, find_namespace_packages
from os import path

DIR = path.dirname(path.abspath(__file__))
INSTALL_PACKAGES = open(path.join(DIR,  'requirements.txt'), encoding='utf-8').read().splitlines()

setup(
    name="DyPhaSe",
    version="0.1",
    packages=find_namespace_packages(),
    include_package_data=True,
    install_requires=INSTALL_PACKAGES,
    python_requires='>=3.9',
    #extras_require={
    #    "R_integration": ["Ensure R >= 4.4.0 and install required R packages: Mfuzz, dplyr, Biobase"]
    #    }
    package_data={
        'DyPhaSe':['data/mouse_uniprot_gene_ensembl.csv','data/human_uniprot_gene_ensembl.csv','data/mouse_LLPS_score.csv','data/human_LLPS_score.csv','data/10090.string.ppi.genename.csv','data/9606.string.ppi.genename.csv']
    }
)
