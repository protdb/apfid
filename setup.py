from setuptools import setup, find_packages

with open("README.MD", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="apfid",
    version="0.1.0",
    author="Laboratory of structural proteomics in the IBMC",
    author_email="kirill.s.nikolsky@yandex.ru",
    description="APFID - Arbitrary Protein Fragment IDentifier parser",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/protdb/apfid",
    project_urls={
        "Bug Tracker": "https://github.com/protdb/apfid/issues",
        "Source Code": "https://github.com/protdb/apfid",
        "Documentation": "https://github.com/protdb/apfid",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=[
    ],
    keywords="protein structure bioinformatics fragment identifier pdb alphafold",
    zip_safe=False,
)
