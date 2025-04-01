from apfid import parse_apfid

# APFID

APFID stands for **A**rbitrary **P**rotein **F**ragment **ID**entifier and serves as unique identificator of molecular structure fragment in
databases such as PDB, AlphaFold or user own database if necessary. It can be found useful to identify and store
data in the researches dealing with structural motifs or domains on a big scale. It is used specifically by
[PSSKB.ORG](https://psskb.org) database as structure identifier and alo in some structural researches performed by
[Laboratory of structural proteomics in the IBMC, Moscow](https://en.ibmc.msk.ru/departments?view=article&id=14:biobanking-group).

This library contains basic parser utility to work with APFIDs using Python 3.10+.

## APFID structure:

commonly, it looks like this: `1YSI_A111_A191` but there are some options described below

Apfid contains of three underscore-separated parts:
* experiment_id: Identificator in Database: PDB ID, Alphafold ID or user-given one
* Start position: chain ID and residue no
* End position: chain ID and residue no

Note: chain ID must represent code used in PDB structure file, not in deposition details.
At rcsb.org the correct ID is described as _Auth ID_.

To identify the whole chain, one can use short form with `{experiment_id}_{chain_id}` and ommit what follows. So `1EDI_A`
is absolutely correct APFID.

In some cases, second underscore can be replaced with "-", so `1YSI_A111_A191` is equal to  `1YSI_A111-A191`

### APFIDv2

Is more forgiving yet not used before in our structures way to identify structure. Basic rules are:

`{experiment_id}[:{model}]_{chain_id}[{start}_[{chain2_id}_]{end}]`

where:
* experiment_id: Identificator in Database: PDB ID, Alphafold ID or user-given one
* model: ID of the model in the same file (can be ommited to fallback to zero, or first model). Added to navigate through NMR and MD
* chain_id: ID of the chain to begin
* chain2_id: ID of the chain to end. if differs, all the chains alphabetically within chain_id and chain2_id are included. _not fully implemented yet as no cases are really imaginable_
* start, end: numbers of residues in chains

examples:

```
1YSI_A111_A191
1YSI:10_A111_191
1YSI:10_A111-191
1YSI:10_A
1YSI_A
```

### Databases ID:

ID is register-insensitive, e.g. APFID `1EDI_A` is equal to `1edi_A`
(but none are equal to `1edi_a` as register can be crucial to chain ID)

User databases experiment_ids (for example, structures uploaded by users in PSSKB services) should have prefix 'USR' and
total length not equal to 4. Letters, numders and hyphens are allowed.

## INSTALLATION

```shell
pip install apfid
```

## USAGE

Parsing from the string:
```python
from apfid import parse_apfid
apf = parse_apfid('1YSI_A111_A191')
print(apf)
# 1YSI_A111_A191
print(apf.upper())
# 1YSI_A111_A191
print(apf.lower())
# 1ysi_A111_A191
print(apf.experiment_id)
# 1YSI
print(apf.chain_id, apf.start, apf.end)
# A 111 191
```

Parameters parsing:
```python
from apfid import Apfid
apf = Apfid('1YSQ', 'A', 111, 191, 10, version=2)
apf.upper()
# '1YSQ:10_A111_191'
```

## Alphafold compatibility

AlphaFold identifier contains three valuable parameters:

- UniProt ID of the protein
- File number for long proteins (e.g. _F1_)
- Structure version (e.g. _v4_)

File number and version can be set to default values, 1 and 4 respectively.

Alphafold ID can be passed to APFID in one of the forms partially specifying parameters:

```python
from apfid import parse_apfid
parse_apfid("AF-A1KHV0-F1-V4_A") # Full PSSKB-like form
parse_apfid("AF-A1KHV0_A")
```

Alphafold ID is processed by separate class `AlphafoldId`. For AF structures, its instance is accessible from Apfid.af_id.
It also provides methods to represent ID in specific forms:

```python
from apfid import parse_apfid
apfid = parse_apfid("AF-A1KHV0-F1-V4_A")
print(apfid.source)
# Alphafold
print(type(apfid.af_id)) # is
# <class 'apfid.AlphafoldId'>
print(apfid.af_id)
# AF-A1KHV0-F1-v4
print(apfid.af_id.get_dl_id()) # can substitute in download
# AF-A1KHV0-F1-model_v4
print(apfid.af_id.get_psskb_id()) # as in PSSKB.ORG database 
# AF-A1KHV0-F1-V4
print(apfid.af_id.uniprot_id)
# A1KHV0
```

## apfid.js

There is a port of APFID module to JS, for now it can be found in the same repository
(not distributed via pip with python module, not distributed via npm or other JS module manager, at least for now).