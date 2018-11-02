# GADES++

## 1.  About

GADES++ is a similarity measure to determine the relatedness between two entities patients based on the evolution of the tumour stages, mutations, and ECOG Performance Status, among others characteristics. 

## 2. License

MIT LICENSE

## 3. Requirements

* Python 3.X

* py_stringmatching

* [UMLS::Similarity](http://www.d.umn.edu/~tpederse/umls-similarity.html)

## Running one sample

GADES++ can be used as a standalone program o library. 
The file "patients1.tsv" contains three entities to compare.
To run GADES++ execute:

`>./gades_plus_plus.py test/patients1.tsv`
