## emodlib

Porting key modules from [EMOD](https://github.com/InstituteforDiseaseModeling/EMOD) version 2.20 --
C++ libraries of demographic, immunological, and behavioral algorithms for disease modeling -- 
to lightweight decoupled libraries with Python bindings to facilitate uptake in other applications.

### dependencies

brew install python
brew install cmake
brew install pybind11

### building C++ source

cmake -H. -Bbuild
cmake --build build -- -j3

### installing python package

cd python
python setup.py develop
