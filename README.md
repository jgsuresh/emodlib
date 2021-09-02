## emodlib

Porting key modules from [EMOD](https://github.com/InstituteforDiseaseModeling/EMOD) version 2.20 --
C++ libraries of demographic, immunological, and behavioral algorithms for disease modeling -- 
to lightweight decoupled libraries with Python bindings to facilitate uptake in other applications.

### dependencies

### Linux
```
apt-get -y install python cmake python-pybind11 libboost
```

#### Mac

```
brew install python
brew install cmake
brew install pybind11
brew install boost
```
* https://developer.apple.com/xcode/ (clang compiler included in XCode command line tools?)

#### Windows
* https://www.python.org/downloads/
* https://cmake.org/download/ (add to system path during MSI install steps so it's callable from commandline?)
* https://pybind11.readthedocs.io/en/stable/installing.html#include-with-pypi (use global install option to ensure pybind11 is findable by cmake?)
* https://docs.microsoft.com/en-us/cpp/build/vscpp-step-0-installation?view=msvc-160 (MSVC C/C++ compiler if not already installed with Visual Studio)
* https://www.boost.org/doc/libs/1_63_0/more/getting_started/windows.html (for boost.preprocessor enum handling)

### building C++ source

```
cmake -H. -Bbuild
cmake --build build -- -j3
```

### installing python package

```
cd python
python setup.py develop
```
