![PyPI - License](https://img.shields.io/pypi/l/pydefect?color=blue)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pydefect)
[![CircleCI](https://circleci.com/gh/kumagai-group/pydefect/tree/master.svg?style=shield)](https://circleci.com/gh/kumagai-group/pydefect/tree/master)

pydefect
=========
Pydefect is a robust, open-source Python library for point-defect calculations in non-metallic solids 
based on first-principle calculations with the [VASP](https://www.vasp.at) code.

**Note: Units used in pydefect are eV for energy and angstrom for length following the vasp convention.**

Installation instructions
---------------------------------------------------------
1. Requirements
  - Python 3.7 or higher
  - vise
  - pymatgen
  - see requirements.txt for others
  

2. Latest stable version is released at PyPI repository, so one can download 
it using `pip install pydefect`.

3. After cloning the repository, it is possible to install `pydefect` using the python package manager pip.
To do so, run this command in the directory containing setup.py:

`pip install ./`

Sometimes errors will be given if a specific version of a package is not 
installed. In this case try installing the exact versions that have been 
tested for use with `pydefect` with the following command:

`pip install -r requirements.txt ./`

This will install the package versions listed in the requirements.txt file.
To prevent interference with other programmes, it is advised that a package 
management system like 
[conda](https://docs.conda.io/projects/conda/en/latest/index.html) is used. 
This allows you to install dependencies in a particular environment so that 
they can be managed and recorded more easily. The commands above can then be 
executed in a conda environment (after installing pip in that environment).

For more information on how dependencies are managed in this branch see this [blog post](https://medium.com/@boscacci/why-and-how-to-make-a-requirements-txt-f329c685181e).

Executing this software
--------------------------

1. Command line method
  - execute ***pydefect -h*** for detailed description of available options

2. Usage as a module
  - pydefect can be imported as a python module

More information in the online manual at: https://kumagai-group.github.io/pydefect/

Files and directories included in pydefect distribution
--------------------------------------------------------
~~~
  README                 : introduction
  LICENSE                : the MIT license 
  setup.py               : installation script
  requirements.txt       : list of required packages

  /pydefect/analyzer     : analysis tools for point-defect calculations
  /pydefect/chem_pot_diag: tools for calculating and drawing the chemical potential diagram
  /pydefect/cli          : command line interfaces
  /pydefect/corrections  : energy and eigenvalue correction related modules
  /pydefect/database     : database related to atoms and symmetries
  /pydefect/input_maker  : tools for generating VASP input files
  /pydefect/tests        : test files used mainly for unittests
  /pydefect/util         : useful tools 
~~~~

License
-----------------------
Python code is licensed under the MIT License.

Development notes
-------------------
### Bugs, requests and questions
Please use the [Issue Tracker](https://github.com/kumagai-group/pydefect/issues) to report bugs, request features.

### Code contributions
Although pydefect is free to use, we sinserely appreciate if you help us to improve this library. 
The simplest but most valuable contribution is to send the feature requests and bug reports.

Please report any bugs and issues at PyDefect's [Github Issues page](https://github.com/oba-group/pydefect).
Please use the ["Fork and Pull"](https://guides.github.com/activities/forking/) workflow to make contributions and stick as closely as possible to the following:

- Code style follows [PEP8](http://www.python.org/dev/peps/pep-0008) and [Google's writing style](https://google.github.io/styleguide/pyguide.html).
- Add unittests wherever possible including scripts for command line interfaces.

### Tests
Run the tests using `pytest pydefect`.
We also use integrated testing on Github via circleCI.

Citing vise
-------------
If pydefect has been used in your research, please cite our Github repo
until we have a publication!

Contact info
--------------
Yu Kumagai<br>
yuuukuma@gmail.co.jp

Tokyo Institute of Technology (Japan)

