This is HPIPM, a high-performance interior-point method solver for dense, optimal control- and tree-structured convex quadratic programs. It provides efficient implementations of dense and structure-exploiting algorithms to solve small to medium scale problems arising in model predictive control and embedded optimization in general and it relies on the high-performance linear algebra package BLASFEO.

--------------------------------------------------

## Getting Started:
The best way to get started with HPIPM is to check out the examples in `/hpipm/examples/`.
HPIPM can be directly used from `C`, but there are also interfaces to Python and Matlab.
Depending on which level you want to use HPIPM, check out the following section below.
The QP notation used in HPIPM can be found in the `doc` folder.

### C
In order to run the C examples in `/hpipm/examples/C/` follow the steps below:
1) clone BLASFEO on your machine: `git clone https://github.com/giaf/blasfeo.git` 
2) from the BLASFEO root folder, run `make static_library & sudo make install_static`
3) from the HPIPM root folder, run `make static_library & make examples`
4) cd to /hpipm/examples/c/ and run getting_started.out to solve a simple OCP-structured QP.

### Python
If you would like to try out the Python interface, check out the examples in `/hpipm/examples/python/` after going through the following steps:
1) clone BLASFEO on your machine: `git clone https://github.com/giaf/blasfeo.git`
2) from the BLASFEO root folder, run `make shared_library -j4 & sudo make install_shared`
3) from the HPIPM root folder, run `make shared_library -j4 & sudo make install_shared`
4) make sure that the location of the installed shared libraries is known to the system by running `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/blasfeo/lib:/opt/hpipm/lib`. If you would like 
to avoid running this command whenever opening a new shell. You can add the commands above to your .bashrc.
5) cd to /hpipm/interfaces/python/hpipm_python and run `pip3 install .`
6) cd to /hpipm/examples/python and run `python3 getting_started.py` to solve a simple OCP-structured QP.


### MATLAB
Since the HPIPM MATLAB interface is based on the HPIPM Python interface, please first follow the steps in the previous paragraph, "Python".
In order to be able to run HPIPM from Matlab the following steps are additionally required:
1) Make sure that the location of the shared libraries `libblasfeo.so` and `libhpipm.so` is known to the system, see step 4) for Python.
2) Check the paths are defined in the terminal using: `echo $LD_LIBRARY_PATH`
3) Start Matlab from a terminal in which the `LD_LIBRARY_PATH` is defined correctly
4) Make sure that hpipm_python has been installed for the Python version used by MATLAB.
You can check which version is used by running the Matlab command `py.sys.path`.
Use the Matlab command `pyversion <path_to_the_python_interpreter>` to change the Python interpreter used by Matlab.
5) Add `hpipm/interfaces/matlab/hpipm_matlab` to your MATLAB path.
6) Get started using HPIPM from Matlab by looking at the examples in `/hpipm/examples/matlab/`.


--------------------------------------------------

References:

- G. Frison, H.H. B. Sørensen, B. Dammann, and J.B. Jørgensen. High-performance
small-scale solvers for linear model predictive control. In
IEEE European Control Conference, pages 128–133. IEEE, 2014 - <https://ieeexplore.ieee.org/document/6981589/>

- BLASFEO: Basic Linear Algebra Subroutines For Embedded Optimization
G. Frison, D. Kouzoupis, T. Sartor, A. Zanelli, M. Diehl
ACM Transactions on Mathematical Software (TOMS) (2018) - <https://arxiv.org/abs/1704.02457>

- <https://github.com/giaf/blasfeo>

--------------------------------------------------

Notes:

- HPIPM is released under the GPL+CE (Gnu General Public License version 3.0 with Classpath Exception).
This license has been chosen to avoid the issues that the LGPL license implies in case of static linking.
Unformally speaking, the GPL+CE license allows the HPIPM library to be linked to both statically and dynamically to produce an executable, and the resulting executable can be released under terms of your choice.
Please refer to the LICENSE.txt file for the legally binding statement of the license.
