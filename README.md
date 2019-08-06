This is HPIPM, a high-performance interior-point method solver for dense, optimal control- and tree-structured convex quadratic programs.
It provides efficient implementations of dense and structure-exploiting algorithms to solve small to medium scale problems arising in model predictive control and embedded optimization in general and it relies on the high-performance linear algebra package BLASFEO.

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
4) make sure that the location of the installed shared libraries is known to the system by running `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/blasfeo/lib:/opt/hpipm/lib`.
If you would like to avoid running this command whenever opening a new shell, you can add the commands above to your .bashrc.
5) cd to `/hpipm/interfaces/python/hpipm_python` and run `pip3 install .`
6) cd to `/hpipm/examples/python` and run `python3 getting_started.py` to solve a simple OCP-structured QP.

### MATLAB and Octave
The interface for Matlab and Octave is based on mex files.
1) Compile and install BLASFEO and HPIPM BLASFEO and HPIPM (see instruction for Python 1)-3)):
2) In a terminal, navigate to the folder `hpipm/interfaces/matlab_octave`, set the needed environment flags by running `source env.sh` in that folder and finally compile the interface by running the function `compile_mex_ocp_qp` in a Matlab or Octave instance started from the same terminal.
3) In a terminal, navigate to the folder `hpipm/examples/matlab_octave`.
Set the needed environment flags by running `source env.sh` in that folder.
Run an instance of Matlab or Octave from the same terminal.
Get started by running the examples in that folder.

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

- HPIPM relies on the high-performance linear algebra library BLASFEO.
BLASFEO provides several implementations optimized for different computer architectures, and it makes heavy use of assembly code.
If you get the error `Illegal instruciton` at running time, you are probably using a BLASFEO version (`TARGET`) unsupported by your CPU.
