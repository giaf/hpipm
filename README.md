This is HPIPM, a high-performance interior-point method solver for dense, optimal control- and tree-structured convex quadratic programs (QP) and quadratically-constrained quadratic programs (QCQP).
It provides efficient implementations of dense and structure-exploiting algorithms to solve small to medium scale problems arising in model predictive control and embedded optimization in general and it relies on the high-performance linear algebra package BLASFEO.

HPIPM (and BLASFEO, which is a dependency), comes with both `make` and `cmake` build systems.
The preferred one is `make`, which can be used to compile and run any library, interface and example in any language.
`make` is also used in the continuous integration `travis` scripts.
`cmake` can only be used to compile the libraries, while the interested user should compile interfaces and run examples by him/herself by taking inspiration form the commands in the various `Makefile`s.

--------------------------------------------------

## Getting Started:
The best way to get started with HPIPM is to check out the examples in `/hpipm/examples/`.
HPIPM can be directly used from `C`, but there are also interfaces to Python and Matlab.
Depending on which level you want to use HPIPM, check out the following section below.
The QP notation used in HPIPM can be found in the `doc` folder.

### C
In order to run the C examples in `/hpipm/examples/C/` follow the steps below:
1) Clone BLASFEO on your machine: `git clone https://github.com/giaf/blasfeo.git` 
2) From the BLASFEO root folder, run `make static_library -j 4 && sudo make install_static` (default installation folder: `/opt/blasfeo`; if a different one is chosen, `BLASFEO_PATH` in HPIPM's `Makefile.rule` should be updated accordingly)
3) From the HPIPM root folder, run `make static_library -j 4 && make examples`
4) In a terminal, navigate to `/hpipm/examples/c/` and run `getting_started.out` to solve a simple OCP-structured QP.

### MATLAB and Octave
#### Linux
The interface for Matlab and Octave is based on mex files.
1) Clone BLASFEO on your machine: `git clone https://github.com/giaf/blasfeo.git`
2) From the BLASFEO root folder, run `make shared_library -j 4 && sudo make install_shared`
3) From the HPIPM root folder, run `make shared_library -j 4 && sudo make install_shared`
4) In a terminal, navigate to the folder `hpipm/interfaces/matlab_octave`.
Set the needed environment flags by running `source env.sh` (you may need to change the `BLASFEO_MAIN_FOLDER`, or to make it equal to the `BLASFEO_PATH`) in that folder.
Compile the interface by running `make all -j 4` (for Octave), or `make compile_mex_with_matlab` (for Matlab).
5) In a terminal, navigate to the folder `hpipm/examples/matlab_octave`.
Set the needed environment flags by running `source env.sh` (you may need to change the `BLASFEO_MAIN_FOLDER`, or to make it equal to the `BLASFEO_PATH`) in that folder.
Run an instance of Matlab or Octave from the same terminal.
Get started by running the examples in that folder.

#### MATLAB on Windows
The interface for Matlab and Octave is based on mex files.
1) Clone BLASFEO on your machine: `git clone https://github.com/giaf/blasfeo.git`
2) Install [Microsoft Visual C++](https://visualstudio.microsoft.com/downloads/)
3) From the BLASFEO root folder, run
```
mkdir build
cd build
cmake ..
cmake --build .
```
Copy `blasfeo.lib` from `build/Debug/` to `lib/`.

4) From the HPIPM root folder, run
```
mkdir build
cd build
cmake ..
cmake --build .
```
(you may need to change the `BLASFEO_PATH` to the actual value matching your BLASFEO installation, e.g. by replacing the second last line above with `cmake -DBLASFEO_PATH=/path/to/your/blasfeo/installation/folder ..`).
Copy `hpipm.lib` from `build/Debug/` to `lib/`.

5) Open Matlab and navigate to the folder `hpipm/interfaces/matlab_octave`.
Set the needed environment flags by running `env.m` (you may need to change the `BLASFEO_MAIN_FOLDER`, or to make it equal to the `BLASFEO_PATH`) in that folder.
Compile the interface by running `compile_mex_all.m`.

(Make sure that Matlab is configured to compile mex files with the same compiler as used in point 2), as sometimes by default it makes use of MinGW64 instead.
If the HPIPM and BLASFEO libraries and the mex files are not all compiled with visual studio, or all compiled with MinGW64, errors about missing dependences from the runtime library happen.)

6) In Matlab, navigate to the folder `hpipm/examples/matlab_octave`.
Get started by running the examples in that folder. You may need to add folder `hpipm/interfaces/matlab_octave` to the Matlab path.

### Simulink
The QP model is read from the file `qp_data.c`, which can be generated using the C, matlab/octave or python interfaces.
1) Follow the steps 1)-4) for the MATLAB interface.
2) In a terminal, navigate to the folder `hpipm/examples/simulink`.
Run `make_sfun.m` to compile the S-function, and `load_paramaters.m` to load some parameters used in the simulink model (e.g. horizon length, number of inputs and states) form `qp_data.c`.
3) Open the simulink model `hpipm_simulink_getting_started.slx` and start the simulation. 

### Python
If you would like to try out the Python interface, check out the examples in `/hpipm/examples/python/` after going through the following steps:
1) Clone BLASFEO on your machine: `git clone https://github.com/giaf/blasfeo.git`
2) From the BLASFEO root folder, run `make shared_library -j 4 && sudo make install_shared`
3) From the HPIPM root folder, run `make shared_library -j 4 && sudo make install_shared`
4) In a terminal, navigate to `/hpipm/interfaces/python/hpipm_python` and run `pip install .` or  `pip3 install .` (depending on your python version).
5) In a terminal, navigate to `/hpipm/examples/python`.
Set the needed environment flags by running `source env.sh` (you may need to change the `BLASFEO_MAIN_FOLDER`, or to make it equal to the `BLASFEO_PATH`) in that folder.
Alternatively you can make sure yourself that the location of the installed shared libraries is known to the system e.g. by running `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/blasfeo/lib:/opt/hpipm/lib` (possibly after updating it to the chosen installation directories).
Finally, run `python example_qp_getting_started.py` or `python3 example_qp_getting_started.py` (depending on your python version) to solve a simple OCP-structured QP.

--------------------------------------------------

## Common mistakes:

### Matrix format
For performance reasons, HPIPM uses an internal matrix (and vector) format based on BLASFEO own matrix (and vector) structs.
The exact layout of the matrix elements in memory is implementation- and architecture-dependent, so it can't be just inspected directly.

However, in the C interface there are setter and getter routines to correctly populate these internal matrix (and vector) representations.
By default these setter and getter routines expect matrices to be passed in column-major order, that is the default order in e.g. Matlab and optimized linear algebra libraries such as BLAS and LAPACK.
There is also a versions of a few setter and getter routines that expect matrices to be passed in row-major order, and these routines are clearly marked using the wording `_rowmaj` at the end of the routine name.
These routines are basically identical to the default (aka column-major) ones with the only difference that the input matrix is transposed while being packed in / unpacked from the internal format.

Notice that the memory layout of the elements in the column-major matrix order corresponds to the layout of the transposed of the same matrix in the row-major matrix order, and the other way around.

Also notice that in C a matrix represented as an array of arrays is stored in row-major order.
Similarly, in Python a matrix represented using numpy as an array of arrays is also stored in row-major order.
On the other hand, in Matlab a dense matrix is stored in column-major order.
This can be source of confusion.

--------------------------------------------------

## References:

- G. Frison, M. Diehl.
*HPIPM: a high-performance quadratic programming framework for model predictive control*.
(2020)
(arXiv preprint <https://arxiv.org/abs/2003.02547>)

- G. Frison, J. Frey, F. Messerer, A. Zanelli, M. Diehl.
*Introducing the quadratically-constrained quadratic programming framework in HPIPM*.
(2022)
(arXiv preprint <https://arxiv.org/abs/2112.11872>)

- G. Frison, H.H. B. Sørensen, B. Dammann, and J.B. Jørgensen.
*High-performance small-scale solvers for linear model predictive control*.
In IEEE European Control Conference, pages 128–133. IEEE, 2014 - <https://ieeexplore.ieee.org/document/6981589/>

- G. Frison, D. Kouzoupis, T. Sartor, A. Zanelli, M. Diehl.
*BLASFEO: Basic Linear Algebra Subroutines For Embedded Optimization*.
ACM Transactions on Mathematical Software (TOMS) (2018)
(arXiv preprint <https://arxiv.org/abs/1704.02457>)

- <https://github.com/giaf/blasfeo>

--------------------------------------------------

## Notes:

- HPIPM relies on the high-performance linear algebra library BLASFEO.
BLASFEO provides several implementations optimized for different computer architectures, and it makes heavy use of assembly code.
If you get the error `Illegal instruction` at run time, you are probably using a BLASFEO version (`TARGET`) unsupported by your CPU.
