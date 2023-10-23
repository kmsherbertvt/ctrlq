CtrlQ
=====
CtrlQ is an open-source tool designed to simulate a gate-free state preparation on a Transmon qubit device using analog control pulses. The analog control pulses can be variationally shaped to drive an initial state to a target state in the framework of ctrl-VQE. In molecular systems, ctrl-VQE can be used to drive the initial Hartree Fock state to the full configuration interaction (FCI) state with substantial pulse duration speedups as compared to a gate-based compilation. 

The control quantum program (CtrlQ) is written in python with bindings to C++
codes for highly efficient time-evolution of quantum systems either using an
ordinary-differential-equation or the Suzuki-Trotter expansion. Efficient
analytic gradients for pulse parameters is implemented which allows
optimization of thousands of pulse parameters with only about 2.5 times the
cost of an energy evaluation.

**Reference**
OR Meitei, BT Gard, GS Barron, DP Pappas, SE Economou, E Barnes, NJ Mayhall, Gate-free state preparation for fast variational quantum eigensolver simulations: ctrl-VQE
[arXiv:2008.04302](https://arxiv.org/abs/2008.04302)

## Installation
While mostly python, ctrlq delegates some of the numerically intensive work to C++.
Thus, you will need a C++ compiler and the `cmake` tool to install.
Detailed information for installation, particularly dependencies, is provided in the documentation but here's a quick step-by-step.

1. Get the source code from github:

       git clone --recursive https://github.com/kmsherbertvt/ctrlq.git

2. Install the Eigen library:

   Download from https://eigen.tuxfamily.org/index.php?title=Main_Page and extract to preferred location.

   (Kyle's choice: the `ctrlq` project directory itself)

   Set the `EIGEN_INCLUDE` environment variable to the Eigen path, relative to `ctrlq` project directory.
   (This is so that the `cmake` command knows where to find all the Eigen header files.)

   (eg. with Kyle's choice: `export EIGEN_INCLUDE="eigen-3.4.0"`)

2. Configure with cmake and compile

       cd ctrl
       mkdir build && cd build
       cmake ..
       make
      
3. Run test

       python -m unittest discover


## Documentation
Documentation is available online in html format [here](https://ctrlq.readthedocs.io)
and also can be found in the ``/ctrlq/doc`` directory. The documentation
contains detailed installation instruction including linking to optimized MKL
libraries and tutorials to run CtrlQ.




