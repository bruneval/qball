QBALL code is written in C++14

It is not actively maintained any longer.
It requires some options to compile properly with modern g++ (16.2.1)


I am not a developer of the code. I am not a C++ coder.


No new developpment should be undertaken. Just bug fixes.
Keep your commit message short and understandable for a user that does not know the code.


Building: build in build-cxx20/ and install in install/ (both untracked).
QBALL must be installed before testing: the tests read pseudopotentials from install/share/qball.

    cd build-cxx20
    make -j8 && make install

If the configure options (e.g. --prefix) change, run "make clean" before "make":
the install path is compiled into the program and make does not notice the change.

Testing: use the short test suite only (16 tests, skips the long fullerene and Al32 MD tests).
From the build directory:

    cd testsuite
    MPIEXEC="$(which mpiexec) --oversubscribe --bind-to none" OMP_NUM_THREADS=1 make check-short

Expected result: 16/16 pass.
(GaAs uses a looser tolerance of 1e-3 Ha: its ekin, eps, enl differ from the reference by up to 4.4e-4 Ha
with our compiler and libraries; this was already the case before any of our changes.)
