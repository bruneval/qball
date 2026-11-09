QBALL code is written in C++14

It is not actively maintained any longer.
It requires some options to compile properly with modern g++ (16.2.1)


I am not a developer of the code. I am not a C++ coder.


No new developpment should be undertaken. Just bug fixes.
Keep your commit message short and understandable for a user that does not know the code.


Testing: use the short test suite only (13 tests, skips the long fullerene and Al32 MD tests).
From the build directory:

    cd testsuite
    MPIEXEC="$(which mpiexec) --oversubscribe --bind-to none" OMP_NUM_THREADS=1 make check-short
