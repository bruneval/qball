#!/bin/sh
# Checks for testsuite/tddft/03_neon_sudden_anl (traveling pseudopotentials).
# Reads the QBALL output of a single Ne atom suddenly set moving at v = 0.5 a.u.
# along x (TD_dt = 0.025 a.u.) and prints one number:
#
#   anl N              A_x(N) - (P(N) - <p_x>(N)), should be 0
#   balance N1 N2      [E_KS + A.v](N2) - [E_KS + A.v](N1) - v (P(N2) - P(N1)), should be 0
#   ks_balance N1 N2   same without A.v: E_KS alone misses -v (A(N2) - A(N1))
#
# N = RT-TDDFT step number (1 = first step). All quantities of step N are
# printed together, for the same wavefunctions:
#   E_KS  <etotal> (the first one in the output is the ground state)
#   A     <anl_velocity> Ne Ax Ay Az </anl_velocity> (nonlocal part of the momentum)
#   <p>   canonical momentum = - total_electronic_current (electron charge -1)
#   P     physical momentum d<sum x>/dt = - d(dipole)/dt, central difference
#         of <electronic_dipole_density> (dipole computed from the density)
# Exits with an error if a needed quantity is missing (e.g. A not printed).

out=$1; mode=$2; n1=$3; n2=$4
awk -v mode="$mode" -v n1="$n1" -v n2="$n2" -v v=0.5 -v dt=0.025 '
function P(n) { if (n-1 < 1 || n+1 > nd) { bad = 1; return 0 } return -(D[n+1]-D[n-1])/(2*dt) }
function A(n) { if (n > na) { bad = 1; return 0 } return AX[n] }
function E(n) { if (n+1 > ne) { bad = 1; return 0 } return EK[n+1] }
function p(n) { if (n > nj) { bad = 1; return 0 } return -J[n] }
/<etotal>/                      { ne++; EK[ne] = $2 }
/total_electronic_current:/     { nj++; J[nj] = $2 }
/<electronic_dipole_density>/   { nd++; D[nd] = $2 }
/<anl_velocity>/                { na++; AX[na] = $3 }
END {
  bad = 0
  if (mode == "anl")             r = A(n1) - (P(n1) - p(n1))
  else if (mode == "balance")    r = (E(n2) + v*A(n2)) - (E(n1) + v*A(n1)) - v*(P(n2) - P(n1))
  else if (mode == "ks_balance") r = E(n2) - E(n1) - v*(P(n2) - P(n1))
  else bad = 1
  if (bad) exit 1
  printf "%.8f", r
}' "$out"
