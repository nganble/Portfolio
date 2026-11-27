#!/usr/bin/env bash
# runBaseline.sh: baseline for the serial and MPI Poisson solvers (red-black SOR)
#   3 cases on (0,1)^2, N = 16, 32, 64, 128, omega = 2.0 (optimal), tol = 1e-10:
#     poly     : TEST=EXACT_SOLUTION SOLUTION=POLY      exact (max error ~ tol)
#     gaussian : TEST=EXACT_SOLUTION SOLUTION=GAUSSIAN  second order (error ratio ~ 4)
#     compute  : TEST=COMPUTE_SOLUTION (no exact solution, residual only)
#   Checks: serial = saved reference; MPI (np = 1, 2, 3, 4) = serial;
#           poly max error < 1e-10; gaussian ratio on the two finest grids > 3.5
#
# Usage:
#   ./runBaseline.sh          compare with ../reference/serial_<case>.txt
#   ./runBaseline.sh save     write ../reference/serial_<case>.txt
#   NP="1 2 4" ./runBaseline.sh

set -u
cd "$(dirname "$0")"

MODE=${1:-check}
NLIST="16 32 64 128"
NP=${NP:-"1 2 3 4"}
ARGS_REST="2.0 1e-10 100000"           # omega tol max_iter
REF=../reference
CC=${CC:-gcc}
MPICC=${MPICC:-mpicc}
MPIRUN=${MPIRUN:-"mpirun --oversubscribe"}

cases=(poly gaussian compute)
flags=("-DTEST=EXACT_SOLUTION -DSOLUTION=POLY"
       "-DTEST=EXACT_SOLUTION -DSOLUTION=GAUSSIAN"
       "-DTEST=COMPUTE_SOLUTION")

BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT
mkdir -p "$REF"

# one line per grid: N iterations residual [error]
summary() {
	local N=$1; shift
	"$@" $N $ARGS_REST | awk -v N=$N '
		/after [0-9]+ iterations/ { for (i=1;i<=NF;i++) if ($i=="after") it=$(i+1) }
		/Max residual:/           { r=$3; sub(",","",r) }
		/Max error:/              { e=$3 }
		END { printf "%-5s %-6s %-9s %s\n", N, it, r, e }'
}

nfail=0
pass() { if [ "$1" = 0 ]; then echo "PASS"; else echo "FAIL"; nfail=$((nfail+1)); fi; }

for k in 0 1 2; do
	c=${cases[$k]}
	$CC    -O3 -DUSE_POWER9=0 ${flags[$k]} poisson2d_serial.c -o "$BUILD/serial_$c" -lm 2>/dev/null || { echo "build failed: serial $c"; exit 1; }
	$MPICC -O3 -DUSE_POWER9=0 ${flags[$k]} poisson2d_mpi.c    -o "$BUILD/mpi_$c"    -lm 2>/dev/null || { echo "build failed: mpi $c"; exit 1; }

	ser="$BUILD/serial_$c.txt"
	{ echo "# N     iter   residual  error"; for N in $NLIST; do summary $N "$BUILD/serial_$c"; done; } > "$ser"
	echo "===== $c ====="
	cat "$ser"

	if [ "$MODE" = save ]; then
		cp "$ser" "$REF/serial_$c.txt"; echo "Saved $REF/serial_$c.txt"; continue
	fi

	# serial vs reference
	printf "%-28s" "serial vs reference:"
	if [ -f "$REF/serial_$c.txt" ]; then diff -q "$ser" "$REF/serial_$c.txt" >/dev/null; pass $?
	else echo "no reference (run with save)"; nfail=$((nfail+1)); fi

	# MPI vs serial
	for np in $NP; do
		out="$BUILD/mpi_${c}_$np.txt"
		{ echo "# N     iter   residual  error"; for N in $NLIST; do summary $N $MPIRUN -np $np "$BUILD/mpi_$c"; done; } > "$out"
		printf "%-28s" "mpi np=$np vs serial:"
		diff -q "$ser" "$out" >/dev/null; pass $?
	done

	# accuracy
	if [ $c = poly ]; then
		printf "%-28s" "poly: max error < 1e-10:"
		awk 'NR>1 && $4+0 >= 1e-10 {bad=1} END {exit bad}' "$ser"; pass $?
	elif [ $c = gaussian ]; then
		ratio=$(awk 'NR>1 {e[n++]=$4} END {printf "%.2f", e[n-2]/e[n-1]}' "$ser")
		printf "%-28s" "gaussian: ratio=$ratio > 3.5:"
		awk -v r=$ratio 'BEGIN {exit !(r > 3.5)}'; pass $?
	fi
done

[ "$MODE" = save ] && exit 0
echo "======================"
if [ $nfail = 0 ]; then echo "All checks PASS"; else echo "$nfail check(s) FAIL"; fi
exit $nfail
