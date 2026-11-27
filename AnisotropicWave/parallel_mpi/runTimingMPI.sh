#!/usr/bin/env bash
# runTimingMPI.sh: timing of the serial and MPI codes (const coefficients, poly)
#   speedup S = T_serial / T_p and efficiency E = S / p, from the fastest of REPEAT runs
#   (the median is also reported, to show the spread)
#
# Usage:
#   ./runTimingMPI.sh
#   MACHINE=workstation NP="1 2 4 8 16" ./runTimingMPI.sh
#   NX=1024 TFINAL=0.1 REPEAT=5 ./runTimingMPI.sh
#   MPIRUN="mpirun" ./runTimingMPI.sh          MPI without --timeout (e.g. MPICH)
#   COEFF=VAR TFINAL=0.01 ./runTimingMPI.sh   variable coefficients (smaller dt: use a smaller tFinal)
# Output: printed and saved to reference/timing_mpi_<MACHINE>.txt (COEFF=VAR: timing_mpi_var_<MACHINE>.txt)

set -u
cd "$(dirname "$0")"

NX=${NX:-3000}   # avoid nx near a power of two (cache conflicts, about 3x slower at 2048)
TFINAL=${TFINAL:-0.1}
NP=${NP:-"1 2 4 5 8 10"}
REPEAT=${REPEAT:-3}
MACHINE=${MACHINE:-mac}
MPIRUN=${MPIRUN:-"mpirun --timeout 300"}
# Same compiler for serial and MPI: macOS Homebrew GCC, Linux system GCC
if [ "$(uname -s)" = Darwin ]; then CXX=${CXX:-g++-16}; else CXX=${CXX:-g++}; fi
export OMPI_CXX=$CXX MPICH_CXX=$CXX     # compiler used by mpicxx (Open MPI / MPICH)
COEFF=${COEFF:-CONS}                    # CONS | VAR
DEFS="-DCOEFF=$COEFF -DSOLUTION=POLY"
REF=../reference
if [ "$COEFF" = CONS ]; then caseName="const poly"; OUTFILE="$REF/timing_mpi_${MACHINE}.txt"
else caseName="var poly"; OUTFILE="$REF/timing_mpi_var_${MACHINE}.txt"; fi

BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT

# Machine description (no host name)
if [ "$(uname -s)" = Darwin ]; then
    cpu="$(sysctl -n machdep.cpu.brand_string) ($(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null)P+$(sysctl -n hw.perflevel1.physicalcpu 2>/dev/null)E cores)"
else
    cpu="$(lscpu 2>/dev/null | sed -n 's/^Model name: *//p' | head -1) ($(nproc) cores)"
fi
compiler="$($CXX --version | head -1)"
mpiver="$(mpirun --version 2>/dev/null | head -1)"
commit="$(git rev-parse --short HEAD 2>/dev/null)$(git diff --quiet -- .. 2>/dev/null || echo '-modified')"

# Build serial and MPI with the same compiler and options
$CXX -O3 $DEFS -o "$BUILD/serial" ../serial/wave2dSerial.C || { echo "serial build failed"; exit 1; }
mpicxx -O3 $DEFS -o "$BUILD/mpi" wave2dMPI.C || { echo "MPI build failed"; exit 1; }

# cpuTime of one run (from the MATLAB-style output line; for MPI the max over ranks)
cpuTime() { "$@" 2>/dev/null | grep -E '^[A-Za-z]+: t=' | sed -E 's/.*cpuTime= *([^ (]+)\(s\).*/\1/'; }
numSteps=$("$BUILD/serial" -nx=$NX -tFinal=$TFINAL | grep '^BASELINE:' | sed -E 's/.*numSteps=([0-9]+).*/\1/')

# min and median of REPEAT runs: prints "min median" (or "nan nan" if a run failed)
timeRuns() {
    for r in $(seq $REPEAT); do cpuTime "$@"; done | sort -g | awk -v n0=$REPEAT '
        { v[NR] = $1 }
        END { n = NR; if (n < n0) { print "nan nan"; exit }
              med = (n % 2) ? v[(n+1)/2] : (v[n/2] + v[n/2+1]) / 2; printf "%.4e %.4e\n", v[1], med }'
}

{
    echo "# MPI timing: $caseName, nx=$NX, tFinal=$TFINAL, numSteps=$numSteps, best of $REPEAT runs"
    echo "# machine=$MACHINE  cpu=$cpu"
    echo "# compiler=$compiler (via mpicxx, $mpiver)  flags=-O3  commit=$commit  date=$(date +%Y-%m-%d)"
    echo "# S_ser, E_ser: relative to the serial code;  S_1, E_1: relative to mpi with np=1"
    printf "%-8s %-8s %-11s %-11s %-7s %-7s %-7s %-7s\n" code np "min(s)" "median(s)" S_ser E_ser S_1 E_1

    read Ts Tsmed < <(timeRuns "$BUILD/serial" -nx=$NX -tFinal=$TFINAL)
    printf "%-8s %-8s %-11.3e %-11.3e %-7.2f %-7.2f %-7s %-7s\n" serial 1 "$Ts" "$Tsmed" 1 1 - -
    T1=""
    for p in $NP; do
        read Tp Tpmed < <(timeRuns $MPIRUN -np $p "$BUILD/mpi" -nx=$NX -tFinal=$TFINAL)
        if [ "$Tp" = nan ]; then
            printf "%-8s %-8s run failed (crash or timeout)\n" mpi "$p"; continue
        fi
        [ "$p" = 1 ] && [ -z "$T1" ] && T1=$Tp     # time with 1 thread / process
        awk -v p="$p" -v ts="$Ts" -v tp="$Tp" -v tm="$Tpmed" -v t1="$T1" 'BEGIN {
            s = ts / tp
            printf "%-8s %-8s %-11.3e %-11.3e %-7.2f %-7.2f ", "mpi", p, tp, tm, s, s / p
            if (t1 == "") printf "%-7s %-7s\n", "-", "-"; else { s1 = t1 / tp; printf "%-7.2f %-7.2f\n", s1, s1 / p } }'
    done
} | tee "$OUTFILE"

echo "Saved $OUTFILE"
