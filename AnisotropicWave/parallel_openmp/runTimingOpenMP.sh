#!/usr/bin/env bash
# runTimingOpenMP.sh: timing of the serial and OpenMP codes (const coefficients, poly)
#   speedup S = T_serial / T_p and efficiency E = S / p, from the fastest of REPEAT runs
#   (the median is also reported, to show the spread)
#
# Usage:
#   ./runTimingOpenMP.sh
#   MACHINE=workstation THREADS="1 2 4 8 16" ./runTimingOpenMP.sh
#   NX=1024 TFINAL=0.1 REPEAT=5 ./runTimingOpenMP.sh
# Output: printed and saved to reference/timing_openmp_<MACHINE>.txt

set -u
cd "$(dirname "$0")"

NX=${NX:-2048}
TFINAL=${TFINAL:-0.1}
THREADS=${THREADS:-"1 2 4 5 8 10"}
REPEAT=${REPEAT:-3}
MACHINE=${MACHINE:-mac}
# Same compiler choice as the Makefile: macOS Homebrew GCC, Linux system GCC
if [ "$(uname -s)" = Darwin ]; then CXX=${CXX:-g++-16}; else CXX=${CXX:-g++}; fi
DEFS="-DCOEFF=CONS -DSOLUTION=POLY"
REF=../reference
OUTFILE="$REF/timing_openmp_${MACHINE}.txt"

BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT

# Machine description (no host name)
if [ "$(uname -s)" = Darwin ]; then
    cpu="$(sysctl -n machdep.cpu.brand_string) ($(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null)P+$(sysctl -n hw.perflevel1.physicalcpu 2>/dev/null)E cores)"
else
    cpu="$(lscpu 2>/dev/null | sed -n 's/^Model name: *//p' | head -1) ($(nproc) cores)"
fi
compiler="$($CXX --version | head -1)"
commit="$(git rev-parse --short HEAD 2>/dev/null)$(git diff --quiet -- .. 2>/dev/null || echo '-modified')"

# Build serial and OpenMP with the same compiler and options
$CXX -O3 $DEFS -o "$BUILD/serial" ../serial/wave2dSerial.C || { echo "serial build failed"; exit 1; }
$CXX -O3 -fopenmp $DEFS -o "$BUILD/openmp" wave2dOpenMP.C || { echo "OpenMP build failed"; exit 1; }

# cpuTime of one run (from the MATLAB-style output line)
cpuTime() { "$@" | grep -E '^[A-Za-z]+: t=' | sed -E 's/.*cpuTime= *([^ (]+)\(s\).*/\1/'; }
numSteps=$("$BUILD/serial" -nx=$NX -tFinal=$TFINAL | grep '^BASELINE:' | sed -E 's/.*numSteps=([0-9]+).*/\1/')

# min and median of REPEAT runs: prints "min median"
timeRuns() {
    for r in $(seq $REPEAT); do cpuTime "$@"; done | sort -g | awk '
        { v[NR] = $1 }
        END { n = NR; med = (n % 2) ? v[(n+1)/2] : (v[n/2] + v[n/2+1]) / 2; printf "%.4e %.4e\n", v[1], med }'
}

{
    echo "# OpenMP timing: const poly, nx=$NX, tFinal=$TFINAL, numSteps=$numSteps, best of $REPEAT runs"
    echo "# machine=$MACHINE  cpu=$cpu"
    echo "# compiler=$compiler  flags=-O3  commit=$commit  date=$(date +%Y-%m-%d)"
    echo "# S_ser, E_ser: relative to the serial code;  S_1, E_1: relative to openmp with threads=1"
    printf "%-8s %-8s %-11s %-11s %-7s %-7s %-7s %-7s\n" code threads "min(s)" "median(s)" S_ser E_ser S_1 E_1

    read Ts Tsmed < <(timeRuns "$BUILD/serial" -nx=$NX -tFinal=$TFINAL)
    printf "%-8s %-8s %-11.3e %-11.3e %-7.2f %-7.2f %-7s %-7s\n" serial 1 "$Ts" "$Tsmed" 1 1 - -
    T1=""
    for p in $THREADS; do
        read Tp Tpmed < <(timeRuns "$BUILD/openmp" -nx=$NX -tFinal=$TFINAL -numThreads=$p)
        [ "$p" = 1 ] && [ -z "$T1" ] && T1=$Tp     # time with 1 thread / process
        awk -v p="$p" -v ts="$Ts" -v tp="$Tp" -v tm="$Tpmed" -v t1="$T1" 'BEGIN {
            s = ts / tp
            printf "%-8s %-8s %-11.3e %-11.3e %-7.2f %-7.2f ", "openmp", p, tp, tm, s, s / p
            if (t1 == "") printf "%-7s %-7s\n", "-", "-"; else { s1 = t1 / tp; printf "%-7.2f %-7.2f\n", s1, s1 / p } }'
    done
} | tee "$OUTFILE"

echo "Saved $OUTFILE"
