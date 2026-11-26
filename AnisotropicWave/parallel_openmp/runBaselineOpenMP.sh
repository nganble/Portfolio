#!/usr/bin/env bash
# runBaselineOpenMP.sh: baseline for the OpenMP code (wave2dOpenMP.C)
#   5 cases with an exact solution, Nx = Ny = 20, 40, 80, 160 on (-1,1)^2, tFinal = 0.5,
#   run with several thread counts and compared with reference/serial_*.txt and reference/matlab_*.txt
#
# Usage:
#   ./runBaselineOpenMP.sh
#   THREADS="1 2 4 8 16" ./runBaselineOpenMP.sh   other thread counts
#   CXX=g++ ./runBaselineOpenMP.sh                other compiler

set -u
cd "$(dirname "$0")"

# Same compiler choice as the Makefile: macOS Homebrew GCC, Linux system GCC
if [ "$(uname -s)" = Darwin ]; then CXX=${CXX:-g++-16}; else CXX=${CXX:-g++}; fi
CXXFLAGS="-O3 -fopenmp"
THREADS=${THREADS:-"1 2 4 8"}
REF=../reference
TOL=1e-12
TFINAL=0.5
NX_LIST="20 40 80 160"
CASES="CONS:POLY CONS:TRIG CONS:PULSE CONS:PULSETRIG VAR:POLY"   # COEFF:SOLUTION

BUILD=$(mktemp -d)
trap 'rm -rf "$BUILD"' EXIT

# Names used in the reference files (same as MATLAB coeffs / ms)
coeffName() { case $1 in CONS) echo const;; VAR) echo var;; esac; }
msName()    { case $1 in POLY) echo poly;; TRIG) echo trig;; PULSE) echo pulse;; PULSETRIG) echo trigPulse;; esac; }

# compare OUT REF: same numSteps for every nx and max|maxErr - ref| < TOL
compare() {
    awk -v tol="$TOL" '
    function val(s, key,   i, n, a) {
        n = split(s, a, " ")
        for (i = 1; i <= n; i++) if (index(a[i], key "=") == 1) return substr(a[i], length(key) + 2)
        return ""
    }
    FNR == NR { if ($1 == "BASELINE:") { nx = val($0, "nx"); rs[nx] = val($0, "numSteps"); re[nx] = val($0, "maxErr") }; next }
    $1 == "BASELINE:" {
        nx = val($0, "nx"); n++
        if (!(nx in re)) { miss = 1 }
        else {
            if (val($0, "numSteps") != rs[nx]) bad = 1
            d = val($0, "maxErr") - re[nx]; if (d < 0) d = -d
            if (d > dmax) dmax = d
        }
    }
    END {
        ok = (n > 0 && !miss && !bad && dmax < tol)
        printf "%8.2e %s", dmax, (ok ? "PASS" : "FAIL")
        if (bad)  printf " (numSteps differ)"
        if (miss) printf " (nx missing)"
        exit (ok ? 0 : 1)
    }' "$2" "$1"
}

allPass=1
summary=$(printf "%-17s" "case"; for t in $THREADS; do printf " %-5s" "$t"; done)$'\n'
for c in $CASES; do
    COEFF=${c%%:*}; SOL=${c##*:}
    cn=$(coeffName "$COEFF"); mn=$(msName "$SOL")
    exe="$BUILD/w_${COEFF}_${SOL}"
    echo
    echo "===== coeffs=$cn, ms=$mn (table: threads=${THREADS%% *}) ====="

    row=$(printf "%-17s" "$cn $mn")
    if ! $CXX $CXXFLAGS -DCOEFF=$COEFF -DSOLUTION=$SOL -o "$exe" wave2dOpenMP.C; then
        echo "BUILD FAILED"; allPass=0; summary+="$row BUILD-FAILED"$'\n'; continue
    fi

    first=""
    for t in $THREADS; do
        out="$BUILD/openmp_${cn}_${mn}_t$t.txt"
        : > "$out"
        prev=""
        for nx in $NX_LIST; do
            res=$("$exe" -nx=$nx -tFinal=$TFINAL -numThreads=$t)
            base=$(echo "$res" | grep '^BASELINE:')
            echo "$base" >> "$out"
            # MATLAB-style table (with ratio) only for the first thread count
            if [ -z "$first" ]; then
                line=$(echo "$res" | grep -E '^[A-Za-z]+: t=')
                err=$(echo "$base" | sed -E 's/.*maxErr=([^ ]+).*/\1/')
                if [ -n "$prev" ]; then
                    echo "$line $(awk -v a="$prev" -v b="$err" 'BEGIN{printf "ratio=%4.2f", a/b}')"
                else
                    echo "$line"
                fi
                prev=$err
            fi
        done

        status=PASS
        msg=$(printf "threads=%-3s" "$t")
        for ref in serial matlab; do
            refFile="$REF/${ref}_${cn}_${mn}.txt"
            if [ -f "$refFile" ]; then
                r=$(compare "$out" "$refFile") || status=FAIL
                msg+="  $ref ref: $r"
            else
                msg+="  $ref ref: not found"; status=FAIL
            fi
        done
        # Results must not depend on the number of threads
        if [ -z "$first" ]; then
            first=$out
        elif ! cmp -s "$first" "$out"; then
            msg+="  (differs from threads=${THREADS%% *})"; status=FAIL
        fi
        echo "$msg"
        [ "$status" = PASS ] || allPass=0
        row+=$(printf " %-5s" "$status")
    done
    summary+="$row"$'\n'
done

echo
echo "===== Summary (threads) ====="
printf "%s" "$summary"
if [ $allPass = 1 ]; then echo "Overall: PASS"; else echo "Overall: FAIL"; exit 1; fi
