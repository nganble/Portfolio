#!/usr/bin/env bash
# runBaselineSerial.sh: baseline for the serial code (wave2dSerial.C)
#   5 cases with an exact solution, Nx = Ny = 20, 40, 80, 160 on (-1,1)^2, tFinal = 0.5
#
# Usage:
#   ./runBaselineSerial.sh          compare with reference/serial_*.txt and reference/matlab_*.txt
#   ./runBaselineSerial.sh save     write reference/serial_*.txt
#   CXX=clang++ ./runBaselineSerial.sh   use another compiler

set -u
cd "$(dirname "$0")"

MODE=${1:-compare}
CXX=${CXX:-g++}
CXXFLAGS="-O3"
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
        printf "max|maxErr-ref| = %8.2e  %s", dmax, (ok ? "PASS" : "FAIL")
        if (bad)  printf " (numSteps differ)"
        if (miss) printf " (nx missing in reference)"
        printf "\n"
        exit (ok ? 0 : 1)
    }' "$2" "$1"
}

allPass=1
summary=""
for c in $CASES; do
    COEFF=${c%%:*}; SOL=${c##*:}
    cn=$(coeffName "$COEFF"); mn=$(msName "$SOL")
    exe="$BUILD/w_${COEFF}_${SOL}"
    out="$BUILD/serial_${cn}_${mn}.txt"
    echo
    echo "===== coeffs=$cn, ms=$mn ====="

    if ! $CXX $CXXFLAGS -DCOEFF=$COEFF -DSOLUTION=$SOL -o "$exe" wave2dSerial.C; then
        echo "BUILD FAILED"; allPass=0; summary+="$cn $mn BUILD-FAILED"$'\n'; continue
    fi

    # Run all grids: MATLAB-style line plus ratio, BASELINE lines to $out
    : > "$out"
    prev=""
    for nx in $NX_LIST; do
        res=$("$exe" -nx=$nx -tFinal=$TFINAL)
        line=$(echo "$res" | grep -E '^[A-Za-z]+: t=')
        base=$(echo "$res" | grep '^BASELINE:')
        echo "$base" >> "$out"
        err=$(echo "$base" | sed -E 's/.*maxErr=([^ ]+).*/\1/')
        if [ -n "$prev" ]; then
            echo "$line $(awk -v a="$prev" -v b="$err" 'BEGIN{printf "ratio=%4.2f", a/b}')"
        else
            echo "$line"
        fi
        prev=$err
    done

    if [ "$MODE" = save ]; then
        cp "$out" "$REF/serial_${cn}_${mn}.txt"
        echo "Saved $REF/serial_${cn}_${mn}.txt"
        summary+="$cn $mn saved"$'\n'
    else
        status="PASS"
        for ref in serial matlab; do
            refFile="$REF/${ref}_${cn}_${mn}.txt"
            if [ -f "$refFile" ]; then
                printf "%-7s ref: " "$ref"; compare "$out" "$refFile" || status="FAIL"
            else
                echo "$ref ref: $refFile not found (run './runBaselineSerial.sh save' first)"; status="FAIL"
            fi
        done
        [ "$status" = PASS ] || allPass=0
        summary+="$cn $mn $status"$'\n'
    fi
done

echo
echo "===== Summary ($MODE) ====="
printf "%s" "$summary" | awk '{printf "%d  %-6s %-10s %s\n", NR, $1, $2, $3}'
if [ "$MODE" != save ]; then
    if [ $allPass = 1 ]; then echo "Overall: PASS"; else echo "Overall: FAIL"; exit 1; fi
fi
