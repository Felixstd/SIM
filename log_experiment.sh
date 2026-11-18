#!/usr/bin/env bash
# Append one experiment (one namelist) as a row to Experiments.md.
#
# Usage:  ./log_experiment.sh <namelist_file> ["short description"]
#
# nx, ny and Deltax are read from src/parameter.h; every name in PARAMS is
# looked up in the namelist first, then in parameter.h.
set -euo pipefail

# ---- settings: edit these -------------------------------------------------
PARAMS=(Rheology regularization solver
        Deltat gamma_nl NLmax OLmax Pstar C e_ratio Tens mu_0 mu_inf delta_c phi_prime uniaxial friction dilatation non_normal repl_press)
PARAM_FILE="src/parameter.h"
INPUT_FILE=$2      # file whose "exp version" line gives the name
OUT="Experiments.md"
# ---------------------------------------------------------------------------

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <namelist_file> [\"description\"]" >&2
    exit 1
fi
NML=$1
DESC=${3:-}
[[ -f $NML ]]        || { echo "Namelist not found: $NML" >&2; exit 1; }
[[ -f $PARAM_FILE ]] || { echo "Parameter file not found: $PARAM_FILE" >&2; exit 1; }
[[ -f $INPUT_FILE ]] || { echo "Input file not found: $INPUT_FILE" >&2; exit 1; }

# get_value KEY FILE -> prints the value of "KEY = value" (case-insensitive),
# ignoring Fortran comments (!), quotes, trailing commas and '&'.
get_value() {
    awk -v key="$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')" '
    {
        line = $0
        sub(/!.*/, "", line)                       # drop comment
        if (match(tolower(line), "^[ \t]*" key "[ \t]*=")) {
            val = substr(line, RLENGTH + 1)
            sub(/[ \t,&]*$/, "", val)              # trailing , & spaces
            sub(/^[ \t]+/, "", val)                # leading spaces
            gsub(/["\047]/, "", val)               # quotes
            print val
            exit
        }
    }' "$2"
}

# lookup KEY -> namelist first, then parameter.h, else "—"
lookup() {
    local v
    v=$(get_value "$1" "$NML")
    [[ -z $v ]] && v=$(get_value "$1" "$PARAM_FILE")
    if [[ -z $v ]]; then
        echo "Warning: '$1' not found in $NML or $PARAM_FILE" >&2
        v="—"
    fi
    printf '%s' "$v"
}

NX=$(lookup nx)
NY=$(lookup ny)
DX_RAW=$(lookup Deltax)
# Fortran 1d03 -> 1e03, then metres -> km
DX_KM=$(printf '%s' "$DX_RAW" | tr 'dD' 'eE' | awk '{ printf "%g", $1 / 1000 }')

# Experiment name = first field of the line labelled "exp version"
NAME=$(awk '/exp version/ { print $1; exit }' "$INPUT_FILE")
if [[ -z $NAME ]]; then
    echo "Error: no 'exp version' line found in $INPUT_FILE" >&2
    exit 1
fi
DATE=$(date +%Y-%m-%d)

# Build header and row
header="| Experiment | Date | nx | ny | dx [km]"
sep="|---|---|---|---|---"
row="| $NAME | $DATE | $NX | $NY | $DX_KM"
for p in "${PARAMS[@]}"; do
    header+=" | $p"
    sep+="|---"
    row+=" | $(lookup "$p")"
done
header+=" | Description |"
sep+="|---|"
row+=" | $DESC |"

if [[ ! -f $OUT ]]; then
    printf '# Experiments\n\n%s\n%s\n' "$header" "$sep" > "$OUT"
else
    # Header of the last (most recent) table in the file
    last_header=$(grep '^| Experiment' "$OUT" | tail -n 1 || true)
    if [[ $last_header != "$header" ]]; then
        # Columns changed: start a new table below the old one
        printf '\n%s\n%s\n' "$header" "$sep" >> "$OUT"
        echo "Columns changed: started a new table in $OUT"
    fi
fi

printf '%s\n' "$row" >> "$OUT"
echo "Added '$NAME' to $OUT"