#!/bin/bash
set -euo pipefail

# 1. Configuration
# Auto-detect cores, default to 4 if detection fails
THREADS=$(nproc 2>/dev/null || echo 4)

# Define directories
# We create them first so we can resolve absolute paths immediately
mkdir -p "FoldX.out1/repaired" "FoldX.out1/out" "FoldX.out1/temp_workspaces" "FoldX.out1/fragments"

# Resolve ABSOLUTE paths to prevent issues when changing directories in subshells
REPAIRED=$(readlink -f "FoldX.out1/repaired")
OUTDIR=$(readlink -f "FoldX.out1/out")
TMP_INPUT=$(readlink -f "FoldX.out1/temp_workspaces")
FRAGMENTS=$(readlink -f "FoldX.out1/fragments")
INITIAL_PDB_DIR="models"

# Export variables for worker functions
export REPAIRED OUTDIR TMP_INPUT FRAGMENTS

# Clean up previous temporary data
rm -rf "$TMP_INPUT"/* "$FRAGMENTS"/*

# 2. Define Worker Function
process_pdb() {
    local pdb="$1"
    
    # Resolve input path
    local abs_pdb
    abs_pdb=$(readlink -f "$pdb")
    local pdb_file
    pdb_file=$(basename "$abs_pdb")
    local base="${pdb_file%.pdb}"

    # Sanitization
    local safe_base
    safe_base=$(echo "$base" | tr '.' '_' | tr ':' '_')
    
    echo ">>> Processing $base (Worker PID: $$)"

    # ISOLATION: Create a unique workspace for this job
    # FoldX is not thread-safe if multiple instances run in the same dir 
    local job_dir="$TMP_INPUT/$safe_base"
    mkdir -p "$job_dir"
    
    # Symlink input pdb to job dir
    ln -sf "$abs_pdb" "$job_dir/${safe_base}.pdb"

    # -- REPAIR --
    # We switch into the job_dir so FoldX writes its temp files there
    pushd "$job_dir" > /dev/null
    
    FoldX --command=RepairPDB \
        --pdb "${safe_base}.pdb" \
        --output-dir "$REPAIRED" > /dev/null 2>&1

    # Verify repair
    local repaired_pdb_name="${safe_base}_Repair.pdb"
    local abs_repaired_path="$REPAIRED/$repaired_pdb_name"
    
    if [[ ! -f "$abs_repaired_path" ]]; then
        echo "Error: Repair failed for $base"
        popd > /dev/null
        return 1
    fi

    # -- STABILITY --
    # Point to the global REPAIRED directory for input
    FoldX --command=Stability \
        --pdb-dir "$REPAIRED" \
        --pdb "$repaired_pdb_name" \
        --output-dir "$OUTDIR" > /dev/null 2>&1

    popd > /dev/null

    # -- PARSING --
    local fxout="${OUTDIR}/${safe_base}_Repair_0_ST.fxout"

    if [[ -f "$fxout" ]]; then
        local dG
        # PARSING FIX: Extract 2nd column (Total Energy) directly using awk.
        # Format: [PDB_PATH] [TOTAL_ENERGY] ...
        dG=$(awk '{print $2}' "$fxout")

        # Basic validation to ensure we captured a value
        if [[ -n "$dG" ]]; then
            # Write to a unique fragment file to avoid write conflicts
            echo "$base,$dG" > "$FRAGMENTS/${safe_base}.csv"
        else
            echo "Warning: Parsed empty energy for $base"
        fi
    else
        echo "Warning: No stability output found for $base"
    fi
    
    # Cleanup workspace
    rm -rf "$job_dir"
}

export -f process_pdb

# 3. Execution
echo ">>> Starting parallel processing on $THREADS cores..."

# find PDBs -> xargs parallel execution -> process_pdb
find "$INITIAL_PDB_DIR" -name "chimera*.pdb" -print0 | \
    xargs -0 -n 1 -P "$THREADS" -I {} bash -c 'process_pdb "$@"' _ {}

# 4. Aggregation
echo ">>> Aggregating results..."
echo "model,dG_kcal_mol" > FoldX.out1/foldx_raw.csv

if ls "$FRAGMENTS"/*.csv 1> /dev/null 2>&1; then
    cat "$FRAGMENTS"/*.csv >> FoldX.out1/foldx_raw.csv
else
    echo "Error: No results generated."
    exit 1
fi

# Cleanup output structure
rm -rf "$TMP_INPUT" "$FRAGMENTS"

###########################################
# 5. Ranking
###########################################
echo ">>> Ranking results"
CSV="FoldX.out1/stability_ranked.csv"

awk -F',' 'NR>1{print $0}' FoldX.out1/foldx_raw.csv | sort -t',' -k2,2n > FoldX.out1/sorted.tmp

if [[ ! -s FoldX.out1/sorted.tmp ]]; then
     echo "Error: No valid energies found."
     exit 1
fi

best=$(head -n1 FoldX.out1/sorted.tmp | cut -d',' -f2)

{
    echo "model,dG_kcal_mol,ddG_vs_best"
    awk -F',' -v best="$best" '{
        ddG = $2 - best
        printf "%s,%.4f,%.4f\n", $1, $2, ddG
    }' FoldX.out1/sorted.tmp
} > "$CSV"

rm FoldX.out1/sorted.tmp
echo ">>> Output written to $CSV"
