#!/bin/bash

# Define cofactor list (regex pattern)
COFACTORS="0WN|0LI|4MK|AXI|B49|BAX|DB8|EMH|EUI|FMM|IRE|LEV|LQQ|MI1|NIL|STI|VGH|YY3|END"

# Initialize output summary file
> prot.pdb

# Get base directory path
base_dir=$(pwd)

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
RESET='\033[0m'

# Header
echo -e "${CYAN}=== Cofactor Removal Summary ===${RESET}"

# Loop over directories
for dir in */; do
    dir=${dir%/}  # remove trailing slash
    full_path="$base_dir/$dir"
    prot_file="$dir/prot.pdb"
    new_prot_file="$dir/prot_new.pdb"
    prot_backupfile="$dir/run_kin_prot.pdb"

    # Skip if backup already exists
    if [[ -f "$prot_backupfile" ]]; then
        echo -e "${YELLOW}[SKIP] Already exists: $prot_backupfile${RESET}"
        continue
    fi

    # Skip if prot.pdb doesn't exist
    if [[ ! -f "$prot_file" ]]; then
        echo -e "${YELLOW}[MISSING] File not found: $prot_file${RESET}"
        continue
    fi

    # Backup original prot.pdb
    cp "$prot_file" "$prot_backupfile"
    echo -e "${GREEN}[OK] Backed up → $prot_backupfile${RESET}"

    # Extract only lines that start with HETATM and match cofactors
    match_lines=$(egrep "^HETATM.*($COFACTORS)" "$prot_backupfile")

    if [[ -n "$match_lines" ]]; then
        # Extract unique residue names (3-letter codes at cols 18-20)
        cofactors=$(echo "$match_lines" | awk '{print substr($0,18,3)}' | sort | uniq)

        # Join cofactors into a comma-separated list
        cof_list=$(echo "$cofactors" | paste -sd ", ")

        # Output to console
        echo -e "${CYAN}Removed [ ${cof_list} ] from → ${full_path}${RESET}"

        # Append raw cofactor lines to master file
        echo "$match_lines" >> prot.pdb
    fi

    # Write filtered file (remove only HETATM cofactors)
    egrep -v "^HETATM.*($COFACTORS)" "$prot_backupfile" > "$new_prot_file"
done

# Footer
echo -e "\n${CYAN}=====================================${RESET}"
echo -e "${GREEN}✅ Done.${RESET}"
echo -e "• Removed cofactor lines saved to: ${YELLOW}prot_new.pdb${RESET}"

