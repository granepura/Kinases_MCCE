#!/bin/bash

# Define cofactor list (regex pattern)
COFACTORS="0WN|0LI|4MK|AXI|B49|BAX|DB8|EMH|EUI|FMM|IRE|LEV|LQQ|MI1|NIL|STI|VGH|YY3|END"

# Initialize output summary file
> step2_out.pdb

# Get base directory path
base_dir=$(pwd)

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
RESET='\033[0m'

# Open logging (both console + file)
exec > >(tee -a rm_cofs_step2_out.log) 2>&1

# Header
echo -e "${CYAN}=== Cofactor Removal Summary ===${RESET}"

# Loop over directories
for dir in */; do
    dir=${dir%/}  # remove trailing slash
    full_path="$base_dir/$dir"
    step2out_file="$dir/step2_out.pdb"
    step2out_backupfile="$dir/run_kin_step2_out.pdb"

    # Skip if backup already exists
    if [[ -f "$step2out_backupfile" ]]; then
        echo -e "${YELLOW}[SKIP] Already exists: $step2out_backupfile${RESET}"
        continue
    fi

    # Skip if step2_out.pdb doesn't exist
    if [[ ! -f "$step2out_file" ]]; then
        echo -e "${YELLOW}[MISSING] File not found: $step2out_file${RESET}"
        continue
    fi

    # Backup original step2_out.pdb
    cp "$step2out_file" "$step2out_backupfile"
    echo -e "${GREEN}[OK] Backed up → $step2out_backupfile${RESET}"

    # Extract matching lines
    match_lines=$(egrep "$COFACTORS" "$step2out_backupfile")

    if [[ -n "$match_lines" ]]; then
        # Extract unique residue names (3-letter codes)
        cofactors=$(echo "$match_lines" | awk '{print substr($0,18,3)}' | sort | uniq)

        # Join cofactors into a comma-separated list
        cof_list=$(echo "$cofactors" | paste -sd ", ")

        # Output to console
        echo -e "${CYAN}Removed [ ${cof_list} ] from → ${full_path}${RESET}"

        # Append raw cofactor lines to master file
        echo "$match_lines" >> step2_out.pdb
    fi

    # Write filtered file (without cofactors)
    egrep -v "$COFACTORS" "$step2out_backupfile" > "$step2out_file"
done

# Footer
echo -e "\n${CYAN}=====================================${RESET}"
echo -e "${GREEN}✅ Done.${RESET}"
echo -e "• Removed cofactor lines saved to: ${YELLOW}step2_out.pdb${RESET}"
echo -e "• Full console log saved to: ${YELLOW}rm_cofs_step2_out.log${RESET}"

