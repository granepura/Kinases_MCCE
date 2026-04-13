#!/bin/bash

#for d in */; do cp "$d/step2_out.pdb" "$d/run_kin_step2_out.pdb"; done

grep HETATM ../run_kin/1XKK/step2_out.pdb | grep -Ev 'NME|ACE' > 1XKK/step2_out.pdb
grep HETATM ../run_kin/2EUF/step2_out.pdb | grep -Ev 'NME|ACE' > 2EUF/step2_out.pdb
grep HETATM ../run_kin/2HYY/step2_out.pdb | grep -Ev 'NME|ACE' > 2HYY/step2_out.pdb
grep HETATM ../run_kin/2ITO/step2_out.pdb | grep -Ev 'NME|ACE' > 2ITO/step2_out.pdb
grep HETATM ../run_kin/2ITY/step2_out.pdb | grep -Ev 'NME|ACE' > 2ITY/step2_out.pdb
grep HETATM ../run_kin/2ITZ/step2_out.pdb | grep -Ev 'NME|ACE' > 2ITZ/step2_out.pdb
grep HETATM ../run_kin/2WGJ/step2_out.pdb | grep -Ev 'NME|ACE' > 2WGJ/step2_out.pdb
grep HETATM ../run_kin/2XP2/step2_out.pdb | grep -Ev 'NME|ACE' > 2XP2/step2_out.pdb
grep HETATM ../run_kin/2YFX/step2_out.pdb | grep -Ev 'NME|ACE' > 2YFX/step2_out.pdb
grep HETATM ../run_kin/3AOX/step2_out.pdb | grep -Ev 'NME|ACE' > 3AOX/step2_out.pdb
grep HETATM ../run_kin/3CS9/step2_out.pdb | grep -Ev 'NME|ACE' > 3CS9/step2_out.pdb
grep HETATM ../run_kin/3LXK/step2_out.pdb | grep -Ev 'NME|ACE' > 3LXK/step2_out.pdb
grep HETATM ../run_kin/3PYY/step2_out.pdb | grep -Ev 'NME|ACE' > 3PYY/step2_out.pdb
grep HETATM ../run_kin/3UE4/step2_out.pdb | grep -Ev 'NME|ACE' > 3UE4/step2_out.pdb
grep HETATM ../run_kin/3UG2/step2_out.pdb | grep -Ev 'NME|ACE' > 3UG2/step2_out.pdb
grep HETATM ../run_kin/3WZD/step2_out.pdb | grep -Ev 'NME|ACE' > 3WZD/step2_out.pdb
grep HETATM ../run_kin/3WZE/step2_out.pdb | grep -Ev 'NME|ACE' > 3WZE/step2_out.pdb
grep HETATM ../run_kin/3ZOS/step2_out.pdb | grep -Ev 'NME|ACE' > 3ZOS/step2_out.pdb
grep HETATM ../run_kin/4AG8/step2_out.pdb | grep -Ev 'NME|ACE' > 4AG8/step2_out.pdb
grep HETATM ../run_kin/4AGC/step2_out.pdb | grep -Ev 'NME|ACE' > 4AGC/step2_out.pdb
grep HETATM ../run_kin/4AGD/step2_out.pdb | grep -Ev 'NME|ACE' > 4AGD/step2_out.pdb
grep HETATM ../run_kin/4AN2/step2_out.pdb | grep -Ev 'NME|ACE' > 4AN2/step2_out.pdb
grep HETATM ../run_kin/4ANQ/step2_out.pdb | grep -Ev 'NME|ACE' > 4ANQ/step2_out.pdb
grep HETATM ../run_kin/4ANS/step2_out.pdb | grep -Ev 'NME|ACE' > 4ANS/step2_out.pdb
grep HETATM ../run_kin/4ASD/step2_out.pdb | grep -Ev 'NME|ACE' > 4ASD/step2_out.pdb
grep HETATM ../run_kin/4G5J/step2_out.pdb | grep -Ev 'NME|ACE' > 4G5J/step2_out.pdb
grep HETATM ../run_kin/4G5P/step2_out.pdb | grep -Ev 'NME|ACE' > 4G5P/step2_out.pdb
grep HETATM ../run_kin/4I22/step2_out.pdb | grep -Ev 'NME|ACE' > 4I22/step2_out.pdb
grep HETATM ../run_kin/4LMN/step2_out.pdb | grep -Ev 'NME|ACE' > 4LMN/step2_out.pdb
grep HETATM ../run_kin/4MKC/step2_out.pdb | grep -Ev 'NME|ACE' > 4MKC/step2_out.pdb
grep HETATM ../run_kin/4WKQ/step2_out.pdb | grep -Ev 'NME|ACE' > 4WKQ/step2_out.pdb
grep HETATM ../run_kin/4ZAU/step2_out.pdb | grep -Ev 'NME|ACE' > 4ZAU/step2_out.pdb
grep HETATM ../run_kin/5AAA/step2_out.pdb | grep -Ev 'NME|ACE' > 5AAA/step2_out.pdb
grep HETATM ../run_kin/5AAB/step2_out.pdb | grep -Ev 'NME|ACE' > 5AAB/step2_out.pdb
grep HETATM ../run_kin/5AAC/step2_out.pdb | grep -Ev 'NME|ACE' > 5AAC/step2_out.pdb
grep HETATM ../run_kin/5L2I/step2_out.pdb | grep -Ev 'NME|ACE' > 5L2I/step2_out.pdb
grep HETATM ../run_kin/5MO4/step2_out.pdb | grep -Ev 'NME|ACE' > 5MO4/step2_out.pdb


