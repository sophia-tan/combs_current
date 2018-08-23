# HAVE TO RUN IN THE ./inputs/ DIRECTORY
#!/bin/bash

for pdb in A756  A758  A772  A774  A852  A911  B1435  B1453  B1495  B1504  B1527  B1559  B1565
    do
        cp -r 3FCS 3FCS$pdb # copy dir to make one for this mutant 
        echo ""${pdb:1:5}" "${pdb:0:1}" PIKAA A " >> 3FCS$pdb/nataa_mutations.resfile # append each res to resfile
    done

