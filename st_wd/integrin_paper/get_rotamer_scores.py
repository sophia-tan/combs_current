from PyrosettaScores import *
from residues_integrin import *

residues = []
for k,v in integrin_res.items():
    residues.append(k[1:])

pyrosetta_scores('cleanpymolintegrin.pdb',residues)
