import optparse # for option sorting
import numpy as np
from pyrosetta import *

init(extra_options = "-constant_seed -ignore_unrecognized_res -mute basic -mute core -mute protocols")
# WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!

def pose_scoring(pose, display_residues = []):

    # create a PyMOLMover for exporting structures directly to PyMOL
    pymover = PyMOLMover()

    full_scorefxn = create_score_function('ref2015')
    pose_score = full_scorefxn(pose)

    # 3. obtain the pose Energies object and all the residue total scores
    energies = pose.energies()
    residue_energies = [energies.residue_total_energy(i)
        for i in range(1, pose.total_residue() + 1)]

    # 4. obtain the non-zero weights of the ScoreFunction, active ScoreTypes
    weights = [rosetta.core.scoring.ScoreType(s)
        for s in range(1, int(rosetta.core.scoring.end_of_score_type_enumeration) + 1)
        if full_scorefxn.weights()[rosetta.core.scoring.ScoreType(s)]]
    # 5. obtain all the pose energies using the weights list
    # Energies.residue_total_energies returns an EMapVector of the unweighted
    #    score values, here they are multiplied by their weights
    #    remember when performing individual investigation, these are the raw
    #    unweighted score!
    residue_weighted_energies_matrix = [
        [energies.residue_total_energies(i)[w] * full_scorefxn.weights()[w]
        for i in range(1, pose.total_residue() + 1)]
        for w in weights]
    return residue_weighted_energies_matrix


################################################################################
def pyrosetta_scores(design_pdb, scaffold_pdb, vdmresnum):
    # create a pose from the desired PDB file
    # -create an empty Pose object
    pose = Pose()
    # -load the data from pdb_file into the pose
    pose_from_file(pose, design_pdb) # score rotamer on lig_vdm pdb
    design_matrix = pose_scoring(pose)
    pose = Pose()
    pose_from_file(pose, scaffold_pdb) # score phi/psi on whole scaffold bc need flanking residues
    scaffold_matrix = pose_scoring(pose)
    # index 18 of matrix is rama_prepro. index 14 is dun
    return design_matrix[14][0], scaffold_matrix[18][vdmresnum]

