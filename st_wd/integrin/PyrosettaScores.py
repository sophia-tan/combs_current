import optparse # for option sorting
import numpy as np
from pyrosetta import *
init()

init(extra_options = "-constant_seed -ignore_unrecognized_res")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!

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

    # Unfortunately, hydrogen bonding scores are NOT stored in the structure
    #    returned by Energies.residue_total_energies
    # 6. hydrogen bonding information must be extracted separately
    pose_hbonds = rosetta.core.scoring.hbonds.HBondSet()
    rosetta.core.scoring.hbonds.fill_hbond_set( pose , False , pose_hbonds )

    # 7. create a dictionary with the pose residue numbers as keys and the
    #    residue hydrogen bonding information as values
    # hydrogen bonding information is stored as test in the form:
    # (donor residue) (donor atom) => (acceptor residue) (accecptor atom) |score
    hbond_dictionary = {}
    for residue in range(1, pose.total_residue() + 1):
        hbond_text = ''
        for hbond in range(1, pose_hbonds.nhbonds() + 1):
            hbond = pose_hbonds.hbond(hbond)
            acceptor_residue = hbond.acc_res()
            donor_residue = hbond.don_res()
            if residue == acceptor_residue or residue == donor_residue:
                hbond_text += str(donor_residue).ljust(4) + ' ' + \
                    str(pose.residue(donor_residue).atom_name(\
                        hbond.don_hatm() )).strip().ljust(4) + \
                    ' => ' + str(acceptor_residue).ljust(4) + ' ' + \
                    str(pose.residue(acceptor_residue).atom_name(\
                        hbond.acc_atm() )).strip().ljust(4) + \
                    ' |score: ' + str(hbond.energy()) + '\n'
        hbond_dictionary[residue] = hbond_text


    # 9. output the pose information
	# the information is not expressed sequentially as it is produced because
	#    several PyRosetta objects and methods output intermediate information
	#    to screen, this would produce and unattractive output
    #######print( '='*80 )
    #######print( 'Loaded from' , pose.pdb_info().name() )
    #######print( pose.total_residue() , 'residues' )
    #######print( 'Total Rosetta Score:' , pose_score )
    full_scorefxn.show(pose)
    # this object is contained in PyRosetta v2.0 and above
    pymover.apply(pose)
    pymover.send_energy(pose)
    # 10. output information on the requested residues
    # loop over the weights, extract the scores from the matrix
    scoresdict = {}
    for w in [13,14,18]: # these are just omega, dun, rama scores
        scoresdict[rosetta.core.scoring.name_from_score_type(weights[w])] = \
            np.round(residue_weighted_energies_matrix[w][0],2)
    return scoresdict

################################################################################
def pyrosetta_scores(design_pdb, vdmresnum):
    # create a pose from the desired PDB file
    # -create an empty Pose object
    pose = Pose()
    # -load the data from pdb_file into the pose
    pose_from_file(pose, design_pdb)
    scoresdict = pose_scoring(pose, [vdmresnum])
    return scoresdict['rama_prepro'], scoresdict['fa_dun']

