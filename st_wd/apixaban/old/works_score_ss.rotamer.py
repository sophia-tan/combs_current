import optparse # for option sorting
from pyrosetta import *
init()

init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!

def pose_scoring(pose, display_residues = []):

    # create a PyMOLMover for exporting structures directly to PyMOL
    pymover = PyMOLMover()


    # b. this method returns a ScoreFunction with its weights set based on
    #    files stored in the database/scoring/weights (.wts files)
    full_scorefxn = create_score_function('ref2015')

    #### c. this method sets the weights based on 'talaris2013.wts' and then
    ####    corrects them based on 'docking.wts_patch'
    ###ws_patch_scorefxn = create_score_function('talaris2013', 'docking')  # create_score_function_ws_patch('talaris2013', 'docking')

    #### d. this method returns a ScoreFunction with its weights set by loading
    ####    weights from 'talaris2013' followed by an adjustment by setting
    ####    weights from 'docking.wts_patch'
    ###patch_scorefxn = create_score_function('talaris2013')
    ###patch_scorefxn.apply_patch_from_file('docking')

    #### e. here an empty ScoreFunction is created and the weights are set manually
    ###scorefxn = ScoreFunction()
    ###scorefxn.set_weight(core.scoring.fa_atr, 0.800)    # full-atom attractive score
    ###scorefxn.set_weight(core.scoring.fa_rep, 0.440)    # full-atom repulsive score
    ###scorefxn.set_weight(core.scoring.fa_sol, 0.750)    # full-atom solvation score
    ###scorefxn.set_weight(core.scoring.fa_intra_rep, 0.004)    # f.a. intraresidue rep. score
    ###scorefxn.set_weight(core.scoring.fa_elec, 0.700)    # full-atom electronic score
    ###scorefxn.set_weight(core.scoring.pro_close, 1.000)    # proline closure
    ###scorefxn.set_weight(core.scoring.hbond_sr_bb, 1.170)    # short-range hbonding
    ###scorefxn.set_weight(core.scoring.hbond_lr_bb, 1.170)    # long-range hbonding
    ###scorefxn.set_weight(core.scoring.hbond_bb_sc, 1.170)    # backbone-sidechain hbonding
    ###scorefxn.set_weight(core.scoring.hbond_sc, 1.100)    # sidechain-sidechain hbonding
    ###scorefxn.set_weight(core.scoring.dslf_fa13, 1.000)    # disulfide full-atom score
    ###scorefxn.set_weight(core.scoring.rama, 0.200)    # ramachandran score
    ###scorefxn.set_weight(core.scoring.omega, 0.500)    # omega torsion score
    ###scorefxn.set_weight(core.scoring.fa_dun, 0.560)    # fullatom Dunbrack rotamer score
    ###scorefxn.set_weight(core.scoring.p_aa_pp, 0.320)
    ###scorefxn.set_weight(core.scoring.ref, 1.000)    # reference identity score

    # ScoreFunction a, b, and e above have the same weights and thus return
    #    the same score for an input pose.  Likewise, c and d should return the
    #    same scores.
    # 2. output the ScoreFunction evaluations
    #ws_patch_scorefxn(pose)    # to prevent verbose output on the next line
    print( '='*80 )
    print( 'ScoreFunction b:', full_scorefxn(pose) )
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
    print( '='*80 )
    print( 'Loaded from' , pose.pdb_info().name() )
    print( pose.total_residue() , 'residues' )
    print( 'Total Rosetta Score:' , pose_score )
    full_scorefxn.show(pose)
    # this object is contained in PyRosetta v2.0 and above
    pymover.apply(pose)
    pymover.send_energy(pose)

    # 10. output information on the requested residues
    for i in display_residues:
        print( '='*80 )
        print( 'Pose numbered Residue' , i )
        print( 'Total Residue Score:' , residue_energies[i-1] )
        print( 'Score Breakdown:\n' + '-'*45 )
        # loop over the weights, extract the scores from the matrix
        for w in range(len(weights)):
            print( '\t' + rosetta.core.scoring.name_from_score_type(weights[w]).ljust(20) + ':\t' ,\
                residue_weighted_energies_matrix[w][i-1] )
        print( '-'*45 )
        # print the hydrogen bond information
        print( 'Hydrogen bonds involving Residue ' + str(i) + ':' )
        print( hbond_dictionary[i][:-1] )
    print( '='*80 )

################################################################################
# COMMANDLINE COMPATIBILITY

# everything below is added to provide commandline usage,
#   the available options are specified below
# this method:
#    1. defines the available options
#    2. loads in the commandline or default values
#    3. calls pose_scoring with these values

# parser object for managing input options
# all defaults are for the example using "test_in.pdb" with reduced
#    cycles/jobs to provide results quickly
parser = optparse.OptionParser()
parser.add_option('--pdb_filename', dest = 'pdb_filename',
    #default = '1wdn_0001.pdb',    # default example PDB
    help = 'the PDB file containing the loop to remodel')
parser.add_option('--residues', dest = 'residues',
    default = '',    # default to the median residue number
    help = 'the (pose numbered) residues to inspect carefully')
(options,args) = parser.parse_args()

# PDB file option
pdb_filename = options.pdb_filename
# create a pose from the desired PDB file
# -create an empty Pose object
pose = Pose()
# -load the data from pdb_file into the pose
pose_from_file(pose, pdb_filename)
# default to the median residue number
residues = options.residues
if not options.residues:
    residues = [int(pose.total_residue()/2)]
elif options.residues == 'all':
    # accept the word 'all' in place of a residue list
    residues = range(1, pose.total_residue() + 1)
else:
    # please provide the residues of interest as , delimited
    residues = [int(r) for r in options.residues.split(',')]

pose_scoring(pose, residues)

