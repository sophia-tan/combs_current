# what's the prevalence of pro-tyr motifs, like in 2df3?
cutoff=0.5
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.apps import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
#from clusterScoring import *
#from cluster_for_sophia import *
sys.path.append('/home/gpu/Sophia/combs/st_wd/Scoring_Function')
from ScoringInteractions import *
                
parsed = pr.parsePDB('pymolintegrin.pdb')
db = pkl.load(open('/home/gpu/Sophia/combs/st_wd/Lookups/AAi_freq/AAi_database_lookups.pkl','rb'))

ifgresn = 'VAL'
vdmresn = 'THR'
ifgatoms = 'CG1 CB CG2'.split(' ')
vdmatoms = 'CG2 CB OG1'.split(' ')
ifginfo = ['A', '760']
vdminfo = ['B', '656']

#vdmresn = 'VAL'
#ifgresn = 'THR'
#vdmatoms = 'CG2 CB CG1'.split(' ')
#ifgatoms = 'OG1 CB CG2'.split(' ')
#vdminfo = ['A', '760']
#ifginfo = ['B', '656']

# filter  for only vdmresn vdms of ifgresn with ifgatoms
# and vdmatoms directly invoved in interactions
num_all_vdms, num_direct, lookupdf, filtered_df = filter_contact(ifgresn, vdmresn, ifgatoms, vdmatoms)
print('num all vdms', num_all_vdms)
print('num direct', len(lookupdf))
print('vdmresn %age', db[1][vdmresn])
print('ifgresn %age', db[1][ifgresn])

query = []
for atom in ifgatoms:
    query.append(
        parsed.select(
            'chain {} and resnum {} and name {}'.format(
                ifginfo[0],
                ifginfo[1],
                atom)).getCoords()[0])

for atom in vdmatoms:
    query.append(
        parsed.select(
            'chain {} and resnum {} and name {}'.format(
                vdminfo[0],
                vdminfo[1],
                atom)).getCoords()[0])
query = np.array(query)
lookupcoords = pkl.load(open(
    '/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/coords_of_{}.pkl'.format(
    constants.AAname[ifgresn]), 'rb'))
    
#lookupcoords = lookupcoords[:50] # delete

ifglists = flip(ifgatoms, ifgresn)
vdmlists = flip(vdmatoms, vdmresn)
rmsds = []
num_atoms = len(query)
coords_ls = [
    item for item in lookupcoords if item[0] in lookupdf.index]
lookupatoms_to_clus = []
lookupatoms_to_clus.append(query) # first element to see how big its cluster is
counter = 0 # to keep count of how many pdbs are being output
for item in coords_ls:
    if len(item) == 3:
        compare_rmsds = []
        ifg_vdm_ind = []
        for ifg_ind, ifgls in enumerate(ifglists):
            for vdm_ind, vdmls in enumerate(vdmlists):
                lookupatoms = get_order_of_atoms(
                    item, ifgresn, vdmresn, ifgls, vdmls)
                moved, transf = pr.superpose(lookupatoms, query)
                temp_rmsd = pr.calcRMSD(moved, query)
                compare_rmsds.append(temp_rmsd)
                ifg_vdm_ind.append([moved, temp_rmsd])
        # item[0] is df index
        if min(compare_rmsds)>2:
            print(compare_rmsds)
        rmsds.append([item[0], min(compare_rmsds)])
        # get index of which one had min rmsd
        for which_ind, each in enumerate(ifg_vdm_ind):
            if each[1] == min(compare_rmsds):
                lookupatoms_to_clus.append(each[0])

                ############################################################ 
                #                output pdb if low rmsd 
                ############################################################ 
                #if each[1] < cutoff and counter < 30 and which_ind==0: 
                #    # this is to ensure rmsd is below cutoff when not flipped
                #    # bc don't want to take care of that in prody to output pdb
                #    row = lookupdf.loc[item[0]]
                #    db_dir = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'
                #    par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
    
                #    ifgchid, ifgresnum = row['chid_ifg'], row['resnum_ifg']
                #    vdmchid, vdmresnum = row['chid_vdm'], row['resnum_vdm']
                #    printout = copy.deepcopy(par)
                #    printout = printout.select(
                #        '(chain {} and resnum {}) or (chain {} and resnum {})'.format(
                #        ifgchid,ifgresnum,vdmchid,vdmresnum))
                #    printout.select('chain {} and resnum {}'.format(ifgchid,ifgresnum)).setChids('Y')
                #    printout.select('chain {} and resnum {}'.format(vdmchid,vdmresnum)).setChids('X')
                #    printout.select('all').setResnums(10)
                #    printout_interactamer = []
                #    integrin_interactamer = []
                #    try: # skip the ones that have segment ids.
                #        for atom in ifgatoms:
                #            integrin_interactamer.append(parsed.select('chain {} and resnum {} and name {}'.format(ifginfo[0],ifginfo[1],atom)))
                #            printout_interactamer.append(printout.select('chain Y and resnum 10 and name {}'.format(atom)))
                #        for atom in vdmatoms:
                #            integrin_interactamer.append(parsed.select('chain {} and resnum {} and name {}'.format(vdminfo[0],vdminfo[1],atom)))
                #            printout_interactamer.append(printout.select('chain X and resnum 10 and name {}'.format(atom)))
                #        integrin_interactamer_prody = []

                #        
                #        integrin_interactamer = sum(integrin_interactamer[1:], integrin_interactamer[0])
                #        printout_interactamer = sum(printout_interactamer[1:], printout_interactamer[0])
                #        try:
                #            assert len(integrin_interactamer) == len(printout_interactamer)

                #            interact_res = printout.select('(chain X and resnum 10) or (chain Y and resnum 10)')
                #            interactamer_transf = pr.applyTransformation(transf, printout_interactamer)
                #            
                #        except:
                #            pass
                #    except:
                #        traceback.print_exc()
                #        pass

    else:
        rmsds.append([int(item[0]), 100000])

pkl.dump(rmsds, open('bleb.pkl','wb'))
# count how many coords are within rmsd cutoff
rmsds = [i for i in rmsds if i[1] < cutoff]
allrmsds = [i[1] for i in rmsds]
print(len(rmsds), 'number of matches')
