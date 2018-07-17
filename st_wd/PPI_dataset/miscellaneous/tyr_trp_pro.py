import pandas as pd, numpy as np
import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
import prody as pr, pickle as pkl
from clusterScoring import *

cutoff = 0.5

tyr_csv = '/home/gpu/Sophia/combs/st_wd/20180626_combed_csvs/{}/{}_vdm_pdb_info.csv'.format('tyrosine', 'tyrosine')
tyr_ifg_csv = '/home/gpu/Sophia/combs/st_wd/20180626_combed_csvs/{}/{}_ifg_pdb_info.csv'.format('tyrosine', 'tyrosine')

parsed = pr.parsePDB('2df3')

ifgresn = 'TYR'
vdmresn = 'TYR'
ifgatoms = 'N CA CB CG CD'.split(' ')
vdm1atoms = 'CG CD1 CD2 CE1 CE2 CZ'.split(' ')
vdm2atoms = 'CG CD1 NE1 CE2 CD2 CE3 CZ3 CH2 CZ2'.split(' ')
ifginfo = ['A', '51']
vdm1info = ['A', '26']
vdm2info = ['A', '132']


query = []
for atom in ifgatoms:
    query.append(
        parsed.select(
            'chain {} and resnum {} and name {}'.format(
                ifginfo[0],
                ifginfo[1],
                atom)).getCoords()[0])

for atom in vdm1atoms:
    query.append(
        parsed.select(
            'chain {} and resnum {} and name {}'.format(
                vdm1info[0],
                vdm1info[1],
                atom)).getCoords()[0])

for atom in vdm2atoms:
    query.append(
        parsed.select(
            'chain {} and resnum {} and name {}'.format(
                vdm2info[0],
                vdm2info[1],
                atom)).getCoords()[0])

query = np.array(query)

rmsds = []
num_atoms = len(query)
print('num atoms', num_atoms)
ifglists = flip(ifgatoms, ifgresn)

ifgcsv = pd.read_csv(tyr_ifg_csv, index_col='iFG_count')
tyr = pd.read_csv(tyr_csv)
tyr = tyr.groupby('iFG_count')
for ifgcount, group in tyr:
    if 'TRP' in list(group['resname']) and 'PRO' in list(group['resname']):
        pdb = ifgcsv.iloc[ifgcount]['pdb']
        db_dir = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'
        par = pr.parsePDB(db_dir + pdb + 'H.pdb')
#        ifgchid,ifgresnum = row['chid_ifg'],row['resnum_ifg']
#        vdmchid,vdmresnum = row['chid_vdm'],row['resnum_vdm']
#        ifgseg,vdmseg = row['segi_ifg'],row['segi_vdm']
#        ifg_all_atoms = constants.atoms_dict[ifgname]
#        vdm_all_atoms = constants.atoms_dict[vdmname]
#        for atom in ifg_all_atoms:
#            sele = 'chain {} and segment {} and resnum {} and name {}'.format(ifgchid,ifgseg,ifgresnum,atom)
#            ifg.append(par.select(sele).getCoords()[0])
#        for atom in vdm_all_atoms:
#            sele = 'chain {} and segment {} and resnum {} and name {}'.format(vdmchid,vdmseg,vdmresnum,atom)
#            vdm.append(par.select(sele).getCoords()[0])
#        coords.append([row.name,np.array(ifg),np.array(vdm)])
#                else:
#                    coords.append([row.name,None])
#            except:
#                traceback.print_exc()
#                coords.append([row.name,None])
#    pkl.dump(coords, open(look_dir+'coords_of_{}_multiple_{}.pkl'.format(aa, multiple_of_20),'wb'))
#
#
#
#
#

        #sele = 
        #print(group.columns)

        # grab coordinates 


#lookupatoms_to_clus = []
#lookupatoms_to_clus.append(query) # first element to see how big its cluster is
#counter = 0 # to keep count of how many pdbs are being output
#for item in coords_ls:
#    if len(item) == 3:
#        compare_rmsds = []
#        ifg_vdm_ind = []
#        for ifg_ind, ifgls in enumerate(ifglists):
#            for vdm_ind, vdmls in enumerate(vdmlists):
#                lookupatoms = get_order_of_atoms(
#                    item, ifgresn, vdmresn, ifgls, vdmls)
#                moved, transf = pr.superpose(lookupatoms, query)
#                temp_rmsd = pr.calcRMSD(moved, query)
#                compare_rmsds.append(temp_rmsd)
#                ifg_vdm_ind.append([moved, temp_rmsd])
#        # item[0] is df index
#        rmsds.append([item[0], min(compare_rmsds)])
#        # get index of which one had min rmsd
#        for which_ind, each in enumerate(ifg_vdm_ind):
#            if each[1] == min(compare_rmsds):
#                lookupatoms_to_clus.append(each[0])
#
#                ############################################################ 
#                #                output pdb if low rmsd 
#                ############################################################ 
#                if each[1] < cutoff and counter < 30 and which_ind==0: 
#                    # this is to ensure rmsd is below cutoff when not flipped
#                    # bc don't want to take care of that in prody to output pdb
#                    row = lookupdf.loc[item[0]]
#                    try:
#                        db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'
#                        par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
#                    except:
#                        db_dir = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'
#                        par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
#    
#                    ifgchid, ifgresnum = row['chid_ifg'], row['resnum_ifg']
#                    vdmchid, vdmresnum = row['chid_vdm'], row['resnum_vdm']
#                    printout = copy.deepcopy(par)
#                    printout = printout.select(
#                        '(chain {} and resnum {}) or (chain {} and resnum {})'.format(
#                        ifgchid,ifgresnum,vdmchid,vdmresnum))
#                    printout.select('chain {} and resnum {}'.format(ifgchid,ifgresnum)).setChids('Y')
#                    printout.select('chain {} and resnum {}'.format(vdmchid,vdmresnum)).setChids('X')
#                    printout.select('all').setResnums(10)
#                    printout_interactamer = []
#                    integrin_interactamer = []
#                    try: # skip the ones that have segment ids.
#                        for atom in ifgatoms:
#                            integrin_interactamer.append(parsed.select('chain {} and resnum {} and name {}'.format(ifginfo[0],ifginfo[1],atom)))
#                            printout_interactamer.append(printout.select('chain Y and resnum 10 and name {}'.format(atom)))
#                        for atom in vdmatoms:
#                            integrin_interactamer.append(parsed.select('chain {} and resnum {} and name {}'.format(vdminfo[0],vdminfo[1],atom)))
#                            printout_interactamer.append(printout.select('chain X and resnum 10 and name {}'.format(atom)))
#                        integrin_interactamer_prody = []
#
#                        
#                        integrin_interactamer = sum(integrin_interactamer[1:], integrin_interactamer[0])
#                        printout_interactamer = sum(printout_interactamer[1:], printout_interactamer[0])
#                        try:
#                            assert len(integrin_interactamer) == len(printout_interactamer)
#
#                            interact_res = printout.select('(chain X and resnum 10) or (chain Y and resnum 10)')
#                            interactamer_transf = pr.applyTransformation(transf, printout_interactamer)
#                            
#                        except:
#                            pass
#                    except:
#                        traceback.print_exc()
#                        pass
#
#    else:
#        rmsds.append([int(item[0]), 100000])
#
## count how many coords are within rmsd cutoff
#rmsds = [i for i in rmsds if i[1] < cutoff]
#print(len(rmsds), 'number of matches')
