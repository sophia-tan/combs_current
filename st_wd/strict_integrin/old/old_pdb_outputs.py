import sys, os
#sys.path.append('/home/gpu/Sophia/combs/src/')
#from combs.apps import *
from residues_integrin import *
import prody as pr
import numpy as np, pickle as pkl, pandas as pd
from Scoring import *

script,inputres = sys.argv

rmsd_outputs = './output_data/'
lookupvdms = '/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/'
integrin = pr.parsePDB('pymolintegrin.pdb')

for resi, resn in integrin_res.items():
    #if 1==1:
    if resi[1:] == inputres:
        for pklf in sorted(os.listdir(rmsd_outputs)):
            if pklf.startswith(resi) and pklf.endswith('intrxn.pkl'):
                pkl_list = pkl.load(open(rmsd_outputs+pklf,'rb'))
                # first row is [# atoms, ifgatoms, vdmatoms]
                # all the other rows are [index in DF, rmsd]
                indices = pkl_list[1:]
                indices = indices[:100] # deletel
                indices = [x[0] for x in indices if x[1]<.5]
                indices = indices[:50]
                targ_res = pklf.split('_')[1][3:]
                int_res = pklf.split('_')[2]
                targ_res = constants.AAname[targ_res]
                vdms_df_path = lookupvdms+'vdms_of_{}.pkl'.format(targ_res)
                pdb_df = pkl.load(open(vdms_df_path,'rb'))
                for ix in indices:
                    row = pdb_df.loc[ix]
                    try:
                        db_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/'
                        par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
                    except:
                        db_dir = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'
                        par = pr.parsePDB(db_dir+row['pdb']+'H.pdb')
                    ifgchid, ifgresnum = row['chid_ifg'], row['resnum_ifg']
                    vdmchid, vdmresnum = row['chid_vdm'], row['resnum_vdm']
                    printout = copy.deepcopy(par)
                    printout = printout.select('(chain {} and resnum {}) or (chain {} and resnum {})'.format(ifgchid,ifgresnum,vdmchid,vdmresnum))
                    printout.select('chain {} and resnum {}'.format(ifgchid,ifgresnum)).setChids('Y')
                    printout.select('chain {} and resnum {}'.format(vdmchid,vdmresnum)).setChids('X')
                    printout.select('all').setResnums(10)
                    printout_interactamer = []
                    printout_interactamer.append(printout.select('chain Y and resnum 10 and not element H D'))
                    printout_interactamer.append(printout.select('chain X and resnum 10 and not element H D'))
                    printout_interactamer = printout_interactamer[0]+printout_interactamer[1]
                    integrin_interactamer = []
                    integrin_interactamer.append(integrin.select('chain A and resnum {}'.format(resi[1:])))
                    try:
                        vdm = integrin.select('chain B and resnum {}'.format(int_res[:3]))
                        assert vdm.getResnames()[0]==int_res[3:]
                        integrin_interactamer.append(vdm)
                    except:
                        vdm = integrin.select('chain A and resnum {}'.format(int_res[:3]))
                        assert vdm.getResnames()[0]==int_res[3:]
                        integrin_interactamer.append(vdm)
                    integrin_interactamer = integrin_interactamer[0]+integrin_interactamer[1]
                    integrin_interactamer = integrin_interactamer.select('not name OXT')
                    printout_interactamer = printout_interactamer.select('not name OXT')
                    for index in range(len(integrin_interactamer)):
                        if list(integrin_interactamer)[index].getResname() != list(printout_interactamer)[index].getResname():
                            print(integrin_interactamer.getNames())
                            print(printout_interactamer.getNames())
                        
                    moved,transf = pr.superpose(printout_interactamer, integrin_interactamer)
                    print(ix,int_res,pr.calcRMSD(moved,integrin_interactamer))
                    print(pr.calcRMSD(moved,integrin_interactamer))
                    
                    outdir = '/home/gpu/Sophia/combs/st_wd/strict_integrin/output_data/pdbfiles/'
                    pr.writePDB(outdir+'{}_{}_{}.pdb'.format(resi,int_res,row.name),moved)
