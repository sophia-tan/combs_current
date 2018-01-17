# run_comb.py equivalent of database non-combing statistics
from pprint import pprint
__all__ = ['database_sasa']

from .Comb import Comb
from .nrPDB import nrPDB
from .ParsedUncombedPDB import ParsedUncombedPDB
from .IntFG import IntFG
from .Vandermer import Vandermer
from .Res import Res
import sys
import os
import time
import traceback
import prody as pr
import pandas as pd

def database_sasa(**kwargs):
    # lookup table! key=aa, value=dict(where key=sasatype, value=list of sasa scores)
    sasa_dict = {}
    AAs = ['ALA', 'ARG', 'VAL', 'MET', 'GLU', 'GLN', 'GLY', \
    'ASP', 'ASN', 'ILE', 'LEU', 'LYS', 'HIS', 'TRP', 'TYR', 'PHE', 'PRO', 'THR', 'SER', 'CYS']
    for aa in AAs:
        sasa_dict[aa] = {'3A': [], '4A': [], '5A': [], 'dssp':[], 'fs':[]}
    
    comb_start = time.time()
    with nrPDB(**kwargs) as cb:
        total_pdbs = sum(len(chains) for chains in cb.pdbs_chains.values())
        for i, pdb_acc in enumerate(cb.pdbs_chains.keys()):
            print('@> Combing PDB ' + pdb_acc + ', ' + str(i + 1) + ' of ' + str(total_pdbs))
            start_time = time.time()
            for chain in cb.pdbs_chains[pdb_acc]:
                try:
                    pdb = ParsedUncombedPDB(cb, pdb_acc, chain)
                    pr_pdb = pdb.prody_pdb
                    # make list of all residues in pr_pdb to iterate
                    sel = pr_pdb.select('chain %s and calpha' % chain)
                    resnum_list = sel.getResnums()
                    for resnum in resnum_list: # make a Res instance (similar to Vandermer)
                        if resnum > 0:
                            res_sel = pr_pdb.select('chain %s resnum %s' %(chain, str(resnum)))
                        else: # negative resnum!
                            res_sel = pr_pdb.select('chain %s resnum `%s`' %(chain, str(resnum)))
                        res = Res(res_sel)
                        res.get_info(pdb) # this calculates all the sasa types!
                        sasa_types = [res.sasa_3A_probe, res.sasa_4A_probe, res.sasa_5A_probe, res.dssp_sasa, res.residue_sasa]
                        
                        ###############
                        # add these sasa scores to the lookup table! 
                        aa_dict = sasa_dict[res.resname]
                        for ix, key in enumerate(aa_dict): # key is the sasa method
                            sasa_dict[res.resname][key].append(sasa_types[ix])
                            
                        
                except Exception:
                    traceback.print_exc()
    return sasa_dict
