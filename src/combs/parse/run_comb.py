__all__ = ['run_comb']

from .comb import Comb
from .parsedPDB import ParsedPDB
from .intFG import IntFG
from .vandermer import Vandermer
import sys
import os
import time
import traceback
# from multiprocessing import Pool


def run_comb(ifg_dict, **kwargs):
    comb_start = time.time()
    with Comb(ifg_dict, **kwargs) as cb:
        ifg_info = cb.ifg_sele_dict.values()
        new_dicts = []
        radius1 = kwargs.get('radius1', 3.5)
        radius2 = kwargs.get('radius2', 4.8)
        for dict_ in ifg_info:
            ifg_parts = []
            for resn, atoms in dict_.items():
                ifg_parts.append(resn + ': ' + atoms)
                if resn == 'ANY':
                    break
            new_dicts.append('(' + ', '.join(ifg_parts) + ')')
        print('@> Combing for ' + '; '.join(new_dicts))
        print()
        total_pdbs = sum(len(chains) for chains in cb.pdbs_chains.values())
        # with Pool(4) as p:
        #     p.starmap(comb_inner_loop, zip([cb]*total_pdbs, [total_pdbs]*total_pdbs, list(range(0,total_pdbs)),
        #               [pdb_acc for pdb_acc in cb.pdbs_chains.keys()]))
        for i, pdb_acc in enumerate(cb.pdbs_chains.keys()):
            print('@> Combing PDB ' + pdb_acc + ', ' + str(i + 1) + ' of ' + str(total_pdbs))
            start_time = time.time()
            for chain in cb.pdbs_chains[pdb_acc]:
                try:
                    pdb = ParsedPDB(cb, pdb_acc, chain)
                    while pdb.possible_ifgs:
                        ifg = IntFG(pdb, cb)
                        ifg.find_contact_atoms(pdb, cb, radius1=radius1, radius2=radius2)
                        if ifg.contact_atoms_protein is not None:
                            ifg.get_hbonds(pdb)
                            ifg.get_ca_hbonds(pdb)
                            while ifg.contact_resnums:
                                vdm = Vandermer(ifg, pdb)
                                vdm.print_pdb(ifg, pdb, cb)
                                vdm.send_info(ifg, pdb)
                            ifg.send_info(pdb)
                    pdb.write_csv(cb)
                    total_time = time.time() - start_time
                    print('@> PDB ' + pdb_acc + ' combed in ' + '{0:.2f}'.format(total_time) + ' sec')
                except Exception:
                    traceback.print_exc()
        cb.write_poss_ifg_csv('total_poss_ifgs', ['total_poss_iFGs'])
        total_comb_time = time.time() - comb_start
        print('@> Database combed in ' + '{0:.2f}'.format(total_comb_time) + ' sec')


# def comb_inner_loop(cb, total_pdbs, i, pdb_acc):
#     print('@> Combing PDB ' + pdb_acc + ', ' + str(i + 1) + ' of ' + str(total_pdbs))
#     start_time = time.time()
#     try:
#         pdb = ParsedPDB(cb, pdb_acc)
#     except Exception:
#         #     print('Parsing Error for PDB ' + pdb_acc)
#         #     # print(type(inst))
#         #     # print(inst.args)
#         #     # print(inst)
#         exc_type, exc_obj, exc_tb = sys.exc_info()
#         fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
#         print(exc_type, fname, exc_tb.tb_lineno, exc_obj)
#     while pdb.possible_ifgs:
#         # try:
#         ifg = IntFG(pdb, cb)
#         ifg.find_contact_atoms(pdb, cb)
#         if ifg.contact_atoms_all is not None:
#             ifg.get_hbonds(pdb)
#             ifg.get_ca_hbonds(pdb)
#             while ifg.contact_resnums:
#                 # try:
#                 vdm = Vandermer(ifg, pdb)
#                 # print(vdm.chid, vdm.resnums)
#                 vdm.print_pdb(ifg, pdb, cb)
#                 vdm.send_info(ifg, pdb)
#                 # except Exception:
#                 #     print('vdM instantiation Error in PDB ' + pdb_acc + ' for iFG ' + str(ifg.count))
#                 #     # print(type(inst))
#                 #     # print(inst.args)
#                 #     # print(inst)
#                 #     exc_type, exc_obj, exc_tb = sys.exc_info()
#                 #     fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
#                 #     print(exc_type, fname, exc_tb.tb_lineno, exc_obj)
#             ifg.send_info(pdb)
#             # except Exception:
#             #     print('iFG instantiation Error in PDB ' + pdb_acc)
#             #     # print(type(inst))
#             #     # print(inst.args)
#             #     # print(inst)
#             #     # exc_type, exc_obj, exc_tb = sys.exc_info()
#             #     exc_type = sys.exc_info()[0]
#             #     exc_tb = sys.exc_info()[-1]
#             #     fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
#             #     # print(exc_type, fname, exc_tb.tb_lineno, exc_obj)
#             #     print(exc_type, fname, exc_tb.tb_lineno)
#     pdb.write_csv(cb)
#     total_time = time.time() - start_time
#     print('@> PDB ' + pdb_acc + ' combed in ' + '{0:.2f}'.format(total_time) + ' sec')