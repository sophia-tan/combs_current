import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
import pickle
import collections
# import time

# class MacOSFile(object):
#
#     def __init__(self, f):
#         self.f = f
#
#     def __getattr__(self, item):
#         return getattr(self.f, item)
#
#     def read(self, n):
#         # print("reading total_bytes=%s" % n, flush=True)
#         if n >= (1 << 31):
#             buffer = bytearray(n)
#             idx = 0
#             while idx < n:
#                 batch_size = min(n - idx, 1 << 31 - 1)
#                 # print("reading bytes [%s,%s)..." % (idx, idx + batch_size), end="", flush=True)
#                 buffer[idx:idx + batch_size] = self.f.read(batch_size)
#                 # print("done.", flush=True)
#                 idx += batch_size
#             return buffer
#         return self.f.read(n)
#
#     def write(self, buffer):
#         n = len(buffer)
#         print("writing total_bytes=%s..." % n, flush=True)
#         idx = 0
#         while idx < n:
#             batch_size = min(n - idx, 1 << 31 - 1)
#             print("writing bytes [%s, %s)... " % (idx, idx + batch_size), end="", flush=True)
#             self.f.write(buffer[idx:idx + batch_size])
#             print("done.", flush=True)
#             idx += batch_size
#
#
# def pickle_dump(obj, file_path):
#     with open(file_path, "wb") as f:
#         return pickle.dump(obj, MacOSFile(f), protocol=pickle.HIGHEST_PROTOCOL)
#
#
# def pickle_load(file_path):
#     with open(file_path, "rb") as f:
#         return pickle.load(MacOSFile(f))

def pickle_load(file_path):
    with open(file_path, "rb") as f:
        return pickle.load(f)

pdb_path = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/1wdn_noligandH.pdb'
sample = combs.Sample()
sample.poi = pr.parsePDB(pdb_path)
sample.bs_residues = list(zip([70, 75, 119, 118, 50, 67, 10, 115, 156, 13, 68, 157, 185],
                              ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A']))
# sample.set_rois()
# sample.set_rois_phipsi()
# sample.set_poi_clash_coords()
# sample.set_roi_bb_coords()
sample.lig_ifg_corr = collections.defaultdict(dict)
sample.lig_ifg_coords = collections.defaultdict(dict)
sample.hotspot_cliques = {}

indir = '/Users/npolizzi/Projects/combs/src/runs/glutamine_binding_protein/20170911/'
relvdm_carboxamide = pickle_load(indir + 'relvdm_carboxamide.pickle')
relvdm_carboxylate = pickle_load(indir + 'relvdm_carboxylate.pickle')
relvdm_amino = pickle_load(indir + 'relvdm_amino.pickle')

relvdm_carboxamide.hotspot_ifg_coords = collections.defaultdict(dict)
relvdm_carboxylate.hotspot_ifg_coords = collections.defaultdict(dict)
relvdm_amino.hotspot_ifg_coords = collections.defaultdict(dict)

relvdm_carboxamide.set_hotspot_ifg_coords()
relvdm_carboxamide.ifg_dict = {'GLN': 'CG CD OE1 NE2'}
relvdm_carboxylate.set_hotspot_ifg_coords()
relvdm_carboxylate.ifg_dict = {'GLU': 'CG CD OE1 OE2'}
relvdm_amino.set_hotspot_ifg_coords()
relvdm_amino.ifg_dict = {'LYS': 'CD CE NZ'}

sample.ligand_conformers = [pr.parsePDB('/Users/npolizzi/Desktop/gln227_mm_good.pdb')]
sample.set_ligand_ifg_correspondence(relvdm_carboxamide, 'NE2', 'NE2')
sample.set_ligand_ifg_correspondence(relvdm_carboxylate, 'O', 'OE2')
sample.set_ligand_ifg_correspondence(relvdm_amino, 'N', 'NZ')
sample.find_ligand_cst_hotspots([relvdm_carboxamide, relvdm_amino, relvdm_carboxylate], tol=0.3)


