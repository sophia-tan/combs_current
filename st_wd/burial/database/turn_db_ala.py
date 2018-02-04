# generate lookup table for burial (min dist to convex hull) in nrPDB database
import sys, os, subprocess
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import numpy as np, pandas as pd, pickle as pkl
from numpy import concatenate, array_equal, dot, cross, linalg, abs
from scipy.spatial import ConvexHull, Delaunay
from prody import *
import traceback

# directory of PDBs in database that have been put through reduce
pdb_dir = '/home/gpu/Sophia/STcombs/20171118/database/reduce/' 

# compile lists of pdbs and chains
pdbnames = []
chains = []

with open('combs_vast_50seqid_blastclust_07252017.txt') as fo:
    for line in fo:
        pdbnames.append(line[:4])
        chains.append(line[4])
zipped = zip(pdbnames, chains)

for pdb, chain in zipped:
    try:
        pdbfile = pdb_dir+pdb+'H.pdb'
        orig_pdb = parsePDB(pdbfile, altloc='A', model=1) # to preserve resnum and resname
        orig_pdb = orig_pdb.select('chain %s and name CA'%chain) # get CA instead of CB bc some glycines
        
        # use Rosetta's fixbb to mutate to all-ala to get C-beta of every residue
        try:
            parsed = parsePDB(pdb_dir+pdb+'H_0001.pdb', altloc='A', model=1) 
        except:
            try:
                subprocess.call(['/home/gpu/src/rosetta3.7/rosetta_src_2016.32.58837_bundle/main/source/bin/fixbb.default.linuxgccrelease', '-overwrite', '-mute', 'basic', '-mute', 'protocols', '-mute', 'core', '-ignore_unrecognized_res', '-resfile', 'alanine.res', '-out:prefix', '%s'%pdb_dir, '-s', '%s'%pdbfile])
                parsed = parsePDB(pdb_dir+pdb+'H_0001.pdb', altloc='A', model=1)
            except: # if rosetta gives error during fixbb, output pdb from prody and do fixbb again
                pr_pdb = parsePDB(pdbfile).select('chain %s'%chain)
                writePDB('temp'+pdb,pr_pdb)
                subprocess.call(['/home/gpu/src/rosetta3.7/rosetta_src_2016.32.58837_bundle/main/source/bin/fixbb.default.linuxgccrelease', '-overwrite', '-mute', 'basic', '-mute', 'protocols', '-mute', 'core', '-ignore_unrecognized_res', '-resfile', 'alanine.res', '-out:prefix', '%s'%pdb_dir, '-s', 'temp%s.pdb'%pdb])
                os.unlink('temp%s.pdb'%pdb)# delete the tempfile  
                os.rename(pdb_dir+'temp%s_0001.pdb'%pdb,pdb_dir+'%sH_0001.pdb'%pdb)
                parsed = parsePDB(pdb_dir+pdb+'H_0001.pdb', altloc='A', model=1)
        
    except Exception:
        print(traceback)
        print(pdb, 'pdb could not be parsed')
