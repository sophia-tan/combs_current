import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import pandas as pd, prody as pr
import pickle as pkl
from sys import argv

script, ifg_res = argv

csv_dir = '/home/gpu/Sophia/combs/st_wd/20180207db_combed_csvs/%s/'%ifg_res
db_pdbs = '/home/gpu/Sophia/combs/st_wd/20180207_db_molprobity_biolassem/'

an = analysis.Analyze(csv_dir)
df = an.get_distant_vdms(seq_distance=10)
print(len(df), 'dist')

df = analysis.Analysis.remove_repeat_proteins(df)
failed = 0

ls=[]
#def parse_vdm(row):
for ix, row in df.iterrows():
    '''get prody selection of just ifg (rename to chain Y) and 
    vdm (chain X). input DF, returns Series of prody selections'''
    parsed = pr.parsePDB('1wdn')
    #chid_ifg,resname_ifg = row['chid_ifg'],row['resname_ifg']
    #resnum_ifg,seg_ifg = row['resnum_ifg'], row['segi_ifg']
    #chid_vdm,resname_vdm = row['chid_vdm'],row['resname_vdm']
    #resnum_vdm,seg_vdm = row['resnum_vdm'], row['segi_vdm']
    #pdb = row['pdb']
    #ifgsele = 'chain %s and resnum %s and segment %s'%(chid_ifg, str(resnum_ifg), seg_ifg)
    #vdmsele = 'chain %s and resnum %s and segment %s'%(chid_vdm, str(resnum_vdm), seg_vdm)
    #parsed = pr.parsePDB(db_pdbs+pdb+'H.pdb').select('({}) or ({})'.format(ifgsele,vdmsele))
    #parsed.select(ifgsele).setChids('Y')
    #parsed.select(vdmsele).setChids('X')
    ls.append(parsed)
    #return parsed
    #try:
    #    parsed.select(ifgsele).setChids('Y')
    #    parsed.select(vdmsele).setChids('X')
    #    return parsed
    #except:
    #    print(pdb, 'FAILED')
    #    print(ifgsele)
    #    print(vdmsele)
    #    #global failed 
    #    #failed += 1
    #    return None
#parsedvdms=pd.Series(df.apply(parse_vdm,axis=1),name='prody_sel')
#print('total failed', failed)
#parsedvdms=pd.concat([df,parsedvdms],axis=1)[['resname_ifg','resname_vdm','prody_sel']]
#pkl.dump(parsedvdms, open('/home/gpu/Sophia/combs/st_wd/Lookups/refinedvdms/\
#    parsedvdms_of_{}.pkl'.format(ifg_res),'wb'))
