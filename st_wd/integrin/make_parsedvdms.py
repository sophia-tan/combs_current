import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
import pandas as pd
import pickle as pkl
from sys import argv

script, ifg, vdm = argv

vdmdict = {}
vdmdict['phenyl'] = ['PHE']
vdmdict['guanidino'] = ['ARG']
vdmdict['tyrCOH'] = ['TYR']
vdmdict['carboxylate'] = ['GLU', 'ASP']

csv_dir = '/home/gpu/Sophia/STcombs/20171118/%s/csv/'%ifg
an = analysis.Analyze(csv_dir)
df = an.get_distant_vdms(seq_distance=10)
print(len(df), 'dist')
def parsedist(row):
    '''helper fx for backboneCO to only include vdms with bbCO within 4A'''
    row = row[1:-1]
    row = row.split(') (')
    row = [x.split(' ') for x in row]
    row = [1 for x in row if x[1] is 'O' and float(x[2]) < 4]
    return sum(row) > 0

if vdm == 'backboneCO':
    print('is bbCO')
    contact = an.ifg_contact_vdm
    # find out which vdMs have bb O that's < 4A away
    print(len(contact), 'len contact')
    contact['dist_info'] = contact['dist_info'].apply(parsedist)
    print(contact['dist_info'].sum(), 'expected len')
    chopped = contact[contact['dist_info']==True]
    print(len(chopped),'len chopped')
    # merge back to main df
    df = pd.merge(df, chopped, how='inner',on=['iFG_count', 'vdM_count'],left_index=True,right_index=True)
    print(len(df),'after merged')
df = analysis.Analysis.remove_repeat_proteins(df)

try:
    dict_parsed = an.parse_vdms_by_aa(df=df, subset=vdmdict[vdm])
    ls_parsed = []
    for key in dict_parsed.keys():
        ls_parsed = ls_parsed + dict_parsed[key]
    pkl.dump(ls_parsed, open('parsedvdms_%s_%s.pkl'%(ifg,vdm), 'wb'))
except:
    parsed_vdms = an.parse_vdms(df)
    pkl.dump(parsed_vdms, open('parsedvdms_%s_%s.pkl'%(ifg,vdm), 'wb'))
#
### add sasa info
##sasadf = an.vdm_sasa_info
##merged = pd.merge(df,sasadf, on=['iFG_count', 'vdM_count'])
##
##pkl.dump(merged, open('noskip_%s_sasa.pkl' %ifg, 'wb'))
