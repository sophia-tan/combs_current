import pandas as pd
import prody as pr

tyr_csv = '/home/gpu/Sophia/combs/st_wd/20180207_combed_csvs/tyrosine/tyrosine_vdm_pdb_info.csv'

tyr = pd.read_csv(tyr_csv)
tyr = tyr.groupby('iFG_count')
for ifgcount, group in tyr:
    print(group)
