import pandas as pd, pickle as pkl

integrin_res = {}
#
integrin_res['A671'] = 0.54
integrin_res['A673'] = 0.23
integrin_res['A753'] = 0.37
integrin_res['A755'] = 0.4
integrin_res['A758'] = 0.64
integrin_res['A760'] = 0.73
integrin_res['A785'] = 0.47
integrin_res['A787'] = 0.17
integrin_res['A900'] = 0.31

integrin_res['B552'] = 0.18
integrin_res['B594'] = 0.46
integrin_res['B603'] = 0.82
integrin_res['B626'] = 0.16
integrin_res['B658'] = 0.44
integrin_res['B664'] = 0.54

list_for_df = []
for key, val in integrin_res.items():
    list_for_df.append(
        pd.Series([key, val], index= ['Mutation', 'expt val']))
df = pd.DataFrame(list_for_df)
pkl.dump(df, open('experimental_data_AiibB3.pkl','wb'))
