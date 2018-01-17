import pickle as pkl

df = pkl.load(open('ss_localfrag_database_lookups.pkl','rb'))
hbonds = pkl.load(open('dist_hbonds_ss.pkl','rb'))

for ix, row in df.iterrows():
    total = sum(row[:7])
    for col in df.columns.values[:7]:
        df.ix[ix, col] = round(row[col] / total,2)

for ix, row in hbonds.iterrows():
    total = sum(row[:7])
    for col in hbonds.columns.values[:7]:
        hbonds.ix[ix, col] = round(row[col] / total, 2)
    if ix == 5:
        x = (row[:5] / df.ix[4][:5])
        x = [round(z,2) for z in x]
        hbonds.ix[ix][:5] = x
    else:
        x = (row[:5] / df.ix[ix][:5])
        x = [round(z,2) for z in x]
        hbonds.ix[ix][:5] = x

print(hbonds)
