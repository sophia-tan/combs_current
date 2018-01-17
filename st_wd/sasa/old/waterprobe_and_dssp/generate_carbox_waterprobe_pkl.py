import pandas as pd
import prody as pr
import numpy as np
import freesasa as fs
import pickle as pkl

# vdm info
df = pkl.load(open('all_vdms_df.pkl','rb'))

uniquepdbs = set(df['pdb'])
indices = df.index.values
newdf = pd.DataFrame(columns=['iFG_count', 'vdM_count', '1.4A'], index=indices)
radius = 1.4
for pdb in uniquepdbs:
    print(pdb)
    try:
        vdms = df[df['pdb'] == pdb ]
        struct = fs.Structure('../reduce/'+pdb+'H.pdb')
        result = fs.calc(struct)
        parsed = pr.parsePDB(pdb)
        cb_sele = parsed.select('protein and (backbone or name CB) and not element H D')
        coords = list(x for y in cb_sele.getCoords() for x in y)
        radii = list(fs.Classifier().radius(x, y) for x, y in zip(cb_sele.getResnames(), cb_sele.getNames()))
        sasa = fs.calcCoord(coords, radii, fs.Parameters({'probe-radius': radius}))

        for ix, row in vdms.iterrows():
            try:
                ifgcount = row['iFG_count']
                vdmcount = row['vdM_count']
                vdmresnum = row['resnum_vdm']
                vdmchid = row['chid_vdm']
                index = row.name

                sele = cb_sele.select('resnum %s' %vdmresnum)
                bb_cb_atom_ind = sele.getIndices()
                prody_pdb_bb_cb_atom_ind = cb_sele.getIndices()

                largeprobesasa = '{0:.2f}'.format(sum(sasa.atomArea(i) for i in np.where(np.in1d(prody_pdb_bb_cb_atom_ind,bb_cb_atom_ind))[0]))

                newdf.ix[index, 'iFG_count'] = ifgcount
                newdf.ix[index, 'vdM_count'] = vdmcount
                newdf.ix[index, '1.4A'] = largeprobesasa
            except:
                pass
    except:
        pass
pkl.dump(newdf, open('tryexceptcarboxamide_vdms_waterprobe_1.4_sasa.pkl','wb'))
