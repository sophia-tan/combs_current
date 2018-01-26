__all__ = ['make_bb_sc_rel_vdms', 'within_threepfive', 'has_bb_or_sc', 'make_rel_vdm_coords']

import prody as pr
import pandas as pd
import numpy as np
import traceback
import os
import pickle
from ..apps.constants import interactamer_atoms, flip_sets, residue_sc_names, flip_residues
from ..analysis import Analyze
from ..analysis import Analysis
from ..analysis.cluster import cluster_adj_mat
from sklearn.neighbors import NearestNeighbors
import shutil

def listdir(path):
    return [file for file in os.listdir(path) if file[0] != '.']

def get_centroids(picklepath, radius=0.2):
    """This will cluster all vdms by iFG location (backbone rel_vdms) or by sidechain + iFG location (sc rel_vdms)
    and output new pickle files to a directory picklepath/clustered."""
    if picklepath[-1] != '/':
        picklepath += '/'
    if os.path.isdir(picklepath):
        for pickletype in listdir(picklepath):
            if pickletype == 'PHI_PSI':
                for phipsi_type in listdir(picklepath + pickletype):
                    for picklefile in listdir(picklepath + pickletype + '/' + phipsi_type + '/pickle/'):
                        with open(picklepath + pickletype + '/' + phipsi_type + '/pickle/'
                                  + picklefile, 'rb') as infile:
                            pick = pickle.load(infile)
                            if len(pick.shape) == 1:
                                with open(outpath + picklefile, 'wb') as outfile:
                                    pickle.dump(pick, outfile)
                            else:
                                ifg_flat = [coords.flatten() for coords in pick[:, -2]]
                                nbrs = NearestNeighbors(metric='euclidean', radius=radius)
                                nbrs.fit(ifg_flat)
                                adj_mat = nbrs.radius_neighbors_graph(ifg_flat)
                                mems, cents = cluster_adj_mat(adj_mat)
                                outpath = picklepath + 'clustered/' + pickletype + '/' + phipsi_type + '/pickle/'
                                try:
                                    os.makedirs(outpath)
                                except:
                                    pass
                                with open(outpath + picklefile, 'wb') as outfile:
                                    pickle.dump(pick[cents, :], outfile)
            else:
                for picklefile in listdir(picklepath + pickletype + '/pickle/'):
                    with open(picklepath + pickletype + '/pickle/' + picklefile, 'rb') as infile:
                        pick = pickle.load(infile)
                        outpath = picklepath + 'clustered/' + pickletype + '/pickle/'
                        try:
                            os.makedirs(outpath)
                        except:
                            pass
                        if len(pick.shape) == 1:
                            with open(outpath + picklefile, 'wb') as outfile:
                                pickle.dump(pick, outfile)
                        else:
                            if pickletype == 'SC':
                                sc_flat = [coords.flatten() for coords in pick[:, -3]]
                                ifg_flat = [coords.flatten() for coords in pick[:, -2]]
                                sc_ifg_flat = np.hstack((sc_flat, ifg_flat))
                                nbrs = NearestNeighbors(metric='euclidean', radius=radius)
                                nbrs.fit(sc_ifg_flat)
                                adj_mat = nbrs.radius_neighbors_graph(sc_ifg_flat)
                                mems, cents = cluster_adj_mat(adj_mat)
                                with open(outpath + picklefile, 'wb') as outfile:
                                    pickle.dump(pick[cents, :], outfile)
                            else:
                                ifg_flat = [coords.flatten() for coords in pick[:, -2]]
                                nbrs = NearestNeighbors(metric='euclidean', radius=radius)
                                nbrs.fit(ifg_flat)
                                adj_mat = nbrs.radius_neighbors_graph(ifg_flat)
                                mems, cents = cluster_adj_mat(adj_mat)
                                with open(outpath + picklefile, 'wb') as outfile:
                                    pickle.dump(pick[cents, :], outfile)


def make_all_rel_vdms(path_to_csv, outpath, comb, seq_distance=5):
    if outpath[-1] != '/':
        outpath += '/'
    an = Analyze(path_to_csv)
    dist_vdms = an.get_distant_vdms(seq_distance=seq_distance)
    dist_vdms_sc = get_vdms_sc(dist_vdms)
    make_rel_vdms(dist_vdms_sc, an, outpath + 'SC', comb, 'CA', 'N', 'C')
    dist_vdms_bb_N_CA = get_vdms_bb_N_CA(dist_vdms)
    make_rel_vdms(dist_vdms_bb_N_CA, an, outpath + 'N_CA', comb, 'N', 'H', 'CA')
    dist_vdms_bb_C_O = get_vdms_bb_C_O(dist_vdms)
    make_rel_vdms(dist_vdms_bb_C_O, an, outpath + 'C_O', comb, 'C', 'O', 'CA')
    dist_vdms_bb_phi_psi = get_vdms_bb_phi_psi(dist_vdms)
    for phi_low, phi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
        for psi_low, psi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
            outpath_phipsi = outpath + 'PHI_PSI/' + str(phi_low) + '_' + str(phi_high) + '_' + str(psi_low) \
                     + '_' + str(psi_high)
            def phi_psi_select(row):
                if type(row['sec_struct_phi_psi_vdm']) is float:
                    return False
                if 'None' in row['sec_struct_phi_psi_vdm']:
                    return False
                spltrow = row['sec_struct_phi_psi_vdm'].split()
                if (phi_low <= float(spltrow[0]) < phi_high) and (psi_low <= float(spltrow[1]) < psi_high):
                    return True
                else:
                    return False
            dist_vdms_phipsi = dist_vdms_bb_phi_psi[dist_vdms_bb_phi_psi.apply(phi_psi_select, axis=1)]
            make_rel_vdms(dist_vdms_phipsi, an, outpath_phipsi, comb, 'CA', 'N', 'C')


def get_sc_hb(row):
    try:
        index = [i for i, rrn in enumerate([n for n in row['rel_resnums'] if n != '-']) if rrn == '0']
        return any([any(
            {y for names in [row['vdM_atom_names'].strip('()').split(') (')[ind]] for y in names.split() if
             y not in ['N', 'O', 'H', 'OXT', 'H1', 'H2', 'H3']}) for ind in index])
    except:
        return False


def get_bb_hb(row):
    try:
        index = [i for i, rrn in enumerate([n for n in row['rel_resnums'] if n != '-']) if rrn == '0']
        return any([any(
            {y for names in [row['vdM_atom_names'].strip('()').split(') (')[ind]] for y in names.split() if
             y in ['N', 'O', 'H']}) for ind in index])
    except:
        return False


def make_all_rel_vdms_hbond(path_to_csv, outpath, comb, seq_distance=5):
    if outpath[-1] != '/':
        outpath += '/'
    an = Analyze(path_to_csv)
    dist_vdms = an.get_distant_vdms(seq_distance=seq_distance)
    hbond_vdms = an.ifg_hbond_vdm[an.ifg_hbond_vdm['rel_resnums'].str.contains('0')]
    hbond_vdms_sc = hbond_vdms[hbond_vdms.apply(get_sc_hb, axis=1)]
    hbond_vdms_bb = hbond_vdms[hbond_vdms.apply(get_bb_hb, axis=1)]
    df_dist_vdms_sc = pd.merge(dist_vdms, hbond_vdms_sc, on=['iFG_count', 'vdM_count'])
    df_dist_vdms_bb = pd.merge(dist_vdms, hbond_vdms_bb, on=['iFG_count', 'vdM_count'])
    dist_vdms_sc = get_vdms_sc(df_dist_vdms_sc)
    make_rel_vdms(dist_vdms_sc, an, outpath + 'SC', comb, 'CA', 'N', 'C')
    dist_vdms_bb_N_CA = get_vdms_bb_N_CA(df_dist_vdms_bb)
    make_rel_vdms(dist_vdms_bb_N_CA, an, outpath + 'N_CA', comb, 'N', 'H', 'CA')
    dist_vdms_bb_C_O = get_vdms_bb_C_O(df_dist_vdms_bb)
    make_rel_vdms(dist_vdms_bb_C_O, an, outpath + 'C_O', comb, 'C', 'O', 'CA')
    dist_vdms_bb_phi_psi = get_vdms_bb_phi_psi(df_dist_vdms_bb)
    for phi_low, phi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
        for psi_low, psi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
            outpath_phipsi = outpath + 'PHI_PSI/' + str(phi_low) + '_' + str(phi_high) + '_' + str(psi_low) \
                     + '_' + str(psi_high)
            def phi_psi_select(row):
                if type(row['sec_struct_phi_psi_vdm']) is float:
                    return False
                if 'None' in row['sec_struct_phi_psi_vdm']:
                    return False
                spltrow = row['sec_struct_phi_psi_vdm'].split()
                if (phi_low <= float(spltrow[0]) < phi_high) and (psi_low <= float(spltrow[1]) < psi_high):
                    return True
                else:
                    return False
            dist_vdms_phipsi = dist_vdms_bb_phi_psi[dist_vdms_bb_phi_psi.apply(phi_psi_select, axis=1)]
            make_rel_vdms(dist_vdms_phipsi, an, outpath_phipsi, comb, 'CA', 'N', 'C')


def make_all_rel_vdms_hbond_partial_ifg(path_to_csv, outpath, comb, ifg_name=None, seq_distance=5):
    if outpath[-1] != '/':
        outpath += '/'
    an = Analyze(path_to_csv)
    dist_vdms = an.get_distant_vdms(seq_distance=seq_distance)
    if ifg_name:
        def func(row):
            # if (ifg_name in row['iFG_atom_names']) and ('0' in row['rel_resnums']):
            #     return True
            # else:
            #     return False

            if (ifg_name in row['iFG_atom_names']) and ('0' in row['rel_resnums']):
                nil_count = row['rel_resnums'].count('0')
                nil_index = row['rel_resnums'].replace('-', '').index('0')
                tests = []
                names_list = row['iFG_atom_names'].split(') (')
                try:
                    for i in range(nil_index, nil_index + nil_count):
                        tests.append(ifg_name in names_list[i])
                except:
                    print('row[iFG_atom_names]=', row['iFG_atom_names'],
                          'row[rel_resnums]=', row['rel_resnums'])
                    tests = [False]
                if any(tests):
                    return True
                else:
                    return False
            else:
                return False
        hbond_vdms = an.ifg_hbond_vdm[an.ifg_hbond_vdm.apply(func, axis=1)]
    else:
        hbond_vdms = an.ifg_hbond_vdm[an.ifg_hbond_vdm['rel_resnums'].str.contains('0')]
    hbond_vdms_sc = hbond_vdms[hbond_vdms.apply(get_sc_hb, axis=1)]
    hbond_vdms_bb = hbond_vdms[hbond_vdms.apply(get_bb_hb, axis=1)]
    df_dist_vdms_sc = pd.merge(dist_vdms, hbond_vdms_sc, on=['iFG_count', 'vdM_count'])
    df_dist_vdms_bb = pd.merge(dist_vdms, hbond_vdms_bb, on=['iFG_count', 'vdM_count'])
    dist_vdms_sc = get_vdms_sc(df_dist_vdms_sc)
    make_rel_vdms(dist_vdms_sc, an, outpath + 'SC', comb, 'CA', 'N', 'C')
    dist_vdms_bb_N_CA = get_vdms_bb_N_CA(df_dist_vdms_bb)
    make_rel_vdms(dist_vdms_bb_N_CA, an, outpath + 'N_CA', comb, 'N', 'H', 'CA')
    dist_vdms_bb_C_O = get_vdms_bb_C_O(df_dist_vdms_bb)
    make_rel_vdms(dist_vdms_bb_C_O, an, outpath + 'C_O', comb, 'C', 'O', 'CA')
    dist_vdms_bb_phi_psi = get_vdms_bb_phi_psi(df_dist_vdms_bb)
    for phi_low, phi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
        for psi_low, psi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
            outpath_phipsi = outpath + 'PHI_PSI/' + str(phi_low) + '_' + str(phi_high) + '_' + str(psi_low) \
                     + '_' + str(psi_high)
            def phi_psi_select(row):
                if type(row['sec_struct_phi_psi_vdm']) is float:
                    return False
                if 'None' in row['sec_struct_phi_psi_vdm']:
                    return False
                spltrow = row['sec_struct_phi_psi_vdm'].split()
                if (phi_low <= float(spltrow[0]) < phi_high) and (psi_low <= float(spltrow[1]) < psi_high):
                    return True
                else:
                    return False
            dist_vdms_phipsi = dist_vdms_bb_phi_psi[dist_vdms_bb_phi_psi.apply(phi_psi_select, axis=1)]
            make_rel_vdms(dist_vdms_phipsi, an, outpath_phipsi, comb, 'CA', 'N', 'C')


def make_all_interactamers_hbond_partial_ifg(path_to_csv, outpath, comb, ifg_name=None, seq_distance=5):
    if outpath[-1] != '/':
        outpath += '/'
    an = Analyze(path_to_csv)
    dist_vdms = an.get_distant_vdms(seq_distance=seq_distance)
    if ifg_name:
        def func(row):
            if (ifg_name in row['iFG_atom_names']) and ('0' in row['rel_resnums']):
                return True
            else:
                return False
        hbond_vdms = an.ifg_hbond_vdm[an.ifg_hbond_vdm.apply(func, axis=1)]
    else:
        hbond_vdms = an.ifg_hbond_vdm[an.ifg_hbond_vdm['rel_resnums'].str.contains('0')]
    hbond_vdms_sc = hbond_vdms[hbond_vdms.apply(get_sc_hb, axis=1)]
    hbond_vdms_bb = hbond_vdms[hbond_vdms.apply(get_bb_hb, axis=1)]
    df_dist_vdms_sc = pd.merge(dist_vdms, hbond_vdms_sc, on=['iFG_count', 'vdM_count'])
    # df_dist_vdms_bb = pd.merge(dist_vdms, hbond_vdms_bb, on=['iFG_count', 'vdM_count'])
    dist_vdms_sc = get_vdms_sc(df_dist_vdms_sc)
    make_interactamers(dist_vdms_sc, an, outpath + 'SC', comb)


def writepdb_and_get_pickle_info(pdb_path, pklout, vdm, comb, origin_atom, plane_atom1, plane_atom2, resn, ifg_count, vdm_count):
    '''helper function for make_rel_vdms and make_interactamers.
    Need vdmelem1dists to choose which vdm flipping to keep.'''

    
    if resn != 'TYR' and resn not in flip_residues:
        bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = \
            make_rel_vdm_coords(vdm, comb, origin_atom, plane_atom1, plane_atom2)
    
    elif resn == 'TYR':
        x = make_rel_vdm_coords(vdm, comb, origin_atom, plane_atom1, plane_atom2)
        y = make_rel_vdm_coords(vdm, comb, origin_atom, 'CE2', plane_atom2, unflipped=False)
        xelem1dist = x[8]
        yelem1dist = y[8]
        if xelem1dist <= yelem1dist:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = x
        else:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = y


    elif resn in flip_residues:
        x = make_rel_vdm_coords(vdm, comb, origin_atom, plane_atom1, plane_atom2)
        y = make_rel_vdm_coords(vdm, comb, origin_atom, plane_atom2, plane_atom1, unflipped=False)
        xelem1dist = x[8]
        yelem1dist = y[8]
        if xelem1dist <= yelem1dist:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = x
        else:
            bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs, vdmelem1avgdists = y

    print('at line 280')
    for i, bbc, scc, sccon, scccs, ic, icon, iccs, pdb in zip(range(len(bb_coords)), bb_coords, sc_coords, \
            sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs):
        
        pklout.append([ifg_count, vdm_count, bbc, scc, sccon, scccs, ic, icon, iccs, resn])
        string = repr(pdb).split()[1].split('_')
        pr.writePDB(pdb_path + '_'.join(string[:-1]) + '_'+string[-1] + '_oriented.pdb.gz', pdb)
    return pklout

def make_rel_vdms(df, an, outpath, comb, origin_atom, plane_atom1, plane_atom2):
    if outpath[-1] != '/':
        outpath += '/'
    resns = set(df.groupby('resname_vdm').groups).intersection(set(interactamer_atoms.keys()))
    for resn in resns:
        vdms = parse_interactamers_aa(df, an, resn)
        if vdms:
            pdb_path = outpath + 'pdbs/' + resn + '/'
            try:
                os.makedirs(pdb_path)
            except:
                pass
            picklepath = outpath + 'pickle/'
            try:
                os.makedirs(picklepath)
            except:
                pass
            pklout = []

            print('Making relative vdMs for ' + resn + '...')
            for vdm in vdms:
                try:
                    ifg_count = int(repr(vdm).split('_')[1])
                    vdm_count = int(repr(vdm).split('_')[3])
                    # if resn == 'PRO':
                    #     bb_coords, sc_coords, ifg_coords, pdbs = make_rel_vdm_coords(vdm, comb, origin_atom,
                    #                                                                  'CD', plane_atom2)
                    # else:
                    
                    pklout=writepdb_and_get_pickle_info(pdb_path,pklout, vdm, comb, origin_atom, plane_atom1, plane_atom2, resn, ifg_count, vdm_count)
                    
                except Exception:
                    traceback.print_exc()
            # output format = [ifg_count, vdm_count,  bb_coords, sc_coords_all, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, vdm_resn]

            if pklout:
                with open(picklepath + resn + '_rel_vdms.pickle', 'wb') as f:
                    pickle.dump(np.array(pklout, dtype=object), f)
            else:
                shutil.rmtree(pdb_path)






def make_interactamers(df, an, outpath, comb):
    if outpath[-1] != '/':
        outpath += '/'
    resns = set(df.groupby('resname_vdm').groups).intersection(set(interactamer_atoms.keys()))
    for resn in resns:
        for key in interactamer_atoms[resn].keys():
            origin_atom, plane_atom1, plane_atom2 = interactamer_atoms[resn][key]
            vdms = parse_interactamers_aa(df, an, resn)
            if vdms:
                pdb_path = outpath + 'pdbs/' + resn + '/' + key +'/'
                try:
                    os.makedirs(pdb_path)
                except:
                    pass
                picklepath = outpath + 'pickle/'
                try:
                    os.makedirs(picklepath)
                except:
                    pass
                pklout = []

                print('Making interactamer vdMs for ' + resn + ', ' + key + '...')
                for vdm in vdms:
                    try:
                        ifg_count = int(repr(vdm).split('_')[1])
                        vdm_count = int(repr(vdm).split('_')[3])
                        # if resn == 'PRO':
                        #     bb_coords, sc_coords, ifg_coords, pdbs = make_rel_vdm_coords(vdm, comb, origin_atom,
                        #                                                                  'CD', plane_atom2)
                        # else:

                        ### added by sophia for flipping vdms ###
                        pklout=writepdb_and_get_pickle_info(pdb_path,pklout, vdm, comb, origin_atom, plane_atom1, plane_atom2, resn, ifg_count, vdm_count)

                    except Exception:
                        traceback.print_exc()
                # output format = [ifg_count, vdm_count, ifg_flip,vdm_flip, bb_coords, sc_coords_all, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, vdm_resn]

                if pklout:
                    with open(picklepath + resn + '_rel_vdms.pickle', 'wb') as f:
                        pickle.dump(np.array(pklout, dtype=object), f)
                else:
                    shutil.rmtree(pdb_path)


def get_vdms_sc(df):
    def sc(row):
        if set(row['atom_names_vdm'].split()) - {'N', 'C', 'CA', 'O', 'OXT'} != set():
            return True
        else:
            return False
    return df[df.apply(sc, axis=1)]


def get_vdms_bb_N_CA(df):
    def bb(row):
        if any(name in row['atom_names_vdm'].split() for name in ['N', 'CA']):
            if all(name not in row['atom_names_vdm'].split() for name in ['C', 'O']):
                return True
            else:
                return False
        else:
            return False
    return df[df.apply(bb, axis=1)]

def get_vdms_bb(df, NCAorCO):
    # added by Sophia. NCAorCO can be 'NCA' or 'CO'
    def bb(row):
        if NCAorCO=='NCA':
            if any(name in row['atom_names_vdm'].split() for name in ['N', 'CA']):
                return True
            else:
                return False
        elif NCAorCO=='CO':
            if any(name in row['atom_names_vdm'].split() for name in ['C', 'O']):
                return True
            else:
                return False
        else:
            return False
    return df[df.apply(bb, axis=1)]

def get_vdms_bb_C_O(df):
    def bb(row):
        if any(name in row['atom_names_vdm'].split() for name in ['C', 'O']):
            if 'N' not in row['atom_names_vdm'].split():
                return True
            else:
                return False
        else:
            return False
    return df[df.apply(bb, axis=1)]


def get_vdms_bb_phi_psi(df):
    def bb(row):
        if 'N' in row['atom_names_vdm'].split():
            if any(name in row['atom_names_vdm'].split() for name in ['C', 'O']):
                return True
            else:
                return False
        else:
            return False
    return df[df.apply(bb, axis=1)]


def parse_interactamers_aa(df, an, resn):
    #use dist_vdms as dataframe (df)
    df = df[df['resname_vdm'] == resn]
    return an.parse_vdms(df)


def make_rel_vdm_coords(pdb, comb, origin_atom, plane_atom1, plane_atom2, unflipped=True, scoring=False, vdmselection=None, ifgselection=None,ifgatoms=None):
    if vdmselection is None:
        vdmselection = 'chain X and resnum 10 and name '
    origin_coords = pdb.select(vdmselection + origin_atom).getCoords()[0]
    pdb_coords = pdb.getCoords()
    pdb_coords_neworigin = pdb_coords - origin_coords
    plane_atom1_coords = pdb.select(vdmselection + plane_atom1).getCoords()[0] - origin_coords
    plane_atom2_coords = pdb.select(vdmselection + plane_atom2).getCoords()[0] - origin_coords
    x_norm = plane_atom1_coords / np.linalg.norm(plane_atom1_coords)
    orthvec = np.cross(plane_atom1_coords, plane_atom2_coords)
    z_norm = -1 * orthvec / np.linalg.norm(orthvec)
    orthvec2 = np.cross(plane_atom1_coords, orthvec)
    y_norm = orthvec2 / np.linalg.norm(orthvec2)
    R = np.array([x_norm, y_norm, z_norm])
    pdb_coords_neworigin_rot = np.dot(pdb_coords_neworigin, R.T)
    pdbcopy = pdb.copy()
    pdbcopy.setCoords(pdb_coords_neworigin_rot)
    new_origin_coords = pdbcopy.select(vdmselection + origin_atom).getCoords()[0]
    new_plane_atom1_coords = pdbcopy.select(vdmselection + plane_atom1).getCoords()[0]
    new_plane_atom2_coords = pdbcopy.select(vdmselection + plane_atom2).getCoords()[0]
    vdm_coords = [y for x in [new_origin_coords, new_plane_atom1_coords, new_plane_atom2_coords] for y in x]
    # format of vdm_coords is [0,0,0, atom1x, atom1y, atom1z, atom2x, atom2y, atom2z]
    bb_coords = np.array([new_origin_coords, new_plane_atom1_coords, new_plane_atom2_coords])
    if pdbcopy.select(vdmselection +'CA').getResnames()[0] != 'GLY':
        # sc_coords = pdbcopy.select('chain X and resnum 10 and sidechain and not element H D').getCoords()
        sc_coords = make_sc_atom_coords(pdbcopy, vdmselection)
        sc_coords_ON = make_sc_atom_coords_ON(pdbcopy, vdmselection)
        sc_coords_CS = make_sc_atom_coords_CS(pdbcopy, vdmselection)
    else:
        sc_coords = None
        sc_coords_ON = None
        sc_coords_CS = None
    ifg_sels = make_ifg_atom_sele(pdbcopy, comb, ifgselection, scoring=scoring, ifgatoms=ifgatoms) # there are 2 if iFG can be flipped
    ### for iFG flipping:
    ### choose only one coord set! take the one where the avg of the distances to the 3 vdm atoms is less
    if len(ifg_sels) == 2:
        ambig_ind = [ind for ind in range(len(ifg_sels[0])) if str(ifg_sels[0][ind].getCoords()) != str(ifg_sels[1][ind].getCoords())] 
        first = ifg_sels[0][ambig_ind[0]].getCoords()[0]
        second = ifg_sels[1][ambig_ind[0]].getCoords()[0]
        def chooseflip(vdm_coords, which):
            one, two, three = np.array(vdm_coords[:3]), np.array(vdm_coords[3:6]), np.array(vdm_coords[6:9])
            distances = []
            for i in [one,two,three]:
                dist = np.linalg.norm(which-i)
                distances.append(dist)
            return np.array(distances).mean()

        first_avgs = chooseflip(vdm_coords, first)
        second_avgs = chooseflip(vdm_coords, second)
        if first_avgs <= second_avgs: 
            ifg_sels = [ifg_sels[0]]
        else:
            ifg_sels = [ifg_sels[1]]
    ifg_coords = []
    ifg_CS_coords = []
    ifg_ON_coords = []
    pdbcopies = []
    orig_sel = ifg_sels[0]
    for ifg_sel in ifg_sels:
        coords = np.array([sel.getCoords() for ifg in ifg_sel for sel in ifg])

        try:
            coords_ON = [ifg.select('element O N').getCoords() for ifg in ifg_sel
                                  if ifg.select('element O N') is not None][0]
            coords_CS = [ifg.select('element C S').getCoords() for ifg in ifg_sel
                                  if ifg.select('element C S') is not None][0]
            ifg_CS_coords.append(coords_CS)
            ifg_ON_coords.append(coords_ON)
        except:
            pass
        ifg_coords.append(coords)
        for sel, coor in zip(orig_sel, coords):
            sel.setCoords(coor)
        pdbcopies.append(pdbcopy.copy())
    
    ### for vdm flipping: report avg distances of vdm flippable atoms (elements 1 and 2)
    ### to all 3 atoms of iFG
    def choosevdm(vdm_coords, ifg_coords, unflipped):
        elem1vdm = np.array(vdm_coords[3:6])
        
        dist = []
        ifg_coords = ifg_coords[0][::-1] # to get terminal 3 atoms of ifg
        for i in range(min(3,len(ifg_coords))):
            d = np.linalg.norm(ifg_coords[i] - elem1vdm)
            dist.append(d)
        return np.array(dist).mean()
    
    vdmelem1avgdists = choosevdm(vdm_coords,ifg_coords, unflipped)
    return [bb_coords]*len(ifg_coords), [sc_coords]*len(ifg_coords), [sc_coords_ON]*len(ifg_coords), \
           [sc_coords_CS]*len(ifg_coords), ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbcopies, vdmelem1avgdists

def make_sc_atom_coords(pdb,vdmselection=None):
    '''Added by Sophia: need to flip vdms if necessary!'''
    resname = pdb.select(vdmselection+'CA').getResnames()[0]
    sc_coords = []
    for name in residue_sc_names[resname]:
        coord = pdb.select(vdmselection + name).getCoords()[0]
        sc_coords.append(coord)
    sc_coords = np.array(sc_coords)
    return sc_coords

def make_sc_atom_coords_ON(pdb, vdmselection=None):
    resname = pdb.select(vdmselection+'CA').getResnames()[0]
    sc_coords = np.array([pdb.select(vdmselection + name).getCoords()[0]
                          for name in residue_sc_names[resname] if name[0] in {'N', 'O'}])
    return sc_coords

def make_sc_atom_coords_CS(pdb, vdmselection=None):
    resname = pdb.select(vdmselection+'CA').getResnames()[0]
    sc_coords = np.array([pdb.select(vdmselection + name).getCoords()[0]
                          for name in residue_sc_names[resname] if name[0] in {'C', 'S'}])
    return sc_coords

# This code is problematic because what if the iFG has a PHE or TYR in it, where more than 2 atoms need to be flipped at once?
# Sophia: don't need to worry about this now bc I'm only looking at the interactamer atoms
def make_ifg_atom_sele(pdb, comb, ifgselection=None, scoring=False, ifgatoms=None):
    """uses iFG definitions in comb object to select iFGs in the parsed protein object that have all atoms
    and occupancies = 1.
    """
    # There is a problem with this code: What if one wants to select atoms from a HEME as an iFG?
    # It is not represented in the one_letter_code dictionary...
    ifgs = []
    if comb.num_res_ifg == 1:
        poss_ifg_sel = pdb.select('chain Y and resnum 10')
        if scoring==True:
            poss_ifg_sel = pdb.select('chain X and resnum 1')
        if ifgselection is not None:
            poss_ifg_sel = pdb.select(ifgselection)
        if poss_ifg_sel is not None:
            ifg_resindices, indices = np.unique(poss_ifg_sel.getResindices(), return_index=True)
            ifg_resnames = poss_ifg_sel.getResnames()[indices]
            for ifg_resindex, ifg_resname in zip(ifg_resindices, ifg_resnames):
                    try:    
                        ifg_selection = pdb.select('resindex ' + str(ifg_resindex) + ' and name '
                                                          + comb.ifg_sele_dict[1][ifg_resname])
                        if scoring==True and ifgatoms is not None and ifg_selection is None and comb.ifg_sele_dict[1][ifg_resname]=='C O':
                            # hack for backboneCO, bc otherwise will return None
                            raise Exception
                    except:
                        ifg_selection = pdb.select('resindex ' + str(ifg_resindex) + ' and name '
                                                          + (' ').join(ifgatoms))
                                                    
                    if ifg_selection is not None:
                        num_atoms = len(ifg_selection)
                        name_list = comb.ifg_sele_dict[1][ifg_resname].split()
                        name_set = set(name_list)
                        if num_atoms == len(name_set):
                            if scoring==True and ifgatoms is not None and comb.ifg_sele_dict[1][ifg_resname]=='C O':
                                # hack for backboneCO, bc otherwise will return None
                                ifgs.append([ifg_selection.select('name ' + name) for name in ifgatoms])
                            elif all(ifg_selection.getResnums() > 0):
                                ifgs.append([ifg_selection.select('name ' + name) for name in name_list])
                            for flip_set in flip_sets:
                                if flip_set.issubset(name_set):
                                    new_names = np.array(name_list)
                                    flip_names = list(flip_set)
                                    i, = np.where(new_names == flip_names[0])
                                    j, = np.where(new_names == flip_names[1])
                                    temp = new_names[i]
                                    new_names[i] = new_names[j]
                                    new_names[j] = temp
                                    ifgs.append([ifg_selection.select('name ' + name) for name in new_names])
        return ifgs
    else:
        poss_ifg_sel = pdb.select('chain Y and resnum 10to11')
        if poss_ifg_sel is not None:
            ifg_resindices_cat_list, indices = np.unique(poss_ifg_sel.getResindices(), return_index=True)
            ifg_resnames_cat_list = poss_ifg_sel.getResnames()[indices]
            ifg_resindex_pairs = [ifg_resindices_cat_list[i:i + 2] for i in range(0, len(ifg_resindices_cat_list), 2)]
            ifg_resname_pairs = [ifg_resnames_cat_list[i:i + 2] for i in range(0, len(ifg_resnames_cat_list), 2)]
            for ifg_resindex_pair, ifg_resname_pair in zip(ifg_resindex_pairs, ifg_resname_pairs):
                resind1, resind2 = ifg_resindex_pair
                resname1, resname2 = ifg_resname_pair
                try:
                    ifg_selection1 = pdb.select('(resindex ' + str(resind1) + ' and name '
                                                      + comb.ifg_sele_dict[1][resname1]+')'
                                                      )
                    ifg_selection2 = pdb.select('(resindex ' + str(resind2) + ' and name '
                                                + comb.ifg_sele_dict[2][resname2] + ')')
                except KeyError:
                    print('Non-canonical residue in iFG, skipping.')
                    ifg_selection1 = None
                if ifg_selection1 is not None and ifg_selection2 is not None:
                    num_atoms = len(ifg_selection1 | ifg_selection2)
                    name_list1 = comb.ifg_sele_dict[1][resname1].split()
                    name_list2 = comb.ifg_sele_dict[2][resname2].split()
                    name_set1 = set(name_list1)
                    name_set2 = set(name_list2)
                    if num_atoms == (len(name_set1) + len(name_set2)):
                        if all(ifg_selection1.getResnums() > 0) and all(ifg_selection2.getResnums() > 0):
                            temp_ifg = [ifg_selection1.select('name ' + name) for name in name_list1]
                            temp_ifg.extend(ifg_selection2.select('name ' + name) for name in name_list2)
                            ifgs.append(temp_ifg)
                            flip1 = False
                            flip2 = False
                            for flip_set in flip_sets:
                                if flip_set.issubset(name_list1):
                                    flip1 = True
                                    new_names1 = np.array(name_list1)
                                    flip_names = list(flip_set)
                                    i, = np.where(new_names1 == flip_names[0])
                                    j, = np.where(new_names1 == flip_names[1])
                                    temp = new_names1[i]
                                    new_names1[i] = new_names1[j]
                                    new_names1[j] = temp
                                    temp_ifg = [ifg_selection1.select('name ' + name) for name in new_names1]
                                    temp_ifg.extend(ifg_selection2.select('name ' + name) for name in name_list2)
                                    ifgs.append(temp_ifg)
                                if flip_set.issubset(name_list2):
                                    flip2 = True
                                    new_names2 = np.array(name_list2)
                                    flip_names = list(flip_set)
                                    i, = np.where(new_names2 == flip_names[0])
                                    j, = np.where(new_names2 == flip_names[1])
                                    temp = new_names1[i]
                                    new_names1[i] = new_names1[j]
                                    new_names1[j] = temp
                                    temp_ifg = [ifg_selection1.select('name ' + name) for name in name_list1]
                                    temp_ifg.extend(ifg_selection2.select('name ' + name) for name in new_names2)
                                    ifgs.append(temp_ifg)
                            if flip1 and flip2:
                                temp_ifg = [ifg_selection1.select('name ' + name) for name in new_names1]
                                temp_ifg.extend(ifg_selection2.select('name ' + name) for name in new_names2)
                                ifgs.append(temp_ifg)
        return ifgs



### SOPHIA ADDED ###
def has_bb_or_sc(row, bb_or_sc,threepfive=False):
    # determine whether there are bb or sc vdmatoms
    # threepfive is an option for calculating freq_aai where you
    # only want to get the bb & sc vdMs if the interacting atoms
    # are within 3.5A
    bb_atoms = ['N', 'CA', 'C', 'O', 'OXT']
    try:
        vdmatoms = row['atom_names_vdm'].split(' ')
    except:
        vdmatoms = row['atom_names'].split(' ')
    num_bb_in_vdm = sum([x in bb_atoms for x in vdmatoms])
    has_bb = num_bb_in_vdm > 0
    has_sc = len(vdmatoms) > num_bb_in_vdm
    bb_bool = np.array([x in bb_atoms for x in vdmatoms]) 
    if bb_or_sc=='bb':
        if threepfive==True:
            close = within_threepfive(row, bb_atoms,bbsc='bb')
            return has_bb and close
        return has_bb 
    elif bb_or_sc=='sc':
        if threepfive==True:
            close = within_threepfive(row, bb_atoms,bbsc='sc')
            return has_sc and close
        return has_sc 

def within_threepfive(row,bb_atoms,bbsc):
    # returns whether or not the interacting atom is within 3.5A
    dist = row['dist_info']
    dist = dist.strip('()').split(') (')
    dist = [x.split(' ') for x in dist] 
    # this is now a list of lists where the inner list is 
    # element0: iFG, element1: vdM, element2: distance
    # check if the vdMs are bb or sc atoms
    if bbsc=='bb':
        vdms = [v for v in dist if v[1] in bb_atoms]
    elif bbsc=='sc':
        vdms = [v for v in dist if v[1] not in bb_atoms]
    distances = [float(v[2]) for v in vdms]
    numclose = sum([d<=3.5 for d in distances])
    return numclose > 0

def make_bb_sc_rel_vdms(ifg_dir, comb,lonepair_imidazole=False):
    ''''makes transformed PDB files and pkl files for clustering'''
    path_to_csv = ifg_dir+'/csv/'
    outpath = ifg_dir+'/clusters/'
    dist_vdms, vdms_bb, vdms_sc, an = Analysis.refine_df(path_to_csv, seq_dist=10, threefive=True,lonepair_imidazole=lonepair_imidazole)

    print(len(vdms_sc), 'len vdms_sc within 3.5A')
    make_interactamers(vdms_sc, an, outpath + 'SC', comb)

    dist_vdms_bb_N_CA = get_vdms_bb(vdms_bb,NCAorCO='NCA')
    dist_vdms_bb_C_O = get_vdms_bb(vdms_bb,NCAorCO='CO')
    print(len(dist_vdms_bb_N_CA), 'N CA')
    print(len(dist_vdms_bb_C_O), 'C O')
    make_rel_vdms(dist_vdms_bb_N_CA, an, outpath + 'N_CA', comb, 'N', 'H', 'CA')
    make_rel_vdms(dist_vdms_bb_C_O, an, outpath + 'C_O', comb, 'C', 'O', 'CA')


