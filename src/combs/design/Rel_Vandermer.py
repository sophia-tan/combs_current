__all__ = ['RelVandermer']

import pickle
import os
import numpy as np
from ..apps.transformation import get_rot_trans
# from scipy.spatial.distance import cdist
from sklearn.neighbors import NearestNeighbors
import collections
from ..apps.constants import bb_type_dict
import networkx as nx
import itertools
import functools
import traceback
import prody as pr
import random
from ..apps.functions import writePDBStream

class RelVandermer:
    def __init__(self, name):
        """name: name of the iFG of the vandermers"""
        self.name = name
        self.rel_vdm_path = None
        self._rois_rot_trans = collections.defaultdict(dict)
        self.rel_vdms_pickle = collections.defaultdict(dict)
        self.rel_vdm_ifg_coords = collections.defaultdict(dict)
        self.rel_vdm_ifg_ON_coords = collections.defaultdict(dict)
        self.rel_vdm_ifg_CS_coords = collections.defaultdict(dict)
        self.rel_vdm_bb_coords = {}
        self.rel_vdm_sc_coords = {}
        self.rel_vdm_sc_ON_coords = {}
        self.rel_vdm_sc_CS_coords = {}
        self.rel_vdm_phipsi_bin = {}
        self.rel_vdm_tags = collections.defaultdict(dict)
        self.rel_vdm_resnames = collections.defaultdict(dict)
        self.hotspots = None
        self.hotspot_graph = None
        self.hotspot_subgraphs = None
        self._ifgs = {}
        self._scs = {}
        self._resn = {}
        self._type = {}
        self._sc_ONs = {}
        self._sc_CSs = {}
        self._vdm_tags = {}
        self._indices = {}
        self._resn_sc = {}
        self._type_sc = {}
        self._vdm_tags_sc = {}
        self._indices_sc = {}
        self.ifg_dict = None
        self.hotspot_ifg_coords = collections.defaultdict(dict)
        self._hs_ifgs = None
        self._hs_nums = None

    def load_rel_vdms_pickle(self, sample, subset=None):
        """Loads rel_vdms at every residue position of the binding site.
        Is this really necessary at every position?"""
        for resnum_chid in sample.rois.keys():
            self.rel_vdms_pickle[resnum_chid] = collections.defaultdict(dict)
            for type_ in bb_type_dict.keys():
                if type_ == 'PHI_PSI':
                    if sample.rois_phipsi[resnum_chid][0] and sample.rois_phipsi[resnum_chid][1]:
                        for phi_low, phi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
                            for psi_low, psi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
                                if (phi_low <= sample.rois_phipsi[resnum_chid][0] < phi_high) \
                                 and (psi_low <= sample.rois_phipsi[resnum_chid][1] < psi_high):
                                    phipsis = [phi_low, phi_high, psi_low, psi_high]
                                    phipsis_str = '_'.join([str(x) for x in phipsis])
                                    type_ = 'PHI_PSI/' + phipsis_str
                                    self.rel_vdm_phipsi_bin[resnum_chid] = phipsis_str
                        picklepath = self.rel_vdm_path + type_ + '/pickle/'

                        if os.path.isdir(picklepath):
                            for picklefile in os.listdir(picklepath):
                                if subset:
                                    if picklefile.split('_')[0] in subset:
                                        with open(picklepath + picklefile, 'rb') as infile:
                                            self.rel_vdms_pickle[resnum_chid][type_.split('/')[0]][picklefile.split('_')[0]] \
                                                = pickle.load(infile)
                                else:
                                    with open(picklepath + picklefile, 'rb') as infile:
                                        self.rel_vdms_pickle[resnum_chid][type_.split('/')[0]][picklefile.split('_')[0]] \
                                            = pickle.load(infile)
                else:
                    picklepath = self.rel_vdm_path + type_ + '/pickle/'
                    if os.path.isdir(picklepath):
                        for picklefile in os.listdir(picklepath):
                            if subset:
                                if picklefile.split('_')[0] in subset:
                                    with open(picklepath + picklefile, 'rb') as infile:
                                        self.rel_vdms_pickle[resnum_chid][type_.split('/')[0]][picklefile.split('_')[0]] \
                                            = pickle.load(infile)
                            else:
                                with open(picklepath + picklefile, 'rb') as infile:
                                    self.rel_vdms_pickle[resnum_chid][type_.split('/')[0]][picklefile.split('_')[0]] \
                                        = pickle.load(infile)

    def set_rel_vdm_bb_coords(self):
        """Stores the backbone coordinates of the first vdm for each residue position
        in the bs site, vdm type, and residue name."""
        for key_ in self.rel_vdms_pickle.keys():
            self.rel_vdm_bb_coords[key_] = collections.defaultdict(dict)
            for type_ in self.rel_vdms_pickle[key_].keys():
                for resn in self.rel_vdms_pickle[key_][type_].keys():
                    # if len(self.rel_vdms_pickle[key_][type_][resn]) > 1:
                    # print(key_, type_, resn, self.rel_vdms_pickle[key_][type_][resn])
                    self.rel_vdm_bb_coords[key_][type_][resn] = self.rel_vdms_pickle[key_][type_][resn][0, 3]
                    # else:
                    #     print(key_, type_, resn, self.rel_vdms_pickle[key_][type_][resn])
                    #     self.rel_vdm_bb_coords[key_][type_][resn] = self.rel_vdms_pickle[key_][type_][resn][3]

    def set_rois_rot_trans(self, sample):
        """will store tuple of (rotation, rel_vdm_bb_com, roi_bb_com) in dictionary"""
        for resnum_chid, roi in sample.rois.items():
            self._rois_rot_trans[resnum_chid] = collections.defaultdict(dict)
            for type_ in self.rel_vdms_pickle[resnum_chid].keys():
                for resn in self.rel_vdms_pickle[resnum_chid][type_].keys():
                    self._rois_rot_trans[resnum_chid][type_][resn] = \
                        get_rot_trans(self.rel_vdm_bb_coords[resnum_chid][type_][resn],
                                      sample.rois_bb_coords[resnum_chid][type_])

    def move_rel_vdms(self, sample):
        """moves the ifg and sidechain coordinates at each roi per vdm type and residue name"""
        for resnum_chid in sample.rois.keys():
            self.rel_vdm_ifg_coords[resnum_chid] = collections.defaultdict(dict)
            self.rel_vdm_ifg_ON_coords[resnum_chid] = collections.defaultdict(dict)
            self.rel_vdm_ifg_CS_coords[resnum_chid] = collections.defaultdict(dict)
            self.rel_vdm_sc_coords[resnum_chid] = collections.defaultdict(dict)
            self.rel_vdm_sc_ON_coords[resnum_chid] = collections.defaultdict(dict)
            self.rel_vdm_sc_CS_coords[resnum_chid] = collections.defaultdict(dict)

            for type_ in self.rel_vdms_pickle[resnum_chid].keys():
                for resn in self.rel_vdms_pickle[resnum_chid][type_].keys():

                    coords = np.array([coord for coord in self.rel_vdms_pickle[resnum_chid][type_][resn][:, -4]])
                    self.rel_vdm_ifg_coords[resnum_chid][type_][resn] = \
                      np.dot((coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                             self._rois_rot_trans[resnum_chid][type_][resn][0]) \
                      + self._rois_rot_trans[resnum_chid][type_][resn][2]

                    coords = np.array([coord for coord in self.rel_vdms_pickle[resnum_chid][type_][resn][:, -3]])
                    if coords[0].shape[0] != 0:
                        self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn] = \
                            np.dot((coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                                   self._rois_rot_trans[resnum_chid][type_][resn][0]) \
                            + self._rois_rot_trans[resnum_chid][type_][resn][2]

                    coords = np.array([coord for coord in self.rel_vdms_pickle[resnum_chid][type_][resn][:, -2]])
                    if coords[0].shape[0] != 0:
                        self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn] = \
                            np.dot((coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                                   self._rois_rot_trans[resnum_chid][type_][resn][0]) \
                            + self._rois_rot_trans[resnum_chid][type_][resn][2]

                    if type_ == 'SC':

                        coords = np.array([coord for coord in self.rel_vdms_pickle[resnum_chid][type_][resn][:, -7]])
                        self.rel_vdm_sc_coords[resnum_chid][type_][resn] = \
                            np.dot((coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                                   self._rois_rot_trans[resnum_chid][type_][resn][0]) \
                            + self._rois_rot_trans[resnum_chid][type_][resn][2]

                        coords = np.array([coord for coord in self.rel_vdms_pickle[resnum_chid][type_][resn][:, -6]])
                        if coords[0].shape[0] != 0:
                            self.rel_vdm_sc_ON_coords[resnum_chid][type_][resn] = \
                                np.dot((coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                                       self._rois_rot_trans[resnum_chid][type_][resn][0]) \
                                + self._rois_rot_trans[resnum_chid][type_][resn][2]

                        coords = np.array([coord for coord in self.rel_vdms_pickle[resnum_chid][type_][resn][:, -5]])
                        if coords[0].shape[0] != 0:
                            self.rel_vdm_sc_CS_coords[resnum_chid][type_][resn] = \
                                np.dot((coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                                       self._rois_rot_trans[resnum_chid][type_][resn][0]) \
                                + self._rois_rot_trans[resnum_chid][type_][resn][2]

    # @staticmethod
    # def calc_clash_score(poi_clash_coords, test_coords, clash_cutoff=2.5):
    #     """calculates clash score via a distance matrix with a cutoff. test_coords is any coordinate array,
    #     such as an ifg or vdm sidechain"""
    #     shape = test_coords.shape
    #     return np.sum(np.sum(cdist(poi_clash_coords, test_coords.flatten().reshape(shape[0] * shape[1],
    #         shape[2])).T.reshape(shape[0], shape[1], len(poi_clash_coords)) < clash_cutoff, axis=1), axis=1)

    @staticmethod
    def calc_clashing(poi_clash_coords, test_coords, cutoff=2.5):
        """calculates clash score via a distance matrix with a cutoff. test_coords is any coordinate array,
        such as an ifg or vdm sidechain"""
        shape = test_coords.shape
        distances, indices = poi_clash_coords.kneighbors(test_coords.flatten().reshape(shape[0] * shape[1], shape[2]))
        cut_distances = distances < cutoff #poi_bb_coords.radius
        return ~cut_distances.reshape(shape[0], shape[1]).any(axis=1)

    def set_rel_vdm_tags(self, sample):
        """Stores the vdm tags per bs residue, vdm type, residue name"""
        for resnum_chid in sample.rois.keys():
            self.rel_vdm_tags[resnum_chid] = collections.defaultdict(dict)
            for type_ in self.rel_vdms_pickle[resnum_chid].keys():
                for resn in self.rel_vdms_pickle[resnum_chid][type_].keys():
                    self.rel_vdm_tags[resnum_chid][type_][resn] = \
                        self.rel_vdms_pickle[resnum_chid][type_][resn][:, :3].astype(int)

    # def remove_clash(self, sample):
    #     """removes vdms with iFG or sidechain clashes with the backbone of the protein of interest"""
    #     for resnum_chid in sample.rois.keys():
    #         for resn in self.rel_vdms_pickle[resnum_chid]['SC'].keys():
    #             unique_rows = self.calc_clashing(sample.poi_clash_coords[resnum_chid],
    #                                                    self.rel_vdm_sc_coords[resnum_chid]['SC'][resn])
    #             self.rel_vdm_sc_coords[resnum_chid]['SC'][resn] = \
    #                 self.rel_vdm_sc_coords[resnum_chid]['SC'][resn][unique_rows]
    #             self.rel_vdm_ifg_coords[resnum_chid]['SC'][resn] = \
    #                 self.rel_vdm_ifg_coords[resnum_chid]['SC'][resn][unique_rows]
    #             self.rel_vdm_tags[resnum_chid]['SC'][resn] = self.rel_vdm_tags[resnum_chid]['SC'][resn][unique_rows]
    #
    #         for type_ in self.rel_vdms_pickle[resnum_chid].keys():
    #             for resn in self.rel_vdms_pickle[resnum_chid][type_].keys():
    #                 if self.rel_vdm_ifg_coords[resnum_chid][type_][resn].shape[0] > 0:
    #                     unique_rows = self.calc_clashing(sample.poi_clash_coords[resnum_chid],
    #                                                            self.rel_vdm_ifg_coords[resnum_chid][type_][resn])
    #                     self.rel_vdm_ifg_coords[resnum_chid][type_][resn] = \
    #                         self.rel_vdm_ifg_coords[resnum_chid][type_][resn][unique_rows]
    #                     self.rel_vdm_tags[resnum_chid][type_][resn] = \
    #                         self.rel_vdm_tags[resnum_chid][type_][resn][unique_rows]
    #                     if type_ == 'SC':
    #                         self.rel_vdm_sc_coords[resnum_chid]['SC'][resn] = \
    #                             self.rel_vdm_sc_coords[resnum_chid]['SC'][resn][unique_rows]

    def remove_clash(self, sample, cutoff_NO_NO=2.5, cutoff_NO_C=2.8, cutoff_C_C=3.2):
        """removes vdms with iFG or sidechain clashes with the backbone of the protein of interest"""
        for resnum_chid in sample.rois.keys():
            for resn in self.rel_vdms_pickle[resnum_chid]['SC'].keys():
                try:
                    sc_ON_coor = self.rel_vdm_sc_ON_coords[resnum_chid]['SC'][resn]
                    unique_rows_ON_ON = self.calc_clashing(sample.poi_clash_coords_NO[resnum_chid],
                                                           sc_ON_coor, cutoff_NO_NO)
                    unique_rows_C_ON = self.calc_clashing(sample.poi_clash_coords_C[resnum_chid],
                                                          sc_ON_coor, cutoff_NO_C)
                except:
                    unique_rows_ON_ON = True
                    unique_rows_C_ON = True

                try:
                    sc_CS_coor = self.rel_vdm_sc_CS_coords[resnum_chid]['SC'][resn]
                    unique_rows_ON_CS = self.calc_clashing(sample.poi_clash_coords_NO[resnum_chid],
                                                       sc_CS_coor, cutoff_NO_C)
                    unique_rows_C_CS = self.calc_clashing(sample.poi_clash_coords_C[resnum_chid],
                                                      sc_CS_coor, cutoff_C_C)
                except:
                    unique_rows_ON_CS = True
                    unique_rows_C_CS = True

                unique_rows = unique_rows_ON_ON * unique_rows_ON_CS * unique_rows_C_ON * unique_rows_C_CS

                self.rel_vdm_sc_coords[resnum_chid]['SC'][resn] = \
                    self.rel_vdm_sc_coords[resnum_chid]['SC'][resn][unique_rows]
                if resn in self.rel_vdm_sc_ON_coords[resnum_chid]['SC']:
                    self.rel_vdm_sc_ON_coords[resnum_chid]['SC'][resn] = \
                        self.rel_vdm_sc_ON_coords[resnum_chid]['SC'][resn][unique_rows]
                if resn in self.rel_vdm_sc_CS_coords[resnum_chid]['SC']:
                    self.rel_vdm_sc_CS_coords[resnum_chid]['SC'][resn] = \
                        self.rel_vdm_sc_CS_coords[resnum_chid]['SC'][resn][unique_rows]
                self.rel_vdm_ifg_coords[resnum_chid]['SC'][resn] = \
                    self.rel_vdm_ifg_coords[resnum_chid]['SC'][resn][unique_rows]
                if resn in self.rel_vdm_ifg_ON_coords[resnum_chid]['SC']:
                    self.rel_vdm_ifg_ON_coords[resnum_chid]['SC'][resn] = \
                        self.rel_vdm_ifg_ON_coords[resnum_chid]['SC'][resn][unique_rows]
                if resn in self.rel_vdm_ifg_CS_coords[resnum_chid]['SC']:
                    self.rel_vdm_ifg_CS_coords[resnum_chid]['SC'][resn] = \
                        self.rel_vdm_ifg_CS_coords[resnum_chid]['SC'][resn][unique_rows]
                self.rel_vdm_tags[resnum_chid]['SC'][resn] = self.rel_vdm_tags[resnum_chid]['SC'][resn][unique_rows]

            for type_ in self.rel_vdms_pickle[resnum_chid].keys():
                for resn_ in self.rel_vdms_pickle[resnum_chid][type_].keys():
                    t1 = self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn_].shape[0]
                    t2 = self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn_].shape[0]

                    if t1 > 0 and t2 > 0:
                        unique_rows_ON_ON = self.calc_clashing(sample.poi_clash_coords_NO[resnum_chid],
                                                         self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn_], cutoff_NO_NO)
                        unique_rows_C_ON = self.calc_clashing(sample.poi_clash_coords_C[resnum_chid],
                                                               self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn_],
                                                               cutoff_NO_C)
                        unique_rows_ON_CS = self.calc_clashing(sample.poi_clash_coords_NO[resnum_chid],
                                                               self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn_],
                                                               cutoff_NO_C)
                        unique_rows_C_CS = self.calc_clashing(sample.poi_clash_coords_C[resnum_chid],
                                                              self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn_],
                                                              cutoff_C_C)
                        unique_rows = unique_rows_ON_ON * unique_rows_ON_CS * unique_rows_C_ON * unique_rows_C_CS
                        self.rel_vdm_ifg_coords[resnum_chid][type_][resn_] = \
                            self.rel_vdm_ifg_coords[resnum_chid][type_][resn_][unique_rows]
                        if resn_ in self.rel_vdm_ifg_ON_coords[resnum_chid][type_]:
                            self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn_] = \
                                self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn_][unique_rows]
                        if resn_ in self.rel_vdm_ifg_CS_coords[resnum_chid][type_]:
                            self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn_] = \
                                self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn_][unique_rows]
                        self.rel_vdm_tags[resnum_chid][type_][resn_] = self.rel_vdm_tags[resnum_chid][type_][resn_][
                            unique_rows]
                        if type_ == 'SC':
                            self.rel_vdm_sc_coords[resnum_chid]['SC'][resn_] = \
                                self.rel_vdm_sc_coords[resnum_chid]['SC'][resn_][unique_rows]
                            if resn_ in self.rel_vdm_sc_ON_coords[resnum_chid]['SC']:
                                self.rel_vdm_sc_ON_coords[resnum_chid]['SC'][resn_] = \
                                    self.rel_vdm_sc_ON_coords[resnum_chid]['SC'][resn_][unique_rows]
                            if resn_ in self.rel_vdm_sc_CS_coords[resnum_chid]['SC']:
                                self.rel_vdm_sc_CS_coords[resnum_chid]['SC'][resn_] = \
                                    self.rel_vdm_sc_CS_coords[resnum_chid]['SC'][resn_][unique_rows]

                    elif t1 > 0:
                        unique_rows_ON_ON = self.calc_clashing(sample.poi_clash_coords_NO[resnum_chid],
                                                               self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn_],
                                                               cutoff_NO_NO)
                        unique_rows_C_ON = self.calc_clashing(sample.poi_clash_coords_C[resnum_chid],
                                                              self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn_],
                                                              cutoff_NO_C)
                        unique_rows = unique_rows_ON_ON * unique_rows_C_ON
                        self.rel_vdm_ifg_coords[resnum_chid][type_][resn_] = \
                            self.rel_vdm_ifg_coords[resnum_chid][type_][resn_][unique_rows]
                        self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn_] = \
                            self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn_][unique_rows]
                        self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn_] = \
                            self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn_][unique_rows]
                        self.rel_vdm_tags[resnum_chid][type_][resn_] = self.rel_vdm_tags[resnum_chid][type_][resn_][
                            unique_rows]
                        if type_ == 'SC':
                            self.rel_vdm_sc_coords[resnum_chid]['SC'][resn_] = \
                                self.rel_vdm_sc_coords[resnum_chid]['SC'][resn_][unique_rows]
                            self.rel_vdm_sc_ON_coords[resnum_chid]['SC'][resn_] = \
                                self.rel_vdm_sc_ON_coords[resnum_chid]['SC'][resn_][unique_rows]
                            self.rel_vdm_sc_CS_coords[resnum_chid]['SC'][resn_] = \
                                self.rel_vdm_sc_CS_coords[resnum_chid]['SC'][resn_][unique_rows]

                    elif t2 > 0:
                        unique_rows_ON_CS = self.calc_clashing(sample.poi_clash_coords_NO[resnum_chid],
                                                               self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn],
                                                               cutoff_NO_C)
                        unique_rows_C_CS = self.calc_clashing(sample.poi_clash_coords_C[resnum_chid],
                                                              self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn],
                                                              cutoff_C_C)
                        unique_rows = unique_rows_ON_CS * unique_rows_C_CS
                        self.rel_vdm_ifg_coords[resnum_chid][type_][resn_] = \
                            self.rel_vdm_ifg_coords[resnum_chid][type_][resn_][unique_rows]
                        self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn_] = \
                            self.rel_vdm_ifg_ON_coords[resnum_chid][type_][resn_][unique_rows]
                        self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn_] = \
                            self.rel_vdm_ifg_CS_coords[resnum_chid][type_][resn_][unique_rows]
                        self.rel_vdm_tags[resnum_chid][type_][resn_] = self.rel_vdm_tags[resnum_chid][type_][resn_][
                            unique_rows]
                        if type_ == 'SC':
                            self.rel_vdm_sc_coords[resnum_chid]['SC'][resn_] = \
                                self.rel_vdm_sc_coords[resnum_chid]['SC'][resn_][unique_rows]
                            self.rel_vdm_sc_ON_coords[resnum_chid]['SC'][resn_] = \
                                self.rel_vdm_sc_ON_coords[resnum_chid]['SC'][resn_][unique_rows]
                            self.rel_vdm_sc_CS_coords[resnum_chid]['SC'][resn_] = \
                                self.rel_vdm_sc_CS_coords[resnum_chid]['SC'][resn_][unique_rows]

    def reshape_ifgs(self):
        """make all ifgs as vectors with 3N column dimension. Must keep track of everything by index.
        This is needed for the graph search for hotspots."""
        for resnum_chid in self.rel_vdms_pickle.keys():
            _ifgs = []
            _resn = []
            _type = []
            _vdm_tags = []
            _indices = []
            _scs = []
            _sc_ONs = []
            _sc_CSs = []
            _resn_sc = []
            _type_sc = []
            _vdm_tags_sc = []
            _indices_sc = []
            for type_ in self.rel_vdms_pickle[resnum_chid].keys():
                for resn in self.rel_vdms_pickle[resnum_chid][type_].keys():
                    _ifgs.append(self.rel_vdm_ifg_coords[resnum_chid][type_][resn])
                    # if type_ == 'SC':
                    #     _scs.append(self.rel_vdm_sc_coords[resnum_chid][type_][resn])
                    #     _sc_ONs.append(self.rel_vdm_sc_ON_coords[resnum_chid][type_][resn])
                    #     _sc_CSs.append(self.rel_vdm_sc_CS_coords[resnum_chid][type_][resn])
                    #     len_sc_ = len(self.rel_vdm_sc_coords[resnum_chid][type_][resn])
                    #     _resn_sc.append([resn] * len_sc_)
                    #     _type_sc.append([type_] * len_sc_)
                    #     _vdm_tags_sc.append(self.rel_vdm_tags[resnum_chid][type_][resn])
                    #     _indices_sc.append(list(range(len_sc_)))
                    len_ = len(self.rel_vdm_ifg_coords[resnum_chid][type_][resn])
                    _resn.append([resn]*len_)
                    _type.append([type_]*len_)
                    _vdm_tags.append(self.rel_vdm_tags[resnum_chid][type_][resn])
                    _indices.append(list(range(len_)))
            # if _scs:
            #     _scs = np.array(functools.reduce(lambda a, b: np.vstack((a, b)), _scs))
            #     _sc_ONs = np.array(functools.reduce(lambda a, b: np.vstack((a, b)), _sc_ONs))
            #     _sc_CSs = np.array(functools.reduce(lambda a, b: np.vstack((a, b)), _sc_CSs))
            #     _resn_sc = np.array(functools.reduce(lambda a, b: np.hstack((a, b)), _resn_sc))
            #     _type_sc = np.array(functools.reduce(lambda a, b: np.hstack((a, b)), _type_sc))
            #     _vdm_tags_sc = np.array(functools.reduce(lambda a, b: np.vstack((a, b)), _vdm_tags_sc))
            #     _indices_sc = np.array(functools.reduce(lambda a, b: np.hstack((a, b)), _indices_sc), dtype=int)
            #     self._scs[resnum_chid] = _scs
            #     self._sc_ONs[resnum_chid] = _sc_ONs
            #     self._sc_CSs[resnum_chid] = _sc_CSs
            #     self._resn_sc[resnum_chid] = _resn_sc.reshape(-1)
            #     self._type_sc[resnum_chid] = _type_sc.reshape(-1)
            #     self._vdm_tags_sc[resnum_chid] = _vdm_tags_sc
            #     self._indices_sc[resnum_chid] = _indices_sc
            _ifgs = np.array(functools.reduce(lambda a, b: np.vstack((a, b)), _ifgs))
            _resn = np.array(functools.reduce(lambda a, b: np.hstack((a, b)), _resn))
            _type = np.array(functools.reduce(lambda a, b: np.hstack((a, b)), _type))
            _vdm_tags = np.array(functools.reduce(lambda a, b: np.vstack((a, b)), _vdm_tags))
            _indices = np.array(functools.reduce(lambda a, b: np.hstack((a, b)), _indices), dtype=int)
            _shape = _ifgs.shape
            self._ifgs[resnum_chid] = _ifgs.reshape(_shape[0], int(np.prod(_shape[1:])))
            self._resn[resnum_chid] = _resn.reshape(-1)
            self._type[resnum_chid] = _type.reshape(-1)
            self._vdm_tags[resnum_chid] = _vdm_tags
            self._indices[resnum_chid] = _indices

    # def find_hotspots(self, rmsd_cutoff=1.0):
    #     """hotspots are comprised of vdms of different bs positions that pairwise place an iFG in space
    #     within the tolerance of rmsd_cutoff. hotspots are more akin to contiguous densities.  they
    #     should not be confused with what would result from greedy clustering of iFGs in space with the same
    #     rmsd cutoff."""
    #     g = nx.Graph()
    #     nbrs = collections.defaultdict(dict)
    #     for i, j in itertools.combinations(self.rel_vdms_pickle.keys(), 2):
    #         try:
    #             if not nbrs[j]:
    #                 nbrs[j] = NearestNeighbors(metric='euclidean', radius=rmsd_cutoff)
    #                 nbrs[j].fit(self._ifgs[j])
    #             adj_mat = nbrs[j].radius_neighbors_graph(self._ifgs[i])
    #             wh = adj_mat.nonzero()
    #             zip_ = list(zip(wh[0], wh[1]))
    #             nodes, edges = self.sort_pairs(zip_, i, j)
    #             g.add_nodes_from(nodes)
    #             g.add_edges_from(edges)
    #             # nx.set_node_attributes(g, 'ifg', node_attr_ifg)
    #         except Exception:
    #             traceback.print_exc()
    #     self.hotspot_graph = g
    #     self.hotspots = list(sorted(nx.connected_components(g), key=len, reverse=True))
    #
    # @staticmethod
    # def sort_pairs(zip_, i, j):
    #     nodes = []
    #     edges = []
    #     for item in zip_:
    #         sorted_ = sorted(item)
    #         if sorted_ in zip_:
    #             nodes.extend([(i, sorted_[0]), (j, sorted_[1])])
    #             edges.append(((i, sorted_[0]), (j, sorted_[1])))
    #         else:
    #             nodes.extend([(i, item[0]), (j, item[1])])
    #             edges.append(((i, item[0]), (j, item[1])))
    #     return nodes, edges

    def find_hotspots(self, rmsd_cutoff=1.0):
        """hotspots are comprised of vdms of different bs positions that pairwise place an iFG in space
        within the tolerance of rmsd_cutoff. hotspots are more akin to contiguous densities.  they
        should not be confused with what would result from greedy clustering of iFGs in space with the same
        rmsd cutoff."""
        g = nx.Graph()
        nbrs = collections.defaultdict(dict)
        for i, j in itertools.combinations(self.rel_vdms_pickle.keys(), 2):
            try:
                if not nbrs[j]:
                    nbrs[j] = NearestNeighbors(metric='euclidean', radius=rmsd_cutoff, n_jobs=-1)
                    nbrs[j].fit(self._ifgs[j])
                adj_mat = nbrs[j].radius_neighbors_graph(self._ifgs[i])
                wh = adj_mat.nonzero()
                zip_ = zip(wh[0], wh[1])
                sortpairs = SortPairs()
                sortpairs.sort(self, i, j, zip_)
                g.add_nodes_from(sortpairs.nodes)
                g.add_edges_from(sortpairs.edges)
                nx.set_node_attributes(g, 'ifg', sortpairs.node_attr_ifg)
                nx.set_node_attributes(g, 'resn', sortpairs.node_attr_resn)
                nx.set_node_attributes(g, 'type', sortpairs.node_attr_type)
                nx.set_node_attributes(g, 'vdm_tags', sortpairs.node_attr_vdm_tags)
                nx.set_node_attributes(g, 'sc', sortpairs.node_attr_sc)
            except Exception:
                traceback.print_exc()
        self.hotspot_graph = g
        self.hotspots = list(sorted(nx.connected_components(g), key=len, reverse=True))
        self.hotspot_subgraphs = list(sorted(nx.connected_component_subgraphs(g), key=len, reverse=True))

    def print_hotspots(self, outdir, number=25):
        "prints pdbs of the top number of hotspots (number) to the output directory (outdir)."
        if outdir[-1] != '/':
            outdir += '/'
        for i in range(number):
            for resnum_chid, index in self.hotspots[i]:
                resn = self._resn[resnum_chid][index]
                type_ = self._type[resnum_chid][index]
                typestr = type_
                vdm_tags = self._vdm_tags[resnum_chid][index]
                if typestr == 'PHI_PSI':
                    typestr = 'PHI_PSI/' + self.rel_vdm_phipsi_bin[resnum_chid]
                pdbpath = self.rel_vdm_path + typestr + '/pdbs/' + resn + '/'
                filename = 'iFG_' + str(vdm_tags[0]) + '_vdM_' + str(vdm_tags[1]) + '_iFlip_' \
                                  + str(vdm_tags[2]) + '_' + self.name + '_' + 'oriented.pdb.gz'
                pdb = pr.parsePDB(pdbpath + filename)
                old_coords = pdb.getCoords()
                new_coords = \
                    np.dot((old_coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                           self._rois_rot_trans[resnum_chid][type_][resn][0]) \
                    + self._rois_rot_trans[resnum_chid][type_][resn][2]
                pdb.setCoords(new_coords)
                newfile_path = outdir + 'hotspots/' + str(i + 1) + '/'
                newfile_name = self.name + '_hotspot_' + str(i + 1) + '_' \
                               + ''.join(str(x) for x in resnum_chid) + '_' + type_ + '_' + filename
                try:
                    os.makedirs(newfile_path)
                except:
                    pass
                pr.writePDB(newfile_path + newfile_name, pdb)

    def print_hotspot(self, hs_subgr_num, outdir):
        """prints pdbs of the top number of hotspots (number) to the output directory (outdir)."""
        if outdir[-1] != '/':
            outdir += '/'
        for node in self.hotspot_subgraphs[hs_subgr_num]:
            resnum_chid = node[0]
            type_ = self.hotspot_subgraphs[hs_subgr_num].node[node]['type']
            typestr = type_
            resn = self.hotspot_subgraphs[hs_subgr_num].node[node]['resn']
            vdm_tags = self.hotspot_subgraphs[hs_subgr_num].node[node]['vdm_tags']

            if typestr == 'PHI_PSI':
                typestr = 'PHI_PSI/' + self.rel_vdm_phipsi_bin[resnum_chid]
            pdbpath = self.rel_vdm_path + typestr + '/pdbs/' + resn + '/'
            filename = 'iFG_' + str(vdm_tags[0]) + '_vdM_' + str(vdm_tags[1]) + '_iFlip_' \
                              + str(vdm_tags[2]) + '_' + self.name + '_' + 'oriented.pdb.gz'
            pdb = pr.parsePDB(pdbpath + filename)
            old_coords = pdb.getCoords()
            new_coords = \
                np.dot((old_coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                       self._rois_rot_trans[resnum_chid][type_][resn][0]) \
                + self._rois_rot_trans[resnum_chid][type_][resn][2]
            pdb.setCoords(new_coords)
            newfile_path = outdir + 'hotspots/' + str(hs_subgr_num + 1) + '/'
            newfile_name = self.name + '_hotspot_' + str(hs_subgr_num + 1) + '_' \
                           + ''.join(str(x) for x in resnum_chid) + '_' + type_ + '_' + filename
            try:
                os.makedirs(newfile_path)
            except:
                pass
            pr.writePDB(newfile_path + newfile_name, pdb)

    def print_cluster(self, mems, label, outdir, number=None):
        """prints pdbs of the top number of hotspots (number) to the output directory (outdir)."""
        if outdir[-1] != '/':
            outdir += '/'
        if number:
            mems = random.sample(mems, number)
        for mem in mems:
            resnum_chid = tuple(self._all_resnum_chid[mem])
            type_ = self._all_type[mem]
            typestr = type_
            resn = self._all_resn[mem]
            vdm_tags = self._all_vdm_tags[mem]

            if typestr == 'PHI_PSI':
                typestr = 'PHI_PSI/' + self.rel_vdm_phipsi_bin[resnum_chid]
            pdbpath = self.rel_vdm_path + typestr + '/pdbs/' + resn + '/'
            filename = 'iFG_' + str(vdm_tags[0]) + '_vdM_' + str(vdm_tags[1]) + '_iFlip_' \
                       + str(vdm_tags[2]) + '_' + self.name + '_' + 'oriented.pdb.gz'
            pdb = pr.parsePDB(pdbpath + filename)
            old_coords = pdb.getCoords()
            new_coords = \
                np.dot((old_coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                       self._rois_rot_trans[resnum_chid][type_][resn][0]) \
                + self._rois_rot_trans[resnum_chid][type_][resn][2]
            pdb.setCoords(new_coords)
            newfile_path = outdir + 'cluster/' + label + '/'
            newfile_name = self.name + '_cluster_' + label + '_' \
                           + ''.join(str(x) for x in resnum_chid) + '_' + type_ + '_' + filename
            try:
                os.makedirs(newfile_path)
            except:
                pass
            pr.writePDB(newfile_path + newfile_name, pdb)

    def print_member(self, poi, mem, label, outdir, ligand, counter=1):
        """prints pdbs of the top number of hotspots (number) to the output directory (outdir)."""
        if outdir[-1] != '/':
            outdir += '/'

        resnum_chid = tuple(self._all_resnum_chid[mem])
        type_ = self._all_type[mem]
        typestr = type_
        resn = self._all_resn[mem]
        vdm_tags = self._all_vdm_tags[mem]

        bb = poi.select('backbone and resnum ' + str(resnum_chid[0]) + ' and chain ' + resnum_chid[1])
        bb.setChids('X')
        bb.setResnums(10)
        num_ind = len(bb.getIndices())
        start = 1
        finish = start + num_ind
        bb.setBetas(list(range(start, finish)))
        bb.setResnames(resn)

        if typestr == 'PHI_PSI':
            typestr = 'PHI_PSI/' + self.rel_vdm_phipsi_bin[resnum_chid]
        pdbpath = self.rel_vdm_path + typestr + '/pdbs/' + resn + '/'
        filename = 'iFG_' + str(vdm_tags[0]) + '_vdM_' + str(vdm_tags[1]) + '_iFlip_' \
                   + str(vdm_tags[2]) + '_' + self.name + '_' + 'oriented.pdb.gz'
        pdb = pr.parsePDB(pdbpath + filename)
        old_coords = pdb.getCoords()
        new_coords = \
            np.dot((old_coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                   self._rois_rot_trans[resnum_chid][type_][resn][0]) \
            + self._rois_rot_trans[resnum_chid][type_][resn][2]
        pdb.setCoords(new_coords)
        newfile_path = outdir
        newfile_name = str(counter) + '_mem_' + str(mem) + '_' + self.name + '_' + label + '_' \
                       + ''.join(str(x) for x in resnum_chid) + '_' + type_ + '_' + filename[:-3]

        sc = pdb.select('sidechain and chain X and resnum 10')
        num_ind = len(sc.getIndices())
        start = finish
        finish = start + num_ind
        sc.setBetas(list(range(start, finish)))

        num_ind = len(ligand.select('all').getIndices())
        start = finish
        finish = start + num_ind
        ligand.setBetas(list(range(start, finish)))

        try:
            os.makedirs(newfile_path)
        except:
            pass
        with open(newfile_path + newfile_name, 'w') as outfile:
            writePDBStream(outfile, bb)
            writePDBStream(outfile, pdb.select('sidechain and chain X and resnum 10'))
            writePDBStream(outfile, ligand)

    def _print_pair(self, poi, mem):
        resnum_chid = tuple(self._all_resnum_chid[mem])
        type_ = self._all_type[mem]
        typestr = type_
        resn = self._all_resn[mem]
        vdm_tags = self._all_vdm_tags[mem]

        bb = poi.select('backbone and resnum ' + str(resnum_chid[0]) + ' and chain ' + resnum_chid[1])

        if typestr == 'PHI_PSI':
            typestr = 'PHI_PSI/' + self.rel_vdm_phipsi_bin[resnum_chid]
        pdbpath = self.rel_vdm_path + typestr + '/pdbs/' + resn + '/'
        filename = 'iFG_' + str(vdm_tags[0]) + '_vdM_' + str(vdm_tags[1]) + '_iFlip_' \
                   + str(vdm_tags[2]) + '_' + self.name + '_' + 'oriented.pdb.gz'
        pdb = pr.parsePDB(pdbpath + filename)
        old_coords = pdb.getCoords()
        new_coords = \
            np.dot((old_coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                   self._rois_rot_trans[resnum_chid][type_][resn][0]) \
            + self._rois_rot_trans[resnum_chid][type_][resn][2]
        pdb.setCoords(new_coords)
        return bb, pdb, filename[:-3]

    def parse_member(self, mem):
        """prints pdbs of the top number of hotspots (number) to the output directory (outdir)."""

        resnum_chid = tuple(self._all_resnum_chid[mem])
        type_ = self._all_type[mem]
        typestr = type_
        resn = self._all_resn[mem]
        vdm_tags = self._all_vdm_tags[mem]


        if typestr == 'PHI_PSI':
            typestr = 'PHI_PSI/' + self.rel_vdm_phipsi_bin[resnum_chid]
        pdbpath = self.rel_vdm_path + typestr + '/pdbs/' + resn + '/'
        filename = 'iFG_' + str(vdm_tags[0]) + '_vdM_' + str(vdm_tags[1]) + '_iFlip_' \
                   + str(vdm_tags[2]) + '_' + self.name + '_' + 'oriented.pdb.gz'
        pdb = pr.parsePDB(pdbpath + filename)
        old_coords = pdb.getCoords()
        new_coords = \
            np.dot((old_coords - self._rois_rot_trans[resnum_chid][type_][resn][1]),
                   self._rois_rot_trans[resnum_chid][type_][resn][0]) \
            + self._rois_rot_trans[resnum_chid][type_][resn][2]
        pdb.setCoords(new_coords)

        return pdb, resn


    # def set_hotspot_ifg_coords(self):
    #     for hs_num, hotspot in enumerate(self.hotspots):
    #         # hotspot = sorted(hotspot, key=lambda x: x[0])
    #         # _inter = [self._ifgs[k][np.array(np.array(list(g))[:, 1], dtype=int)]
    #         #                                    for k, g in itertools.groupby(list(hotspot), key=lambda x: x[0])]
    #         # self.hotspot_ifg_coords[hs_num] = functools.reduce(lambda a, b: np.vstack((a, b)), _inter)
    #
    #         self.hotspot_ifg_coords[hs_num] = [self._ifgs[resnum_chid][index] for resnum_chid, index in hotspot]
    #     _hs_nums = np.array([[hs_num]*len(coords) for hs_num, coords in self.hotspot_ifg_coords.items()])
    #     # _hs_ifgs = np.array([coords for coords in self.hotspot_ifg_coords.values()])
    #     # self._hs_ifgs = np.array(functools.reduce(lambda a, b: np.vstack((a, b)), _hs_ifgs))
    #     self._hs_ifgs = np.array(functools.reduce(lambda a, b: np.vstack((a, b)), list(self.hotspot_ifg_coords.values())))
    #     self._hs_nums = np.array(functools.reduce(lambda a, b: np.hstack((a, b)), _hs_nums))
    #     # self._hs_ifgs = np.array()

    def make_pymol_session(self):
        pass

class SortPairs:
    """SortPairs will sort a zipped index object from the nn adjacency matrix, such that
    duplicates are removed.  It will also make attribute dictionaries to be
    added to the graph nodes in a subsequent step"""
    def __init__(self):
        self.nodes = []
        self.edges = []
        self.node_attr_ifg = {}
        self.node_attr_resn = {}
        self.node_attr_type = {}
        self.node_attr_vdm_tags = {}
        self.node_attr_sc = {}

    def sort(self, rel_vdm, i, j, zip_):
        func = functools.partial(self.sort_internal, rel_vdm, i, j)
        list(map(func, zip_))

    def sort_internal(self, rel_vdm, i, j, item):
        ind_i = item[0]
        ind_j = item[1]
        self.node_attr_ifg[(i, ind_i)] = rel_vdm._ifgs[i][ind_i]
        self.node_attr_ifg[(j, ind_j)] = rel_vdm._ifgs[j][ind_j]
        resn_i = rel_vdm._resn[i][ind_i]
        resn_j = rel_vdm._resn[j][ind_j]
        self.node_attr_resn[(i, ind_i)] = resn_i
        self.node_attr_resn[(j, ind_j)] = resn_j
        type_i = rel_vdm._type[i][ind_i]
        type_j = rel_vdm._type[j][ind_j]
        self.node_attr_type[(i, ind_i)] = type_i
        self.node_attr_type[(j, ind_j)] = type_j
        self.node_attr_vdm_tags[(i, ind_i)] = rel_vdm._vdm_tags[i][ind_i]
        self.node_attr_vdm_tags[(j, ind_j)] = rel_vdm._vdm_tags[j][ind_j]
        self.nodes.extend([(i, ind_i), (j, ind_j)])
        self.edges.append(((i, ind_i), (j, ind_j)))
        if type_i == 'SC':
            self.node_attr_sc[(i, ind_i)] = \
                rel_vdm.rel_vdm_sc_coords[i][type_i][resn_i][rel_vdm._indices[i][ind_i]]
        if type_j == 'SC':
            self.node_attr_sc[(j, ind_j)] = \
                rel_vdm.rel_vdm_sc_coords[j][type_j][resn_j][rel_vdm._indices[j][ind_j]]





