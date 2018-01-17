# equivalent to Vandermer.py

__all__ = ['Res']

import prody as pr
import numpy as np
import gzip
import freesasa
from ..rotamer.rotamer import calc_rotamer
from ..apps.renumber import renumber_chids_resnums
from ..apps.constants import one_letter_code


class Res:
    def __init__(self, parsed_pdb):
        self.resnum = parsed_pdb.getResnums()[0]
        self.chid = parsed_pdb.getChids()[0]
        self.resname = parsed_pdb.getResnames()[0]
        self.sele = parsed_pdb
        
        self.residue_sasa = None
        self.dssp_sasa = None
        self.sasa_3A_probe = None
        self.sasa_4A_probe = None
        self.sasa_5A_probe = None
        self.bb_cb_atom_ind = self.get_bb_cb_atom_indices()

    def get_bb_cb_atom_indices(self):
        sele = self.sele.select('protein and (backbone or name CB) and resnum `' + str(self.resnum)
                                                          + '` and not element H D')
        if sele is not None:
            return sele.getIndices()
        else:
            return None

    def calc_large_probe_sasa(self, parsed_pdb, fs_result):
        if self.bb_cb_atom_ind.any():
            return '{0:.2f}'.format(sum(fs_result.atomArea(i) for i in np.where(np.in1d(parsed_pdb.prody_pdb_bb_cb_atom_ind,
                                                                   self.bb_cb_atom_ind))[0]))

    def get_dssp_sasa(self):
        ri, _ind = np.unique(self.sele.select('resnum `' + str(self.resnum) + '`').getResindices(), return_index=True)
        sasa = self.sele.getData('dssp_acc')[_ind]
        self.dssp_sasa = ' '.join('{0:.2f}'.format(res) for res in sasa)

    def calc_sasa(self, parsed_pdb):
        """Calculates the per atom solvent accessible surface area of the iFG and the sasa of the residue containing
        the iFG.  Needs FreeSASA module to be imported.  Takes as argument an instance of ParsedPDB class, which
        contains the iFG.  Right now this function isn't optimized, in the sense that the iFG atoms must be in the
        same residue.  Need better general way to select iFG atoms...

        parsed_pdb: an instance of class ParsedPDB having attributes .contacts, .fs_struct, .fs_result, .dssp,
        .prody_pdb
        """

        # assert isinstance(self.pdb_name, str), 'pdb name of iFG is not defined'
        # assert self.pdb_name + '.pdb' in os.listdir(path_to_pdbs), 'pdb file is not in directory'
        assert isinstance(parsed_pdb.fs_struct,
                          freesasa.Structure), 'parsed_pdb object must have attribute freesasa structure obj'
        assert isinstance(parsed_pdb.fs_result,
                          freesasa.Result), 'parsed_pdb object must have attribute freesasa result obj'

        selections = freesasa.selectArea(('vdm_residue, chain ' + self.chid + ' and resi ' + str(self.resnum),),
                                         parsed_pdb.fs_struct, parsed_pdb.fs_result)

        self.residue_sasa = '{0:.2f}'.format(selections['vdm_residue'])
        self.sasa_3A_probe = self.calc_large_probe_sasa(parsed_pdb, parsed_pdb.fs_result_cb_3A)
        self.sasa_4A_probe = self.calc_large_probe_sasa(parsed_pdb, parsed_pdb.fs_result_cb_4A)
        self.sasa_5A_probe = self.calc_large_probe_sasa(parsed_pdb, parsed_pdb.fs_result_cb_5A)


    def get_info(self, parsed_pdb):
        self.get_dssp_sasa()
        self.calc_sasa(parsed_pdb)

        #_vdm_sasa_info = [ifg.count,
        #                  ifg.vdm_count,
        #                  self.residue_sasa,
        #                  self.dssp_sasa,
        #                  self.sasa_3A_probe,
        #                  self.sasa_4A_probe,
        #                  self.sasa_5A_probe
        #                  ]

        #parsed_pdb._vdm_sasa_info.append(_vdm_sasa_info)
