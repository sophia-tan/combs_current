# Comb.py equivalent for database non-combing statistics 
__all__ = ['nrPDB']

import os
import re
from contextlib import ContextDecorator
import csv
from collections import defaultdict
from ..apps.constants import one_letter_code


class nrPDB(ContextDecorator):
    """
    :arg path_to_pdb_chain_file: path to a txt file that has 5 letter pdb accession and unique chain on each line
    :type path_to_pdb_chain_file: str
    :arg file_tag: user-defined tag on the pdb and csv files that get printed
    :type file_tag: str
    :arg input_dir_pdb: path to input directory for pdbs to be combed.
    :type input_dir_pdb: str
    :arg output_dir_pdb: path to output directory for vandermer pdbs. Default is current working directory.
    :type output_dir_pdb: str
    :arg output_dir_csv: path to output directory for comb csv database files. Default is current working 
        directory.
    :type output_dir_csv: str
    :type increment: int
    :attribute pdb_chains: each unique chain associated with a pdb to be searched through to find iFGs. This is a built
    attribute that requires the path_to_pdb_chain_file and needs no other user input.
    :type pdb_chains: dict with keys=pdb accession code, values=unique chain
    :arg add_non_canonical: dictionary with keys as three letter resname strings, values as a one-letter string.
    :type add_non_canonical: dictionary
    :arg path_to_reduce: path to the reduce program
    :arg reduce: name of reduce executable
    
    """

    def __init__(self, **kwargs):
        """param kwargs: """
        _cwd = os.getcwd()

        path_to_pdb_chain_file = kwargs.get('path_to_pdb_chain_file', _cwd)
        

        try:
            self.pdbs_chains = self.make_pdb_chain_dict(path_to_pdb_chain_file)
        except:
            self.pdbs_chains = None

        self.file_tag = kwargs.get('file_tag', 'comb')
        self.input_dir_pdb = kwargs.get('input_dir_pdb', _cwd)
        if self.input_dir_pdb[-1] != '/':
            self.input_dir_pdb += '/'
        self.input_dir_dssp = kwargs.get('input_dir_dssp', _cwd)
        if self.input_dir_dssp[-1] != '/':
            self.input_dir_dssp += '/'
        self.output_dir_pdb = kwargs.get('output_dir_pdb', _cwd)
        if self.output_dir_pdb[-1] != '/':
            self.output_dir_pdb += '/'
        self.output_dir_csv = kwargs.get('output_dir_csv', _cwd)
        if self.output_dir_csv[-1] != '/':
            self.output_dir_csv += '/'
        self.path_to_reduce = kwargs.get('path_to_reduce', _cwd)
        if self.path_to_reduce[-1] != '/':
            self.path_to_reduce += '/'
        self.reduce = kwargs.get('reduce', 'reduce')


    @staticmethod
    def make_pdb_chain_dict(path_to_pdb_chain_file_):
        with open(path_to_pdb_chain_file_) as infile:
            pdb_chain_dict = defaultdict(list)
            for line in infile:
                try:
                    pdb = line[0:4].lower()
                    chain = line[4]
                    pdb_chain_dict[pdb].append(chain)
                except Exception:
                    print('This pdb was not included from txt file: ', line)
        return pdb_chain_dict

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False
