#!/usr/bin/env python

import os
import logging
from itertools import product
from argparse import ArgumentParser

import numpy as np

DEFAULT_MOLS_PER_FILE=1

from configparser import ConfigParser

logger = logging.getLogger(__file__)

class AbscoConfigSplitter(object):

    def __init__(self, config_filename, chunk_size=None, mols_per_file=DEFAULT_MOLS_PER_FILE):

        # Store so we can use it to name output files
        self.in_filename = config_filename

        self.config = ConfigParser()
        self.config.read(config_filename)

        self.chunk_size = chunk_size
        self.mols_per_file = mols_per_file

    def channel_windows(self):

        in_beg_list = np.array([float(v) for v in self.config['channels']['wn1'].split()])
        in_end_list = np.array([float(v) for v in self.config['channels']['wn2'].split()])
        lblres_list = np.array([float(v) for v in self.config['channels']['lblres'].split()])
        outres_list = np.array([float(v) for v in self.config['channels']['outres'].split()])

        out_beg_list = []
        out_end_list = []
        for beg, end, lblres, outres in zip(in_beg_list, in_end_list, lblres_list, outres_list):
            if self.chunk_size is None or self.chunk_size == 0:
                num_chunks = 1
            else:
                num_chunks = int(np.round(abs(beg-end)/self.chunk_size))

            # Existing window is small enough
            if num_chunks <= 1:
                yield beg, end, lblres, outres
            else:
                for chunk_idx in range(num_chunks):
                    # Make sure ending chunk does not go beyond initial ending
                    sub_beg = beg + self.chunk_size * chunk_idx 
                    sub_end = min(sub_beg + self.chunk_size, end)
                    yield sub_beg, sub_end, lblres, outres

    def molecule_lists(self):

        mol_names = self.config['molecules']['molnames'].split()

        if self.mols_per_file is None or self.mols_per_file == 0:
            yield tuple(mol_names)
            return

        curr_set = []
        for mol in mol_names:
            if len(curr_set) == self.mols_per_file:
                yield tuple(curr_set)
                curr_set = []
            
            curr_set.append(mol)

        if len(curr_set) > 0:
            yield tuple(curr_set)


    def write_configs(self):
        
        config_prod = product(self.channel_windows(), self.molecule_lists())
        fn_base, fn_ext = os.path.splitext(self.in_filename)
        print(self.mols_per_file)

        for idx, (channel_info, mol_names) in enumerate(config_prod):
            chan_beg, chan_end, lbl_res, out_res = channel_info
            
            output_filename = "{base}-{num:03}-{chan_beg:.0f}_{chan_end:.0f}-{molecules}{ext}".format(
                base=fn_base, ext=fn_ext, num=(idx+1), chan_beg=chan_beg, chan_end=chan_end, molecules="_".join(mol_names))

            logger.debug("Creating configuration file: {}".format(output_filename))

            # Reread original file name to create a copy for writing
            out_config = ConfigParser()
            out_config.read(self.in_filename)

            # Modify sections
            out_config['channels']['wn1'] = "{:0.8f}".format(chan_beg)
            out_config['channels']['wn2'] = "{:0.8f}".format(chan_end)
            out_config['channels']['lblres'] = "{:0.10e}".format(lbl_res)
            out_config['channels']['outres'] = "{:0.10e}".format(out_res)

            out_config['molecules']['molnames'] = " ".join(mol_names)

            # Add a subdirectory to output path so that multiple runs 
            # do not clobber each other
            out_sub_dir = "{num:02}".format(num=idx+1)
            out_config['output']['intdir'] = os.path.join(out_config['output']['intdir'], out_sub_dir)

            with open(output_filename, "w") as out_file:
                out_config.write(out_file)

def main():

    parser = ArgumentParser("Split a ABSCO configuration file into smaller pieces to make computation easier")

    parser.add_argument("config_file", 
        help="Source configuration filename")

    parser.add_argument("-c", "--chunk_size", type=float,
        help="Size of the spectral chunks in the units specified by the configuration, defaults are unit specific")

    parser.add_argument("-m", "--molecules_per_file", default=DEFAULT_MOLS_PER_FILE, type=int,
        help="Number of molecules per each file, default: {}. Set to 0 to use all".format(DEFAULT_MOLS_PER_FILE))
    
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG, format="%(message)s")

    splitter = AbscoConfigSplitter(args.config_file, args.chunk_size, args.molecules_per_file)
    splitter.write_configs()

if __name__ == "__main__":
    main()
