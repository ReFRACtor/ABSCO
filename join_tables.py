#!/usr/bin/env python

import os
import logging
from argparse import ArgumentParser

import netCDF4
import numpy as np

logger = logging.getLogger(__file__)

# Datasets to check for consistency between files
DIMS_ENSURE_SAME = ("nlay", "nlev", "nranges_lims", "ntemp", "nvmr")
VARIABLES_ENSURE_SAME = ("/H2O_VMR", "/P_layer", "/P_level", "/Temperature")

# How many units of the cross section table to copy at a time to reduce memory overhead
# Number here ~ 4GB
XSECT_CHUNK_SIZE = 180000 

def copy_attrs(src_var, dst_var):
    for attr_k, attr_v in src_var.__dict__.items():
        # Skip attributes starting with _, such as _FillValue
        if attr_k[0] != "_":
            setattr(dst_var, attr_k, attr_v)

class AbscoFileJoiner(object):

    def __init__(self, inp_filenames):

        # Open the netCDF4 object for each file
        self.inp_objects = []
        for fn in inp_filenames:
            if not os.path.exists(fn):
                raise Exception("Input ABSCO file not found:", fn)

            inp_nc_file = netCDF4.Dataset(fn, "r") 
            self.inp_objects.append(inp_nc_file)

        self._sort_tables()

    def _sort_tables(self):
        "Sort open table objects based on extent ranges"

        def by_extent(obj):
            return obj['Extent_Ranges'][0, 0]

        self.inp_objects = sorted(self.inp_objects, key=by_extent)

    def gas_name(self):

        all_gn = [ os.path.basename(o.filepath()).split("_")[0] for o in self.inp_objects ]

        gn0 = all_gn[0]
        for curr_gn in all_gn[1:]:
            if curr_gn != gn0:
                raise Exception("Not all input files have the same gas name: {gn0} {curr_gn}".format(gn0, curr_gn))

        return gn0

    def check_properties(self):

        logger.debug("Checking for consistency among input ABSCO tables")

        base_obj = self.inp_objects[0]
        for comp_obj in self.inp_objects[1:]:
            # Check dimensions
            for dim_name in DIMS_ENSURE_SAME:
                # Some dimensions not present in all files, such as nvmr
                if dim_name not in base_obj.dimensions and dim_name not in comp_obj.dimensions:
                    continue

                base_size = base_obj.dimensions[dim_name].size 
                comp_size = comp_obj.dimensions[dim_name].size
                if base_size != comp_size:
                    raise Exception("Dimension {dim_name} does not match ({d1} vs {d2}) between files {fn1} and {fn2}".format(
                        dim_name=dim_name, d1=base_size, d2=comp_size, fn1=base_obj.filepath(), fn2=comp_obj.filepath()))

            # Check variables using the masked array comparison routine
            for var_name in VARIABLES_ENSURE_SAME:
                # H2O_VMR not in all files
                if not var_name in base_obj.variables and not var_name in comp_obj.variables:
                    continue

                if not np.ma.allequal(base_obj[var_name], comp_obj[var_name], fill_value=True):
                    raise Exception("Dataset {var_name} is not the same between files {fn1} and {fn2}".format(
                        var_name=var_name, fn1=base_obj.filepath(), fn2=comp_obj.filepath()))

    def total_extent(self):
        return np.min(self.inp_objects[0]['Extent_Ranges']), np.max(self.inp_objects[-1]['Extent_Ranges'])

    def write(self, output_filename=None):

        if output_filename is None:
            _, ext = os.path.splitext(self.inp_objects[0].filepath())
            extent_1, extent_2 = self.total_extent()
            output_filename = "{gas_name}_{extent_1:.0f}-{extent_2:.0f}_joined{ext}".format(
                gas_name=self.gas_name(), extent_1=extent_1, extent_2=extent_2, ext=ext)

        logger.debug("Creating joined ABSCO table file: {}".format(output_filename))
        
        output_fil = netCDF4.Dataset(output_filename, "w", zlib=True, complevel=4)

        self._write_dimensions(output_fil)
        self._write_common_variables(output_fil)
        self._write_updated_variables(output_fil)

    def _write_dimensions(self, output_fil):
        # Copy the dimensions and variables that are the same from the first file
        logger.debug("Copying unchanged dimensions")
        for dim_name in DIMS_ENSURE_SAME:
            if dim_name not in self.inp_objects[0].dimensions:
                continue

            dim_obj = self.inp_objects[0].dimensions[dim_name]
            output_fil.createDimension(dim_name, dim_obj.size)

        # Create updated dimensions, accumulate the new sizes
        nranges_val = 0
        nfreq_val = 0
        for curr_fil in self.inp_objects:
            nranges_val += curr_fil.dimensions["nranges"].size
            nfreq_val += curr_fil.dimensions["nfreq"].size

        logger.debug("Creating changed dimensions: nranges = {}, nfreq = {}".format(nranges_val, nfreq_val))
        nranges_dim = output_fil.createDimension("nranges", nranges_val)
        nfreq_dim = output_fil.createDimension("nfreq", nfreq_val)

    def _write_common_variables(self, output_fil):
        inp_file = self.inp_objects[0]
        
        # Turn off auto masking to avoid having issue when copying masked values
        inp_file.set_auto_mask(False)

        # Copy common variables from the first object
        for var_name in VARIABLES_ENSURE_SAME:
            if var_name not in inp_file.variables:
                continue

            logger.debug("Copying unchanged variable: {}".format(var_name))
            src_var = inp_file[var_name]
            dst_var = output_fil.createVariable(var_name, src_var.dtype, src_var.dimensions, fill_value=src_var._FillValue)

            dst_var[:] = src_var[:]

            copy_attrs(src_var, dst_var)

        inp_file.set_auto_mask(True)

    def _write_updated_variables(self, output_fil):
        # Create new extent and spectral grid information
        logger.debug("Creating updated variables: Spectral_Grid, Extent_Ranges, Extent_Indices, Cross_Section")

        dst_dtype = self.inp_objects[0]["Spectral_Grid"].dtype
        dst_spec_grid = output_fil.createVariable("Spectral_Grid", dst_dtype, ("nfreq",), fill_value=np.nan)
        copy_attrs(self.inp_objects[0]["Spectral_Grid"], dst_spec_grid)

        dst_dtype = self.inp_objects[0]["Extent_Ranges"].dtype
        dst_extent_ranges = output_fil.createVariable("Extent_Ranges", dst_dtype, ("nranges", "nranges_lims"), fill_value=np.nan)
        copy_attrs(self.inp_objects[0]["Extent_Ranges"], dst_extent_ranges)

        dst_dtype = self.inp_objects[0]["Extent_Indices"].dtype
        dst_extent_indicies = output_fil.createVariable("Extent_Indices", dst_dtype, ("nranges", "nranges_lims"))
        copy_attrs(self.inp_objects[0]["Extent_Indices"], dst_extent_indicies)

        dst_dtype = self.inp_objects[0]["Cross_Section"].dtype

        # Some files are not broadened by another gas
        if "nvmr" in output_fil.dimensions:
            cross_section_dims = ("nfreq", "ntemp", "nlay", "nvmr")
        else:
            cross_section_dims = ("nfreq", "ntemp", "nlay")

        dst_cross_section = output_fil.createVariable("Cross_Section", dst_dtype, cross_section_dims, fill_value=np.nan)
        copy_attrs(self.inp_objects[0]["Cross_Section"], dst_cross_section)

        logger.debug("Copying values for updated variables")
        dst_freq_beg = 0
        dst_ranges_beg = 0
        for curr_fil in self.inp_objects:
            # Turn off auto masking to avoid having issue when copying masked values
            curr_fil.set_auto_mask(False)

            logger.debug("Copying from {}".format(curr_fil.filepath()))

            src_nfreq = curr_fil.dimensions["nfreq"].size
            src_nranges = curr_fil.dimensions["nranges"].size

            src_spec_grid = curr_fil["Spectral_Grid"]
            src_extent_ranges = curr_fil["Extent_Ranges"]
            src_extent_indicies = curr_fil["Extent_Indices"]
            src_cross_section = curr_fil["Cross_Section"]
    
            # These refer to spectral points and can be copied over as is
            dst_spec_grid[dst_freq_beg:dst_freq_beg+src_nfreq] = src_spec_grid[:]
            dst_extent_ranges[dst_ranges_beg:dst_ranges_beg+src_nranges, :] = src_extent_ranges[:, :]

            # Update indicies based on the range sizes of the input
            dst_range_idx = dst_freq_beg
            for range_num in range(src_nranges):
                range_nfreq = src_extent_indicies[range_num, 1] - src_extent_indicies[range_num, 0] + 1
                
                dst_extent_indicies[dst_ranges_beg+range_num, 0] = dst_range_idx
                dst_extent_indicies[dst_ranges_beg+range_num, 1] = dst_range_idx + range_nfreq - 1

                dst_range_idx += range_nfreq

            # Break into chunks to avoid exhausting memory
            dst_chunk_beg = dst_freq_beg
            src_chunk_beg = 0
            chunk_num = 1
            while dst_chunk_beg < (dst_freq_beg+src_nfreq):
                dst_chunk_end = min(dst_chunk_beg + XSECT_CHUNK_SIZE, dst_freq_beg+src_nfreq)
                src_chunk_end = min(src_chunk_beg + XSECT_CHUNK_SIZE, src_nfreq)

                logger.debug(".. Cross_Section #{} {}-{} -> {}-{}".format(chunk_num, src_chunk_beg, src_chunk_end, dst_chunk_beg, dst_chunk_end))
                dst_cross_section[dst_chunk_beg:dst_chunk_end, ...] = src_cross_section[src_chunk_beg:src_chunk_end , ...]

                dst_chunk_beg += XSECT_CHUNK_SIZE
                src_chunk_beg += XSECT_CHUNK_SIZE
                chunk_num += 1

            dst_freq_beg += src_nfreq
            dst_ranges_beg += src_nranges

            curr_fil.set_auto_mask(True)

def main():

    parser = ArgumentParser("Join multiple ABSCO tables for the same molecule type together their spectral axis. All tables must have been generated with the same configuration parameters.")

    parser.add_argument("inp_filenames", metavar="FILENAME", nargs="+", 
        help="Input table filenames to be combined, order does not matter")

    parser.add_argument("-o", "--out_filename", metavar="FILENAME",
        help="Optional output filename, otherwise files are named with common prefix with output extent range")

    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG, format="%(message)s")

    if len(args.inp_filenames) <= 1:
        parser.error("Not enough input ABSCO table filenames specified, must provide more than one")

    joiner = AbscoFileJoiner(args.inp_filenames)
    joiner.check_properties()
    joiner.write(args.out_filename)

if __name__ == "__main__":
    main()
