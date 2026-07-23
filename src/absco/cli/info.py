#!/usr/bin/env python

import sys
import argparse
import numpy as np
import xarray as xa

from absco._common import utils

def format_value(val):
    """Format a value with 2 decimals or exponential notation if needed."""
    if abs(val) < 1.0 or abs(val) > 100000.0:
        return f'{val:.2e}'
    else:
        return f'{val:.2f}'

def format_range(arr, unit=''):
    """Format a range with appropriate notation."""
    min_val = np.nanmin(arr)
    max_val = np.nanmax(arr)
    unit_str = f' {unit}' if unit else ''
    return f'[{format_value(min_val)}, {format_value(max_val)}]{unit_str}'

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Display information about an ABSCO netCDF file including '
                    'spectral, pressure, temperature, and water vapor ranges.')
    parser.add_argument('ncFile', type=str,
                        help='ABSCO netCDF file to read.')
    args = parser.parse_args()

    utils.file_check(args.ncFile)

    print(f'ABSCO File: {args.ncFile}\n')

    with xa.open_dataset(args.ncFile) as ds:
        # Get molecule name from filename or attributes
        import os
        mol_name = os.path.basename(args.ncFile).split('_')[0]
        print(f'Molecule: {mol_name}')

        # Check if this is an H2O-affected molecule
        h2o_affected = mol_name in ['H2O', 'HDO', 'O2', 'CO2', 'N2']
        o2_special = mol_name == 'O2'

        # Spectral range
        wn = np.array(ds['Spectral_Grid'])
        print(f'\nSpectral Range: {format_range(wn, "cm-1")}')
        print(f'  Number of points: {len(wn)}')

        # Pressure range
        p = np.array(ds['P_level'])
        print(f'\nPressure Range: {format_range(p, "mbar")}')
        print(f'  Number of levels: {len(p)}')
        print(f'  Number of layers: {len(p)-1}')

        # Temperature range
        t = np.array(ds['T_level'])
        print(f'\nTemperature Range: {format_range(t, "K")}')
        # Get unique temperatures across all levels
        unique_temps = np.unique(t[~np.isnan(t)])
        print(f'  Number of temperature points: {len(unique_temps)}')

        # Water vapor VMR range (if applicable)
        if h2o_affected:
            h2o = np.array(ds['H2O_VMR'])
            print(f'\nH2O VMR Range: {format_range(h2o, "ppmv")}')
            print(f'  Number of H2O VMR points: {len(h2o)}')

        # O2 VMR range (if applicable)
        if o2_special:
            o2 = np.array(ds['O2_VMR'])
            print(f'\nO2 VMR Range: {format_range(o2, "ppmv")}')
            print(f'  Number of O2 VMR points: {len(o2)}')

        # Cross section shape and fill value info
        xs = ds['Cross_Section']
        print(f'\nCross Section Shape: {xs.shape}')
        if hasattr(xs, '_FillValue'):
            print(f'  Fill Value: {xs._FillValue}')
        elif hasattr(xs, 'missing_value'):
            print(f'  Fill Value: {xs.missing_value}')

        # Check for any metadata attributes
        if hasattr(ds, 'title'):
            print(f'\nTitle: {ds.title}')
        if hasattr(ds, 'history'):
            print(f'History: {ds.history}')

if __name__ == '__main__':
    main()
