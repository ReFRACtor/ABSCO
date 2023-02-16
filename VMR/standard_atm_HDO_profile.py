#!/usr/bin/env python

USSCSV = 'USS_AIRS_profile.csv'
HDOCSV = 'HDO_example_profile.csv'
SCALE = 3.107e-4

def makeProfHDO(inCSV=USSCSV, scale=SCALE, outCSV=HDOCSV):
  """
  HDO profile generation
  """

  import pandas as PD

  uss = PD.read_csv(inCSV)
  hdo = PD.DataFrame(uss['H2O'] * scale).rename(columns={'H2O': 'HDO'})
  hdo.to_csv(outCSV)

  print('Wrote {}'.format(outCSV))
# end makeProfHDO()

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(
    description='Scale a US standard atmosphere with a ' + \
    'non-standard HDO scale factor and write CSV that can be ' + \
    'used ABSCO_preprocess.py and ABSCO_compute.py',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--uss_csv', '-i', type=str, default=USSCSV, 
    help='US Standard Atmosphere CSV')
  parser.add_argument('--hdo_scale', '-s', type=float, default=SCALE, 
    help='Scale factor used with H2O profile in uss_csv to ' + \
    'compute corresponding HDO profile')
  parser.add_argument('--out_csv', '-o', type=str, default=HDOCSV,
    help='Name of CSV to which HDO profile is written')
  args = parser.parse_args()

  makeProfHDO(
    inCSV=args.uss_csv, scale=args.hdo_scale, outCSV=args.out_csv)
# endif main()
