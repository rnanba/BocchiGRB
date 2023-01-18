#!/usr/bin/env python
import sys
import csv
from argparse import ArgumentParser
from math import log10
from astropy import units as u
from astropy.cosmology import LambdaCDM, z_at_value

cosmo = LambdaCDM(H0=67.3, Om0=0.315, Ode0=0.685)
bocchi_star_lyr = 100000000.0

def get_bocchi_grb_mag(grb):
    name = grb['name']
    z = grb['z']
    mag = grb['mag']
    # print(name, z, mag)
    bocchi_grb_d = bocchi_star_lyr * u.lyr
    bocchi_grb_z = z_at_value(cosmo.luminosity_distance, bocchi_grb_d)
    d = cosmo.lookback_distance(z)
    ld = cosmo.luminosity_distance(z)
    bocchi_grb_ld = cosmo.luminosity_distance(bocchi_grb_z)
    abs_mag = mag - (5 * log10(ld.to(u.parsec).value) - 5)
    bocchi_grb_mag = 5 * log10(bocchi_grb_ld.to(u.parsec).value) - 5 + abs_mag
    return {
        "name": name,
        "z": z,
        "mag": mag,
        "lookback_distance": d.to(u.lyr).value,
        "luminosity_distance": ld.to(u.lyr).value,
        "abs_mag": abs_mag,
        "bocchi_mag": bocchi_grb_mag,
    }

def test_grb_080319B():
    grb = {
        "name": "GRB 080319B",
        "z": 0.9382,
        "mag": 5.8
    }
    print(get_bocchi_grb_mag(grb))


argparser = ArgumentParser(description='Calculate magnitude of Bocchi\'s '\
                           'first star at 0.1Gly')
argparser.add_argument('grb_data', metavar='grb_data.csv',
                       help="data file of GRB (generated with '\
                       '\'simbad-get-grb-data.py\').")
args = argparser.parse_args()

with open(args.grb_data) as f:
    records = csv.DictReader(f)
    results = []
    for r in records:
        z = float(r['Z_VALUE']) if r['Z_VALUE'] else None
        vmag = float(r['FLUX_V']) if r['FLUX_V'] else None
        bmag = float(r['FLUX_B']) if r['FLUX_B'] else None
        if z is None or (vmag is None and bmag is None):
            continue
        main_id = r['MAIN_ID']
        ids = r['IDS'].split('|')
        grb_id = None
        for id in ids:
            if id.startswith('GRB '):
                grb_id = id
                break
        name = grb_id if grb_id else main_id
        mag = vmag if vmag else bmag
        if z < 0:
            print("skip " + name + ": negative red-shift (z=" + str(z) + ")",
                  file=sys.stderr)
            continue
        results.append(get_bocchi_grb_mag({ 'name': name, 'z': z, 'mag': mag }))

cols = [
    'name',
    'z',
    'mag',
    'lookback_distance', 
    'luminosity_distance',
    'abs_mag',
    'bocchi_mag'
]
w = csv.DictWriter(sys.stdout, cols)
w.writeheader()
for r in results:
    w.writerow(r)
