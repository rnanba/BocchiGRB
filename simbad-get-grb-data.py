#!/usr/bin/env python
import datetime
from astroquery.simbad import Simbad

Simbad.TIMEOUT = 120
Simbad.remove_votable_fields('coordinates')
Simbad.add_votable_fields('ids', 'z_value', 'flux(B)', 'flux(V)')
now = datetime.datetime.now()
table = Simbad.query_criteria(otype='gB')
date_str = now.strftime('%Y%m%d_%H%M%S')
table.write('simbad-grb-data-'+date_str+'.csv', format='csv')
