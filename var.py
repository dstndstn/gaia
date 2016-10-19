import numpy as np
import astropy.io.votable
from astrometry.util.fits import fits_table

# Nasty.  From Gaia archive, go to ADQL interface, execute:
#
#  select * from gaiadr1.phot_variable_time_series_gfov;
#  -> download to var.vot
#
#  select * from gaiadr1.gaia_source as s
#  where s.source_id in
#  (select source_id from gaiadr1.phot_variable_time_series_gfov)
#
#  -> download to var-sources.vot


T = astropy.io.votable.parse('var-sources.vot')
# Grab the table
r = T.resources[0]
t = r.tables[0]
# This is a numpy MaskedArray object.
arr = t.array
dt = arr.dtype
# Copy fields into a fits_table object
T = fits_table()
for k,v in dt.fields.items():
    T.set(k, arr[k].filled())

# Convert strings into regular character arrays, not random groups
for col in ['phot_variable_flag',]:
    T.set(col, np.array([t for t in T.get(col)]))

T.writeto('variables-sources.fits')

sourcemap = dict([(sid,i) for i,sid in enumerate(T.source_id)])

V = astropy.io.votable.parse('var.vot')
# Grab the table
r = V.resources[0]
t = r.tables[0]
# This is a numpy MaskedArray object.
arr = t.array
dt = arr.dtype
# Copy fields into a fits_table object
V = fits_table()
for k,v in dt.fields.items():
    V.set(k, arr[k].filled())

V.source_row = np.array([sourcemap[sid] for sid in V.source_id])
    
# Convert strings into regular character arrays, not random groups
#for col in ['type_best_classification',
#            'type2_best_sub_classification',
#            'mode_best_classification']:
#    T.set(col, np.array([t for t in T.get(col)]))
V.writeto('variables.fits')

