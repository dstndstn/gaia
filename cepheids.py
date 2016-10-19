import astropy.io.votable
from astrometry.util.fits import fits_table

# Nasty.  From Gaia archive, go to ADQL interface, execute:
#   select * from gaiadr1.cepheid
# and download as VOTable. -> cepheids.vot

T = astropy.io.votable.parse('cepheids.vot')
# Grab the table
r = T.resources[0]
t = r.tables[0]
# This is a numpy MaskedArray object.
arr = t.array
# Copy fields into a fits_table object
T = fits_table()
for k,v in dt.fields.items():
    T.set(k, arr[k].filled())

# Convert strings into regular character arrays, not random groups
for col in ['type_best_classification',
            'type2_best_sub_classification',
            'mode_best_classification']:
    T.set(col, np.array([t for t in T.get(col)]))
T.writeto('cepheids.fits')
