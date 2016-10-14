from __future__ import print_function
import numpy as np
from astrometry.util.fits import *


T = fits_table('data/apass-dr9-2.fits') #, rows=np.arange(10))
#print('columns:', T.get_columns())

print('mobs nobs', T.mobs.min(), T.mobs.max(), T.nobs.min(), T.nobs.max())
print('field recno', T.field.min(), T.field.max(), T.recno.min(), T.recno.max())

T.mobs = T.mobs.astype(np.int16)
T.nobs = T.mobs.astype(np.int16)

assert(np.all(np.unique(T.recno) == (np.arange(len(T))+1)))

assert(len(np.unique(T.recno)) == len(T))

# Reorder by recno
T.cut(np.argsort(T.recno))

assert(np.all(T.recno == np.arange(len(T))+1))

#eq = (T.recno == (np.arange(len(T))+1))
#print(np.sum(eq), 'recno == #row')
#print(np.sum(np.logical_not(eq)), 'recno != #row')

T.rename('raj2000', 'ra')
T.rename('dej2000', 'dec')
T.rename('e_raj2000', 'e_ra')
T.rename('e_dej2000', 'e_dec')
T.rename('b-v', 'bv')
T.rename('e_b-v', 'e_bv')

for band in ['b', "g'", "i'", "r'", 'v']:
    col = 'u_e_%smag' % band
    u = T.get(col)
    T.delete_column(col)
    b2 = band.replace("'", '')
    ui = np.zeros(len(u), np.uint8)
    ui[u == 1] = 1
    ui[np.logical_not(np.isfinite(u))] = 2
    T.set('etype_%smag' % b2, ui)

    if b2 != band:
        T.rename('%smag'   % band, '%smag'   % b2)
        T.rename('e_%smag' % band, 'e_%smag' % b2)

T.about()

cols = ('ra dec e_ra e_dec '
        + 'bmag vmag gmag rmag imag '
        + 'e_bmag e_vmag e_gmag e_rmag e_imag '
        + 'etype_bmag etype_vmag etype_gmag etype_rmag etype_imag '
        + 'field mobs nobs '
        + 'bv e_bv')
T.writeto('data/apass-dr9-3.fits', columns=cols.split())
        
