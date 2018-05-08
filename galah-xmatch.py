import matplotlib
matplotlib.rcParams['figure.figsize'] = (10,8)
import pylab as plt
from astrometry.util.fits import *
from astrometry.util.plotutils import *
import numpy as np
import fitsio
from glob import glob
from astrometry.libkd.spherematch import *
from astrometry.util.starutil_numpy import deg2dist

galah = fits_table('GALAH_DR2_catalog.fits')
print(len(galah), 'GALAH')
galah_kd = tree_build_radec(galah.raj2000, galah.dej2000)

gaia_kd_fn = '/global/cscratch1/sd/dstn/gaia-all.kd.fits'
gaia_kd = tree_open(gaia_kd_fn)
print('Read Gaia KD')

print('Matching...')
rad = 2./3600.
#I,J,d = trees_match(galah_kd, gaia_kd, deg2dist(rad), nearest=True)

### NOTE, this does not seem to be working correctly!!

I,J,D = trees_match(galah_kd, gaia_kd, deg2dist(rad))
print(len(I), 'matches')

# ugh, nearest
bestdist = np.empty(len(galah), np.float32)
bestdist[:] = 1000.
bestmatch = np.empty(len(galah), np.int32)
bestmatch[:] = -1
for i,j,d in zip(I,J,D):
    if d < bestdist[i]:
        bestdist[i] = d
        bestmatch[i] = j

I, = np.nonzero(bestmatch > -1)
J = bestmatch[I]

gaia = fits_table(gaia_kd_fn, rows=J)
Gfiles = fits_table('/global/cscratch1/sd/dstn/gaia-all-filenames.fits')
file_inds = np.unique(gaia.file_index)
TT = []
for f_ind in file_inds:
    ii = np.flatnonzero(gaia.file_index == f_ind)
    rows = gaia.index[ii]
    fn = '/global/project/projectdirs/cosmo/staging/gaia/dr2/gaia-source/' + Gfiles.gaia_filenames[f_ind].strip()
    print('Reading', len(rows), 'from', fn)
    gg = fits_table(fn, rows=rows)
    TT.append(gg)
del gaia
print('Merging...')
gaia = merge_tables(TT)
del TT

print('Re-matching...')
I,J,d = match_radec(galah.raj2000, galah.dej2000, gaia.ra, gaia.dec, rad,
                    nearest=True)
print('Matched', len(I))
gaia.cut(J)

galah.matched_gaia = np.zeros(len(galah), bool)
galah.matched_gaia[I] = True
galah.gaia_match_distance = np.zeros(len(galah), np.float32)
galah.gaia_match_distance[I] = d

for c in gaia.get_columns():
    if c in galah.get_columns():
        c = 'gaia_'+c
    X = gaia.get(c)
    sh = X.shape
    if len(sh) == 1:
        M = np.zeros(len(galah), X.dtype)
        M[I] = X[J]
        galah.set(c, M)
    elif len(sh) == 2:
        M = np.zeros([len(galah)] + sh[1:], X.dtype)
        M[I,:] = X[J,:]
        galah.set(c, M)
    else:
        assert(False)

galah.writeto('/global/cscratch1/sd/dstn/galah-gaia-match.fits')



