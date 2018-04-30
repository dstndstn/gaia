import pylab as plt
from astrometry.libkd.spherematch import *
from astrometry.util.fits import *
import numpy as np
from astrometry.util.starutil_numpy import *
from astrometry.util.plotutils import *
from glob import glob
from astrometry.util.util import *
import fitsio
from astrometry.util.file import *

ra, dec, radius = 290.5, 44., 11.5

#kd = tree_open('/global/cscratch1/sd/dstn/gaia-all.kd.fits')
# xyz = radectoxyz(290.5, 44.)
# radius = deg2dist(11.5)
# rows = kd.search(xyz, radius, 0, 0)
# print(len(rows), 'in field')
# 
# T = fits_table('/global/cscratch1/sd/dstn/gaia-all.kd.fits', rows=rows)
# print('Matching...')
# I = match_radec(T.ra, T.dec, T.ra, T.dec, 32./3600., notself=True, indexlist=True)
# print('Found', len(I), 'matches')
# 
# Igood, = np.nonzero(I)
# print(len(Igood), 'stars with matches')
# J = np.hstack((I[igood] for igood in Igood))
#Trows = np.unique(np.append(Igood,J))
#Tgood = T[Trows]
#files = np.unique(Tgood.file_index)
#F = fits_table('/global/cscratch1/sd/dstn/gaia-all-filenames.fits')
# TT = []
# for fnum in files:
#     fn = filenames[fnum]
#     print('Reading', fn)
#     ii = np.flatnonzero(Tgood.file_index == fnum)
#     print(len(ii), 'rows')
#     #TT.append(fits_table(fn, rows=Tgood.index[ii]))
#     TT.append(fits_table(fn))
# TT = merge_tables(TT)

F = fits_table('/global/cscratch1/sd/dstn/gaia-all-filenames.fits')
filenames = ['/global/project/projectdirs/cosmo/staging/gaia/dr2/fits-flat/' + fn.strip() for fn in F.gaia_filenames]

print('Finding healpixes in range..')
#hps = healpix_rangesearch_radec(290.5, 44., 11.5, 1024)
hps = healpix_rangesearch_radec_approx(ra, dec, radius, 1024)
print('Found', len(hps), 'in range')
# Convert to nested healpix index
hpnest = [healpix_xy_to_nested(hp, 1024) for hp in hps]

#HEALpix level 10 = source_id / 549755813888
# sourceid_ranges = []
# for hp in hpnest:
#     sourceid_min = hp * 2**39
#     sourceid_max = (hp+1) * 2**39
#     sourceid_ranges.append((sourceid_min, sourceid_max))
file_ranges = []
for fn in filenames:
    words = os.path.basename(fn.replace('.fits','')).split('_')
    lo = int(words[1])
    hi = int(words[2])
    file_ranges.append((fn,lo,hi))

TT = []
hpset = set(hpnest)
keepmin = min(hpnest)
keepmax = max(hpnest)
for i,(fn,lo,hi) in enumerate(file_ranges):
    hplo = lo//(2**39)
    hphi = hi//(2**39)
    print(i, 'hp range', hplo, hphi, 'vs', keepmin, 'to', keepmax)
    if hplo < keepmin or hphi > keepmax:
        continue
    keep = False
    for hp in range(hplo, hphi+1):
        if hp in hpset:
            keep = True
            break
    if not keep:
        continue
    print('Reading', i, fn)
    TT.append(fits_table(fn))
TT = merge_tables(TT)
print('Read total of', len(TT), 'sources')

I,J,d = match_radec(TT.ra, TT.dec, ra, dec, radius)
TT.cut(I)
print('Cut to', len(TT), 'really in range')

print('Matching...')
I = match_radec(TT.ra, TT.dec, TT.ra, TT.dec, 32./3600.,
                notself=True, indexlist=True)
Igood, = np.nonzero(I)
Jgood = np.hstack([I[i] for i in Igood])
print(len(Igood), 'stars have matches')
print('Total of', len(Jgood), 'matches')
Ikeep = np.unique(np.hstack((Igood, Jgood)))
print('Unique stars with matches:', len(Ikeep))
maxmatches = max([len(I[i]) for i in Igood])
print('Max number of matches:', maxmatches)

matchinds = np.empty((len(I), maxmatches), np.int32)
matchinds[:,:] = -1
for j,i in enumerate(Igood):
    inds = I[i]
    matchinds[i,:len(inds)] = inds

fitsio.write('/global/cscratch1/sd/dstn/kepler-matches-2.fits', matchinds)
pickle_to_file(I, '/global/cscratch1/sd/dstn/kepler-matches-2.pickle')
TT.writeto('/global/cscratch1/sd/dstn/gaia-kepler-2.fits')
