from __future__ import print_function
import numpy as np
import pylab as plt
from astrometry.util.fits import fits_table
from collections import Counter

V = fits_table('variables.fits')

S = np.unique(V.source_id)
print(len(S), 'source ids')

G = fits_table('variables-sources.fits')
G.nepochs = np.zeros(len(G), int)

n = Counter(V.source_row)
for i,nn in n.most_common():
    G.nepochs[i] = nn

plt.clf()
# plt.hist(G.phot_g_mean_flux, 100, range=(0, 20000), log=True)
plt.scatter(G.ra, G.dec, c=G.nepochs.astype(float), edgecolors='none')
plt.colorbar()
plt.title('Gaia variable sources: Nepochs')
plt.xlabel('RA (deg)')
plt.ylabel('Dec (deg)')
plt.savefig('g.png')

nplotted = 0
plt.clf()
for s in S:
    I = np.flatnonzero(V.source_id == s)
    print(len(I), 'observations of', s)
    mean_g = np.mean(V.g_flux[I])
    if mean_g < 1000:
        print('Mean flux', mean_g, '-- skipping')
        continue
    # t = V.observation_time[I]
    # dt = np.diff(t)
    # print('dt', dt)
    # I = I[:-1][dt <= 1]
    plt.plot(V.observation_time[I], V.g_flux[I], '.-')
    nplotted += 1
    if nplotted == 25:
        break

# plt.xlim(1660, 1690)
#plt.xlim(1666, 1673)
plt.xlim(1665, 1685)
    
plt.xlabel('Obs time')
plt.ylabel('g flux')
plt.yscale('symlog')
plt.title('Gaia variable objects')
plt.savefig('v.png')
