import os
import matplotlib
matplotlib.rcParams['figure.figsize'] = (10,8)
import pylab as plt
from astrometry.libkd.spherematch import *
from astrometry.util.fits import *
import numpy as np
from astrometry.util.starutil_numpy import *
from astrometry.util.plotutils import *
from glob import glob
from collections import Counter
from scipy.interpolate import CubicSpline

S = fits_table('/global/cscratch1/sd/dstn/sdss-specObj-dr14-gaia-match.fits')
S.r_mag = -2.5 * (np.log10(S.calibflux[:,2]) - 9.)
S.star = (S.get('class') == 'STAR  ')
S.qso  = (S.get('class') == 'QSO   ')
S.gal  = (S.get('class') == 'GALAXY')
S.dm = 5.*np.log10(1./(S.parallax/1000.))-5.

if not os.path.exists('GaiaDR2_Passbands.fits'):
    # https://www.cosmos.esa.int/documents/29201/1645651/GaiaDR2_Passbands_ZeroPoints.zip/49cdce41-8eee-655d-7ed2-4e7a83598c1d
    f = open('GaiaDR2_Passbands.dat')
    curves = fits_table()
    curves.wavelength = []
    curves.g = []
    curves.g_err = []
    curves.bp = []
    curves.bp_err = []
    curves.rp = []
    curves.rp_err = []
    
    for line in f.readlines():
        words = line.strip().split()
        curves.wavelength.append(float(words[0]))
        curves.g.append(float(words[1]))
        curves.g_err.append(float(words[2]))
        curves.bp.append(float(words[3]))
        curves.bp_err.append(float(words[4]))
        curves.rp.append(float(words[5]))
        curves.rp_err.append(float(words[6]))
    curves.to_np_arrays()
    curves.g[curves.g > 99] = 0.
    curves.bp[curves.bp > 99] = 0.
    curves.rp[curves.rp > 99] = 0.
    curves.g_err[curves.g_err > 99] = 0.
    curves.bp_err[curves.bp_err > 99] = 0.
    curves.rp_err[curves.rp_err > 99] = 0.
    curves.writeto('GaiaDR2_Passbands.fits')

curves = fits_table('GaiaDR2_Passbands.fits')

fbp = CubicSpline(curves.wavelength, curves.bp)
frp = CubicSpline(curves.wavelength, curves.rp)
fg  = CubicSpline(curves.wavelength, curves.g)

Q = S[S.matched_gaia * S.qso * (S.zwarning==0)]
print(len(Q), 'quasars')

Q.cut(np.argsort(Q.phot_g_mean_mag))
Q.synflux_bp = np.zeros(len(Q))
Q.synflux_rp = np.zeros(len(Q))
Q.synflux_g = np.zeros(len(Q))
for iq in range(len(Q)):
    q = Q[iq]
    fn = '/global/project/projectdirs/cosmo/data/sdss/dr14/sdss/spectro/redux/%s/spectra/lite/%04i/spec-%04i-%i-%04i.fits' % (q.run2d.strip(), q.plate, q.plate, q.mjd, q.fiberid)
    #print(fn)
    print('.', sep='', end='', flush=True)
    spec = fits_table(fn)
    spec.nm = 0.1 * 10.**spec.loglam
    Q.synflux_bp[iq] = np.sum(spec.flux * fbp(spec.nm))
    Q.synflux_rp[iq] = np.sum(spec.flux * frp(spec.nm))
    Q.synflux_g [iq] = np.sum(spec.flux * fg (spec.nm))
    
    if iq and (iq % 1000 == 0):
        fn = '/global/cscratch1/sd/dstn/sdss-qso-synflux-interim.fits'
        Q.writeto(fn)
        print()
        print('Wrote', fn)

fn = '/global/cscratch1/sd/dstn/sdss-qso-synflux.fits'
Q.writeto(fn)
