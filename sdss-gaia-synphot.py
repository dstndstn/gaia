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

print('Reading spec table...')
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

S.synflux_bp     = np.zeros(len(S), np.float32)
S.synflux_bp_err = np.zeros(len(S), np.float32)
S.synflux_rp     = np.zeros(len(S), np.float32)
S.synflux_rp_err = np.zeros(len(S), np.float32)
S.synflux_g      = np.zeros(len(S), np.float32)
S.synflux_g_err  = np.zeros(len(S), np.float32)

I = np.flatnonzero((S.zwarning == 0) * (S.matched_gaia))
print('Synphot for', len(I), 'spectra...')
for nn,ii in enumerate(I):
    q = S[ii]
    fn = '/global/project/projectdirs/cosmo/data/sdss/dr14/sdss/spectro/redux/%s/spectra/lite/%04i/spec-%04i-%i-%04i.fits' % (q.run2d.strip(), q.plate, q.plate, q.mjd, q.fiberid)
    #print(fn)
    #print('.', sep='', end='', flush=True)
    spec = fits_table(fn)
    nm = 0.1 * 10.**spec.loglam

    hplanck = 6.626e-34 # J/s
    c = 3.00e8 # m/s
    Agaia = 0.7278 # m^2  -- Area of Gaia's primary mirror(s)
    # width of frequency bins in spectrum
    dnm = np.zeros_like(nm)
    dnm[0] = nm[1]-nm[0]
    dnm[-1] = nm[-1]-nm[-2]
    dnm[1:-1] = (nm[2:] - nm[:-2])/2.
    flux_err = np.sqrt(1./ spec.ivar)
    def fluxdensity(flux, nm):
        return flux * 1e-19 * 1e-9*nm / (hplanck * c)
    fluxdens = fluxdensity(spec.flux, nm)  # in photons / (s nm m^2)
    fluxdens_err = fluxdensity(flux_err, nm)
    f_b = fbp(nm)
    f_r = frp(nm)
    f_g = fg (nm)
    S.synflux_bp    [ii] = np.sum(fluxdens     * f_b * dnm) * Agaia
    S.synflux_bp_err[ii] = np.sum(fluxdens_err * f_b * dnm) * Agaia
    S.synflux_rp    [ii] = np.sum(fluxdens     * f_r * dnm) * Agaia
    S.synflux_rp_err[ii] = np.sum(fluxdens_err * f_r * dnm) * Agaia
    S.synflux_g     [ii] = np.sum(fluxdens     * f_g * dnm) * Agaia
    S.synflux_g_err [ii] = np.sum(fluxdens_err * f_g * dnm) * Agaia

    print(nn, '%.1f %.1f (%5.3f) %.1f %.1f (%5.3f) %.1f %.1f (%5.3f) ' %
          (S.synflux_bp[ii], S.phot_bp_mean_flux[ii],
           S.synflux_bp[ii] /S.phot_bp_mean_flux[ii],
           S.synflux_rp[ii], S.phot_rp_mean_flux[ii],
           S.synflux_rp[ii] /S.phot_rp_mean_flux[ii],
           S.synflux_g [ii], S.phot_g_mean_flux [ii],
           S.synflux_g [ii] /S.phot_g_mean_flux [ii]))

    if nn and (nn % 100000 == 0):
        fn = '/global/cscratch1/sd/dstn/sdss-gaia-synflux-interim2.fits'
        S[:ii].writeto(fn)
        print()
        print('Wrote', fn)

fn = '/global/cscratch1/sd/dstn/sdss-gaia-synflux.fits'
S.writeto(fn)

