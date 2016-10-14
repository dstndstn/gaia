from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from collections import Counter
from glob import glob
import numpy as np
import pylab as plt
from astrometry.util.fits import *
from astrometry.util.util import *
from astrometry.libkd.spherematch import *
from astrometry.util.plotutils import *
from astrometry.util.starutil_numpy import *
import os

def main():

    plt.figure(figsize=(12,6))

    wsubfn = 'wsubs/dr3-sub.fits'
    if not os.path.exists(wsubfn):
        Wsub = decals_sub()
        Wsub.writeto(wsubfn)
    T = fits_table(wsubfn, columns=['ra','dec','decam_flux', 'decam_flux_ivar',
                                    'phot_g_mean_mag', 'type',
                                    'shapeexp_r', 'shapeexp_r_ivar'])

    print('Read', len(T), 'matches')

    T.G = T.phot_g_mean_mag
    T.gflux = T.decam_flux[:,1]
    T.rflux = T.decam_flux[:,2]
    T.zflux = T.decam_flux[:,4]
    T.typ = np.array([t[0] for t in T.type])

    for k,n in Counter(T.typ).most_common():
        print('  type', k, ':', n)

    T.cut((T.decam_flux_ivar[:,1] > 0) *
          (T.decam_flux_ivar[:,2] > 0) *
          (T.decam_flux_ivar[:,4] > 0) *
          (T.gflux > 0) *
          (T.rflux > 0) *
          (T.zflux > 0))
    print(len(T), 'with g,r,z')

    ivg = np.median(T.decam_flux_ivar[:,1])
    ivr = np.median(T.decam_flux_ivar[:,2])
    ivz = np.median(T.decam_flux_ivar[:,4])
    print('Median ivars:', ivg, ivr, ivz)

    T.cut((T.decam_flux_ivar[:,1] > 0.25*ivg) *
          (T.decam_flux_ivar[:,2] > 0.25*ivr) *
          (T.decam_flux_ivar[:,4] > 0.25*ivz))
    print(len(T), 'with g,r,z ivar > 0.25 median')

    T.g = -2.5 * (np.log10(T.gflux) - 9)
    T.r = -2.5 * (np.log10(T.rflux) - 9)
    T.z = -2.5 * (np.log10(T.zflux) - 9)

    plt.clf()
    loghist(T.g - T.r, T.r - T.z, 200, range=((-2,5),(-2,5)))
    plt.xlabel('g - r (mag)')
    plt.ylabel('r - z (mag)')
    plt.title('%i Gaia-DECaLS DR3 matches' % (len(T)))
    plt.savefig('dr3-1.png')

    I = np.flatnonzero(T.typ == 'P')

    plt.clf()
    loghist(T.g[I] - T.r[I], T.r[I] - T.z[I], 200, range=((-2,5),(-2,5)))
    plt.xlabel('g - r (mag)')
    plt.ylabel('r - z (mag)')
    plt.title('%i Gaia-DECaLS DR3 matches: PSFs' % (len(I)))
    plt.savefig('dr3-2.png')

    I = np.flatnonzero(T.typ != 'P')

    plt.clf()
    loghist(T.g[I] - T.r[I], T.r[I] - T.z[I], 200, range=((-2,5),(-2,5)))
    plt.xlabel('g - r (mag)')
    plt.ylabel('r - z (mag)')
    plt.title('%i Gaia-DECaLS DR3 matches: galaxies' % (len(I)))
    plt.savefig('dr3-3.png')

    I = np.flatnonzero(T.typ == 'E')

    plt.clf()
    loghist(T.g[I] - T.r[I], T.r[I] - T.z[I], 200, range=((-2,5),(-2,5)))
    plt.xlabel('g - r (mag)')
    plt.ylabel('r - z (mag)')
    plt.title('%i Gaia-DECaLS DR3 matches: EXP galaxies' % (len(I)))
    plt.savefig('dr3-4.png')

    T[I].writeto('dr3-exp.fits', columns=['ra','dec'])

    f = open('gals1.html', 'w')
    f.write('<html><body>\n')
    for i in I[:50]:
        url = ('http://legacysurvey.org/viewer-dev/jpeg-cutout/?ra=%.5f&dec=%.5f'
               % (T.ra[i], T.dec[i]))
        url2 = ('http://legacysurvey.org/viewer-dev/?ra=%.5f&dec=%.5f'
               % (T.ra[i], T.dec[i]))
        f.write('<a href="%s"><img src="%s"></a>\n' % (url2,url))
    f.write('</body></html>\n\n')

    I = np.argsort(-(T.shapeexp_r * np.sqrt(T.shapeexp_r_ivar) * (T.typ == 'E')))

    f = open('gals2.html', 'w')
    f.write('<html><body>\n')
    for i in I[:50]:
        url = ('http://legacysurvey.org/viewer-dev/jpeg-cutout/?ra=%.5f&dec=%.5f'
               % (T.ra[i], T.dec[i]))
        url2 = ('http://legacysurvey.org/viewer-dev/?ra=%.5f&dec=%.5f'
               % (T.ra[i], T.dec[i]))
        f.write('<a href="%s"><img src="%s"></a>\n' % (url2,url))
    f.write('</body></html>\n\n')

        


    # # Make row-matched table
    # 
    # blanks = fits_table()
    # blanks.matched = np.zeros(len(Tgas) - len(I), bool)
    # 
    # M = Wsub[J]
    # M.matched = np.ones(len(M), bool)
    # M.matchdist = d
    # 
    # MB = merge_tables([M, blanks], columns='fillzero')
    # 
    # rows = np.zeros(len(Tgas), int)
    # rows[:] = -1
    # rows[I] = np.arange(len(I))
    # 
    # MB = MB[rows]
    # MB.writeto('tgas-matched-decals-dr3.fits')



def decals_sub():
    subs = []

    fns = glob('/project/projectdirs/cosmo/data/legacysurvey/dr3/sweep/3.0/sweep-*.fits')
    fns.sort()

    for fn in fns:
        dtag = os.path.basename(fn).replace('sweep-','').replace('.fits','')

        wfn = 'wsubs/decals-dr3-%s.fits' % (dtag)
        if os.path.exists(wfn):
            print('Reading', wfn)
            D = fits_table(wfn)
            subs.append(D)
            continue

        # 
        print('Reading DECaLS catalog', fn)
        D = fits_table(fn)
        print(len(D), 'sources')

        dlo = D.dec.min()
        dhi = D.dec.max()
        rlo = D.ra.min()
        rhi = D.ra.max()
        margin = 0.002

        print('DECaLS Dec range', dlo,dhi, 'RA range', rlo, rhi)

        r,d = (rlo + rhi) / 2., (dlo + dhi) / 2.

        dist = max(degrees_between(r, d, rlo, dlo),
                   degrees_between(r, d, rhi, dhi))
        rad = deg2dist(dist)

        # Finkbeiner's Gaia catalogs
        Nside = 32

        hps = healpix_rangesearch_radec_approx(r, d, rad, Nside)
        print('Healpixes in range:', hps)

        GG = []
        for hpxy in hps:
            hp = healpix_xy_to_ring(hpxy, Nside)

            G = fits_table('/project/projectdirs/cosmo/work/gaia/chunks-gaia_rel1/chunk-%05i.fits' % hp)
            print(len(G), 'Gaia stars in hp', hp)
            I = np.flatnonzero((G.dec > (dlo - margin)) * (G.dec < (dhi + margin)) *
                               (G.ra  > (rlo - margin)) * (G.ra  < (rhi + margin)))
            if len(I) == 0:
                continue
            G.cut(I)
            print(len(G), 'Gaia stars in RA,Dec box')
            GG.append(G)
        G = merge_tables(GG)

        print('Matching...')
        I,J,d = match_radec(G.ra, G.dec, D.ra, D.dec, 1./3600.,
                            nearest=True)
        print(len(I), 'matches')
    
        if len(I) == 0:
            continue

        plt.clf()
        plt.subplot(1,2,1)
        plothist(G.ra, G.dec, doclf=False)
        plt.subplot(1,2,2)
        plothist(G.ra[I], G.dec[I], doclf=False)
        plt.savefig('matched-%s.png' % dtag)

        G.cut(I)
        Dj = D[J]
        Dj.rename('ra',  'ra_decals')
        Dj.rename('dec', 'dec_decals')
        G.add_columns_from(Dj)
        G.match_dist = d
        G.writeto(wfn)
        print('Wrote', wfn)
        subs.append(G)

    subs = merge_tables(subs)
    return subs

if __name__ == '__main__':
    main()
