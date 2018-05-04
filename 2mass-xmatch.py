from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pylab as plt
from astrometry.util.fits import *
from astrometry.util.util import *
from astrometry.libkd.spherematch import *
from astrometry.util.plotutils import *
from astrometry.util.starutil_numpy import *
import os

def main():

    Tgas = fits_table('/project/projectdirs/cosmo/staging/gaia/tgas-source/tgas-source.fits')


    # T = Tgas.copy()
    # # Run the proper motions back to the approximate epoch of 2MASS
    # # (1991) from approx epoch of Gaia DR1 (2015).
    # dyr = 2015. - 1991.
    # 
    # T.ra  -= dyr * Tgas.pmra  / 3600. / np.cos(np.deg2rad(Tgas.dec))
    # T.dec -= dyr * Tgas.pmdec / 3600.

    tmsubfn = '2mass-subs/2mass-sub.fits'
    if not os.path.exists(tmsubfn):
        Tmsub = twomass_sub(T)
        Tmsub.writeto(tmsubfn)
    else:
        Tmsub = fits_table(tmsubfn)


    # Proper motions...

    I,J,d = match_radec(Tgas.ra, Tgas.dec, Tmsub.ra, Tmsub.dec, 4./3600., nearest=True)
    print(len(I), 'matches')
    print(len(np.unique(I)), 'unique TGAS stars matched')

    P    = Tgas.parallax[I]
    Perr = Tgas.parallax_error[I]
    dec  = Tgas.dec[I]
    G    = Tgas.phot_g_mean_mag[I]

    Jmag = Tmsub.j_mag[J]

    plt.clf()
    plt.plot((Tgas.ra[I] - Tmsub.ra[J]) * np.cos(np.deg2rad(Tgas.dec[I])),
             Tgas.dec[I] - Tmsub.dec[J], 'k.')
    plt.xlabel('dRA (deg)')
    plt.ylabel('dDec (deg)')
    plt.savefig('dradec1.png')

    plt.clf()
    plt.plot(Tgas.pmra[I], Tgas.pmdec[I], 'k.')
    plt.xlabel('pmra (mas/yr)')
    plt.ylabel('pmdec (mas/yr)')
    plt.savefig('dradec2.png')

    dyrs = (Tgas.ref_epoch[I] - 2000.) - (Tmsub.jdate[J] - 2451545.0) / 365.25
    print('dyrs: median', np.median(dyrs))

    mra  = dyrs * Tgas.pmra[I] / 1000. / 3600. / np.cos(np.deg2rad(Tgas.dec[I]))
    dra  = Tgas.ra[I] - (Tmsub.ra[J] + mra)
    mdec = dyrs * Tgas.pmdec[I] / 1000. / 3600.
    ddec = Tgas.dec[I] - (Tmsub.dec[J] + mdec)

    plt.clf()
    plt.plot(dra, ddec, 'k.')
    plt.xlabel('dRA (deg)')
    plt.ylabel('dDec (deg)')
    plt.savefig('dradec3.png')


    sys.exit(0)


    K = np.flatnonzero((P > 0) * (P > Perr*10.))
    print(len(K), 'stars with parallax > 0')
    
    plt.clf()
    loghist(G[K] - Jmag[K], G[K] + 5.*np.log10(P[K]), 200)
    plt.xlabel('G - J (mag)')
    plt.ylabel('G + 5 log(parallax)')
    plt.title('%i Gaia-2MASS matches, Parallax SN>10' %
              (len(K)))
    plt.ylim(20, 8)
    plt.xlim(-1, 4)
    plt.savefig('cmd2.png')

    plt.clf()
    loghist(G[K] - Jmag[K], Jmag[K] + 5.*np.log10(P[K]), 200)
    plt.xlabel('G - J (mag)')
    plt.ylabel('J + 5 log(parallax)')
    plt.title('%i Gaia-2MASS matches, Dec %.2f to %.2f, P SN>10' %
              (len(K), dec.min(), dec.max()))
    plt.ylim(16, 6)
    plt.xlim(-1, 4)
    plt.savefig('cmd3.png')


    # cut
    # kk = K[(G[K] - J[K]) < -4]
    # print(len(kk), 'with G-J < -4')
    # html = '<html><body>'
    # for r,d in zip(Tgas.ra[I][kk], Tgas.dec[I][kk]):
    #     html += ('<img src="%s"><img src="%s><img src="%s"><br>' %
    #              ('http://legacysurvey.org/viewer/jpeg-cutout/?ra=%.4f&dec=%.4f&layer=decals-dr3',
    #               'http://legacysurvey.org/viewer/jpeg-cutout/?ra=%.4f&dec=%.4f&layer=sdssco',
    #               'http://legacysurvey.org/viewer/jpeg-cutout/?ra=%.4f&dec=%.4f&layer=unwise-neo1'))
    # html += '</body></html>'
    # f = open('2mass.html', 'w')
    # f.write(html)
    # f.close()


    # Make row-matched table

    blanks = fits_table()
    blanks.matched = np.zeros(len(Tgas) - len(I), bool)

    M = Tmsub[J]
    M.matched = np.ones(len(M), bool)
    M.matchdist = d

    MB = merge_tables([M, blanks], columns='fillzero')

    rows = np.zeros(len(Tgas), int)
    rows[:] = -1
    rows[I] = np.arange(len(I))

    MB = MB[rows]
    MB.writeto('tgas-matched-2mass.fits')




def twomass_sub(T):
    subs = []

    for icat in range(972):
        tmfn = '2mass-subs/2mass-%03i-sub2.fits' % icat
        if os.path.exists(tmfn):
            print('Reading', tmfn)
            W = fits_table(tmfn)
            subs.append(W)
            continue

        print('Reading 2MASS healpix', icat)
        TM = fits_table('/project/projectdirs/cosmo/data/2mass/healpix/2mass_hp%03i.fits' % icat)
        print('Dec range', TM.dec.min(), TM.dec.max())

        dlo = TM.dec.min()
        dhi = TM.dec.max()
        margin = 0.002
        TI = np.flatnonzero((T.dec > (dlo - margin)) * (T.dec < (dhi + margin)))
        print(len(TI), 'stars in Dec range', dlo,dhi)

        print('Matching...')
        I,J,d = match_radec(T.ra[TI], T.dec[TI], TM.ra, TM.dec, 4./3600.)
        print(len(I), 'matches')
    
        if len(I) == 0:
            continue

        subs.append(TM[J])

        TM[J].writeto(tmfn)
        print('Wrote', tmfn)

    subs = merge_tables(subs)
    return subs

if __name__ == '__main__':
    main()
