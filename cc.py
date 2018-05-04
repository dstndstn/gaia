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

if __name__ == '__main__':

    Tgas = fits_table('/project/projectdirs/cosmo/staging/gaia/tgas-source/tgas-source.fits')


    # T2 = fits_table('tycho2.fits', columns=['ra','dec','mag_bt','mag_vt', 'mag_hp',
    #                                         'sigma_mag_bt', 'sigma_mag_vt'])
    # print('Median sigma mag bt:', np.median(T2.sigma_mag_bt))
    # print('Median sigma mag vt:', np.median(T2.sigma_mag_vt))
    # I,J,d = match_radec(Tgas.ra, Tgas.dec, T2.ra, T2.dec, 1./3600., nearest=True)
    # print(len(I), 'matches to Tycho-2')
    # Tgas.cut(I)
    # T2.cut(J)
    # Tgas.add_columns_from(T2)
    # B = Tgas.mag_bt
    # V = Tgas.mag_vt

    afn = 'tgas-matched-apass-dr9.fits'

    T2 = fits_table('data/apass-dr9-3.fits') #, columns=['ra','dec', 'bmag','vmag','gmag','rmag','imag',
    #                                                   'e_bmag','e_vmag','e_gmag','e_rmag','e_imag',
    #                                                   ])
    # for band in 'bvgri':
    #     e = T2.get('e_%smag' % band)
    #     I = np.flatnonzero(np.isfinite(e))
    #     print('Mag', band, '# finite:', len(I), 'median', np.median(e[I]))
    # 
    # # print('Median sigma mag vt:', np.median(T2.sigma_mag_vt))
    I,J,d = match_radec(Tgas.ra, Tgas.dec, T2.ra, T2.dec, 1./3600., nearest=True)
    print(len(I), 'matches to APASS')
    # 
    # # Make row-matched table for APASS
    blanks = fits_table()
    blanks.matched = np.zeros(len(Tgas) - len(J), bool)
    M = T2[J]
    M.matched = np.ones(len(M), bool)
    M.matchdist = d
    MB = merge_tables([M, blanks], columns='fillzero')
    rows = np.zeros(len(Tgas), int)
    rows[:] = -1
    rows[I] = np.arange(len(I))
    MB = MB[rows]
    MB.writeto(afn)
    # 
    # Tgas.cut(I)
    # T2.cut(J)
    # Tgas.add_columns_from(T2)
    # 
    # for band in 'bvgri':
    #     e = T2.get('e_%smag' % band)
    #     I = np.flatnonzero(np.isfinite(e))
    #     print('Matched: mag', band, '# finite:', len(I), 'median', np.median(e[I]))


    TA = fits_table(afn)
    TA.rename('matched', 'matched_a')

    #TM = fits_table('tgas-matched-2mass.fits.gz', columns=['j_mag', 'matched'])
    #TM.writeto('tgas-matched-2mass-cut.fits')
    TM = fits_table('tgas-matched-2mass-cut.fits')
    TM.rename('matched', 'matched_t')
    #W  = fits_table('tgas-matched-wise.fits.gz', columns=['w1mpro', 'matched'])
    #W.writeto('tgas-matched-wise-cut.fits')
    W = fits_table('tgas-matched-wise-cut.fits')
    W.rename('matched', 'matched_w')

    Tgas.add_columns_from(TM)
    Tgas.add_columns_from(W)
    Tgas.add_columns_from(TA)
    print(len(Tgas), 'TGAS sources')

    Tgas.cut(Tgas.matched_w * Tgas.matched_a)
    print(len(Tgas), 'matches from Tgas to WISE and APASS')

    # Tgas.cut(Tgas.matched_t * Tgas.matched_w * Tgas.matched_a)
    # print(len(Tgas), 'matches from Tgas to 2MASS and WISE and APASS')

    #Tgas.cut(Tgas.matched_t * Tgas.matched_w)
    #print(len(Tgas), 'matches from Tgas to 2MASS and WISE')

    B = Tgas.bmag
    V = Tgas.vmag

    G = Tgas.phot_g_mean_mag
    J = Tgas.j_mag
    W1 = Tgas.w1mpro

    P    = Tgas.parallax
    Perr = Tgas.parallax_error

    plt.clf()
    loghist(G - J, J - W1, 200, range=((-1,4),(-1,4)))
    plt.xlabel('G - J (mag)')
    plt.ylabel('J - W1 (mag)')
    plt.title('%i Gaia-AllWISE-2MASS matches' % len(G))
    plt.savefig('cc.png')

    plt.clf()
    loghist(G - J, G - W1, 200, range=((-1,4),(-1,6)))
    plt.xlabel('G - J (mag)')
    plt.ylabel('G - W1 (mag)')
    plt.title('%i Gaia-AllWISE-2MASS matches' % len(G))
    plt.savefig('cc3.png')

    plt.clf()
    I = np.flatnonzero(np.isfinite(B) * np.isfinite(V))
    loghist(B[I] - V[I], G[I] - W1[I], 200, range=((-1,3),(-1,6)))
    plt.xlabel('B - V (mag)')
    plt.ylabel('G - W1 (mag)')
    plt.title('%i Gaia-AllWISE-APASS matches' % len(I))
    plt.savefig('cc4.png')

    K = np.flatnonzero((P > 0) * (P > Perr*30.))
    print(len(K), 'stars with parallax > 0')

    plt.clf()
    loghist(B[K] - V[K], G[K] - W1[K], 200, range=((-1,3),(-1,6)))
    plt.xlabel('B - V (mag)')
    plt.ylabel('G - W1 (mag)')
    plt.title('%i Gaia-AllWISE-APASS matches, parallax S/N > 30' % len(I))
    plt.savefig('cc5.png')

    plt.clf()
    loghist(B[K] - V[K], G[K] + 5.*np.log10(P[K]), 200, range=((-0.5,2.5),(11,20)))
    plt.xlabel('B - V (mag)')
    plt.ylabel('G + 5 log(parallax)')
    plt.title('%i Gaia-APASS matches, parallax SN>30' % len(K))
    plt.ylim(20, 11)
    plt.savefig('cc6.png')

    K = np.flatnonzero((P > 0) * (P > Perr*10.))

    plt.clf()
    loghist(B[K] - V[K], G[K] + 5.*np.log10(P[K]), 200, range=((-0.5,2.5),(8,20)))
    plt.xlabel('B - V (mag)')
    plt.ylabel('G + 5 log(parallax)')
    plt.title('%i Gaia-APASS matches, parallax SN>10' % len(K))
    plt.ylim(20, 8)
    plt.savefig('cc7.png')


    K = np.flatnonzero((P > 0) * (P > Perr*10.))
    print(len(K), 'stars with parallax > 0')

    plt.clf()
    loghist(G[K] - J[K], J[K] - W1[K], 200, range=((-1,4),(-1,4)))
    plt.xlabel('G - J (mag)')
    plt.ylabel('J - W1 (mag)')
    plt.title('%i Gaia-AllWISE-2MASS matches, Parallax SN>10' % len(K))
    plt.savefig('cc2.png')
