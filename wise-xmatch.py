from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pylab as plt
from astrometry.util.fits import *
from astrometry.libkd.spherematch import *
from astrometry.util.plotutils import *
from astrometry.util.starutil_numpy import *
import os

def main():

    wsubfn = 'wsubs/wsub.fits'
    if not os.path.exists(wsubfn):
        Wsub = wise_sub()
        Wsub.writeto(wsubfn)
    else:
        Wsub = fits_table(wsubfn)

    Tgas = fits_table('/project/projectdirs/cosmo/staging/gaia/tgas-source/tgas-source.fits')
    I,J,d = match_radec(Tgas.ra, Tgas.dec, Wsub.ra, Wsub.dec, 4./3600., nearest=True)
    print(len(I), 'matches')
    print(len(np.unique(I)), 'unique TGAS stars matched')
    
    P    = Tgas.parallax[I]
    Perr = Tgas.parallax_error[I]
    dec  = Tgas.dec[I]
    G    = Tgas.phot_g_mean_mag[I]

    W1 = Wsub.w1mpro[J]

    K = np.flatnonzero((P > 0) * (P > Perr*10.))
    print(len(K), 'stars with parallax > 0')
    
    # plt.clf()
    # plt.plot(G[K] - W1[K], G[K] + 5.*np.log10(P[K]), 'b.', alpha=0.2)
    # plt.xlabel('G - W1 (mag)')
    # plt.ylabel('G + 5 log(parallax)')
    # plt.savefig('cmd.png')

    plt.clf()
    loghist(G[K] - W1[K], G[K] + 5.*np.log10(P[K]), 200)
    plt.xlabel('G - W1 (mag)')
    plt.ylabel('G + 5 log(parallax)')
    plt.title('%i Gaia-AllWISE matches, Parallax SN>10' %
              (len(K)))
    plt.ylim(20, 8)
    plt.savefig('cmd2.png')

    plt.clf()
    loghist(G[K] - W1[K], W1[K] + 5.*np.log10(P[K]), 200)
    plt.xlabel('G - W1 (mag)')
    plt.ylabel('W1 + 5 log(parallax)')
    plt.title('%i Gaia-AllWISE matches, Dec %.2f to %.2f, P SN>10' %
              (len(K), dec.min(), dec.max()))
    plt.ylim(16, 4)
    plt.savefig('cmd3.png')


    # Make row-matched table

    blanks = fits_table()
    blanks.matched = np.zeros(len(Tgas) - len(I), bool)

    M = Wsub[J]
    M.matched = np.ones(len(M), bool)
    M.matchdist = d

    MB = merge_tables([M, blanks], columns='fillzero')

    rows = np.zeros(len(Tgas), int)
    rows[:] = -1
    rows[I] = np.arange(len(I))

    MB = MB[rows]
    MB.writeto('tgas-matched-wise.fits')



def wise_sub():
    Wsub = []
    # NW = 0
    # Iall = []
    # Jall = []
    # dall = []
    T = fits_table('~/legacypipe-dir/tycho2.fits.gz')
    #for iwise, (dlo,dhi) in enumerate(wise_dec_range):
    for iwise in range(48):
        wfn = 'wsubs/wsub-%02i.fits' % (iwise+1)
        if os.path.exists(wfn):
            print('Reading', wfn)
            W = fits_table(wfn)
            Wsub.append(W)
            continue

        # kdfn = '/project/projectdirs/cosmo/data/wise/allwise-catalog/wise-allwise-cat-part%02i-radec.kd' % (iwise+1)
        # if not os.path.exists(kdfn):
        #     print('Reading WISE catalog')
        #     Wrd = fits_table('/project/projectdirs/cosmo/data/wise/allwise-catalog/wise-allwise-cat-part%02i-radec.fits' % (iwise+1))
        #     print(len(Wrd), 'WISE sources in Dec range')
        #     print('Building tree...')
        #     kd = tree_build_radec(Wrd.ra, Wrd.dec)
        #     print('Writing tree...')
        #     kdfn = 'wise-allwise-cat-part%02i-radec.kd' % (iwise+1)
        #     tree_save(kd, kdfn)
        #     print('Wrote', kdfn)
        #     tree_free(kd)
        # print('Reading', kdfn)
        # Wkd = tree_open(kdfn)
        # 
        # Tkd = tree_build_radec(T.ra[TI], T.dec[TI])
        # print('Matching...')
        # I,J,d = trees_match(Tkd, Wkd, arcsec2dist(4.))
        # print(len(I), 'matches')
        # tree_close(Wkd)
        # tree_free(Tkd)

        print('Reading WISE catalog')
        Wrd = fits_table('/project/projectdirs/cosmo/data/wise/allwise-catalog/wise-allwise-cat-part%02i-radec.fits' % (iwise+1))
        print(len(Wrd), 'WISE sources in Dec range')
        print('Dec range', Wrd.dec.min(), Wrd.dec.max())

        dlo = Wrd.dec.min()
        dhi = Wrd.dec.max()
        margin = 0.002
        TI = np.flatnonzero((T.dec > (dlo - margin)) * (T.dec < (dhi + margin)))
        print(len(TI), 'Tycho-2 stars in Dec range', dlo,dhi)

        print('Matching...')
        I,J,d = match_radec(T.ra[TI], T.dec[TI], Wrd.ra, Wrd.dec, 4./3600.)
        print(len(I), 'matches')
    
        if len(I) == 0:
            continue

        print('Reading full WISE catalog')
        W = fits_table('/project/projectdirs/cosmo/data/wise/allwise-catalog/wise-allwise-cat-part%02i.fits' % (iwise+1), rows=J)
    
        # ?? Re-match just to avoid remapping J
        print('Re-matching')
        I,J,d = match_radec(T.ra[TI], T.dec[TI], W.ra, W.dec, 4./3600.)
        print(len(I), 'matches,', len(np.unique(I)), 'unique Tycho-2')

        # Record full Tycho-2 indices
        #Iall.append(TI[I])
        #Jall.append(NW + J)
        #dall.append(d)
    
        Wsub.append(W[J])
        #NW += len(J)

        W[J].writeto(wfn)
        print('Wrote', wfn)

        #if iwise == 9:
        #    break
    
    Wsub = merge_tables(Wsub)
    return Wsub


wise_dec_range = [
    (-90.,        -74.4136),
    (-74.413600,  -68.5021),
    (-68.502100,  -63.9184),
    (-63.918400,  -59.9494),
    (-59.949400,  -56.4176),
    (-56.417600,  -53.2041),
    (-53.204100,  -50.1737),
    (-50.173700,  -47.2370),
    (-47.237000,  -44.3990),
    (-44.399000,  -41.6214),
    (-41.621400,  -38.9014),
    (-38.901400,  -36.2475),
    (-36.247500,  -33.6239),
    (-33.623900,  -31.0392),
    (-31.039200,  -28.4833),
    (-28.483300,  -25.9208),
    (-25.920800,  -23.3635),
    (-23.363500,  -20.8182),
    (-20.818200,  -18.2736),
    (-18.273600,  -15.7314),
    (-15.731400,  -13.1899),
    (-13.189900,  -10.6404),
    (-10.640400,  -8.09100),
    (-8.091000 ,  -5.52990),
    (-5.529900 ,  -2.94460),
    (-2.944600 ,  -0.36360),
    (-0.363600 ,  2.225200),
    (2.225200  ,  4.838400),
    (4.838400  ,  7.458000),
    (7.458000  ,  10.09010),
    (10.090100 ,  12.73250),
    (12.732500 ,  15.40770),
    (15.407700 ,  18.10780),
    (18.107800 ,  20.83530),
    (20.835300 ,  23.59230),
    (23.592300 ,  26.38530),
    (26.385300 ,  29.21470),
    (29.214700 ,  32.08660),
    (32.086600 ,  35.00620),
    (35.006200 ,  37.98660),
    (37.986600 ,  41.04520),
    (41.045200 ,  44.18640),
    (44.186400 ,  47.38530),
    (47.385300 ,  50.67800),
    (50.678000 ,  54.08870),
    (54.088700 ,  57.69760),
    (57.697600 ,  61.79110),
    (61.791100 ,  66.58320),
    (66.583200 ,  73.30780),
    (73.307800 ,  90.     ),]

if __name__ == '__main__':
    main()
