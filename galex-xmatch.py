import matplotlib
matplotlib.rcParams['figure.figsize'] = (10,8)
import pylab as plt
from astrometry.util.fits import *
from astrometry.util.plotutils import *
import numpy as np
import fitsio
from glob import glob
from astrometry.libkd.spherematch import *

galex_msc = None
galex_asc = None

def match_fns(X):
    chunki,gaia_fns = X
    gaia = merge_tables([fits_table(fn, columns=['ra','dec','pmra','pmra_error','pmdec','pmdec_error','parallax','parallax_error','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag']) for fn in gaia_fns])

    # MI,MJ,d = match_radec(galex_msc.ra, galex_msc.dec, gaia.ra, gaia.dec, 3./3600., nearest=True)
    # gaia[MJ].writeto('/global/cscratch1/sd/dstn/gaia-galex/gaia-msc-matches-%i.fits' % chunki)

    AI,AJ,d = match_radec(galex_asc.ra, galex_asc.dec, gaia.ra, gaia.dec, 3./3600., nearest=True)
    gaia[AJ].writeto('/global/cscratch1/sd/dstn/gaia-galex/gaia-asc-matches-%i.fits' % chunki)

    print('Chunk', chunki, 'done')

    #return MI,AI

def main1():
    gaia_fns = glob('/global/homes/d/dstn/cosmo/staging/gaia/dr2/fits-flat/Gaia*.fits')
    gaia_fns.sort()
    print(len(gaia_fns), 'Gaia files')
    N = len(gaia_fns)
    args = []
    NC = 100
    while len(gaia_fns):
        args.append((len(args), gaia_fns[:N//NC]))
        gaia_fns = gaia_fns[N//NC:]
    print(len(args), 'chunks')
    
    global galex_msc
    global galex_asc

    # galex_fns = glob('/global/homes/d/dstn/cosmo/data/galex/catalogs/msc/SP_*')
    # galex_msc = merge_tables([fits_table(fn, columns=['ra','dec','e_bv','mag_nuv','magerr_nuv','mag_fuv','magerr_fuv']) for fn in galex_fns])
    # print('GALEX MSC:', len(galex_msc))

    galex_fns = glob('/global/homes/d/dstn/cosmo/data/galex/catalogs/asc/SP_*')
    galex_asc = merge_tables([fits_table(fn, columns=['ra','dec','e_bv','mag_nuv','magerr_nuv','mag_fuv','magerr_fuv']) for fn in galex_fns])
    print('GALEX ASC:', len(galex_asc))

    from astrometry.util.multiproc import multiproc

    mp = multiproc(16)
    mp.map(match_fns, args)

    #galex_kd = tree_build_radec(galex)
    #MI,MJ,d = match_radec(galex.ra, galex.dec, gaia.ra, gaia.dec, 3./3600., nearest=True)
    #X = spherematch_c.nearest(kd1, kd2, maxradius, notself)

def main2():
    #main1()

    if True:
        fns = glob('/global/cscratch1/sd/dstn/gaia-galex/gaia-asc-matches-*.fits')
        fns.sort()
        gaia_asc = merge_tables([fits_table(fn) for fn in fns])
        print(len(gaia_asc), 'Gaia/ASC')
    
        galex_fns = glob('/global/homes/d/dstn/cosmo/data/galex/catalogs/asc/SP_*')
        galex_asc = merge_tables([fits_table(fn, columns=['ra','dec','e_bv','mag_nuv','magerr_nuv','mag_fuv','magerr_fuv']) for fn in galex_fns])
        print('GALEX ASC:', len(galex_asc))
    
        AI,AJ,d = match_radec(galex_asc.ra, galex_asc.dec, gaia_asc.ra, gaia_asc.dec, 3./3600., nearest=True)
        mgalex = galex_asc[AI]
        mgaia = gaia_asc[AJ]
        mgalex.rename('ra','ra_galex')
        mgalex.rename('dec','dec_galex')
        mgaia.add_columns_from(mgalex)
        mgaia.writeto('/global/cscratch1/sd/dstn/gaia-galex/galex-gaia-asc.fits')

        unmatched = np.ones(len(galex_asc), bool)
        unmatched[AI] = False
        galex_asc[unmatched].writeto('/global/cscratch1/sd/dstn/gaia-galex/galex-asc-unmatched.fits')

    import sys
    sys.exit(0)

    fns = glob('/global/cscratch1/sd/dstn/gaia-galex/gaia-msc-matches-*.fits')
    fns.sort()
    gaia_msc = merge_tables([fits_table(fn) for fn in fns])
    print(len(gaia_msc), 'Gaia/MSC')

    galex_fns = glob('/global/homes/d/dstn/cosmo/data/galex/catalogs/msc/SP_*')
    galex_msc = merge_tables([fits_table(fn, columns=['ra','dec','e_bv','mag_nuv','magerr_nuv','mag_fuv','magerr_fuv']) for fn in galex_fns])
    print('GALEX MSC:', len(galex_msc))

    AI,AJ,d = match_radec(galex_msc.ra, galex_msc.dec, gaia_msc.ra, gaia_msc.dec, 3./3600., nearest=True)
    mgalex = galex_msc[AI]
    mgaia = gaia_msc[AJ]
    mgalex.rename('ra','ra_galex')
    mgalex.rename('dec','dec_galex')
    mgaia.add_columns_from(mgalex)
    mgaia.writeto('/global/cscratch1/sd/dstn/gaia-galex/galex-gaia-msc.fits')



if __name__ == '__main__':
    galex_msc_fns = glob('/global/project/projectdirs/cosmo/data/galex/catalogs/msc/SP_*')
    galex_msc = merge_tables([fits_table(fn) for fn in galex_msc_fns])
    print('GALEX MSC:', len(galex_msc))

    galex_asc_fns = glob('/global/project/projectdirs/cosmo/data/galex/catalogs/asc/SP_*')
    galex_asc = merge_tables([fits_table(fn) for fn in galex_asc_fns])
    print('GALEX ASC:', len(galex_asc))

    G = fits_table('/global/cscratch1/sd/dstn/gaia-all.fits',
                   columns=['ra','dec'])
    print('Read', len(G), 'Gaia stars')
    G.rows = np.arange(len(G))

    allrows = []
    #allrows_msc = []
    decs = np.linspace(-90, 90, 19)
    for dlo,dhi in zip(decs, decs[1:]):
        K = np.flatnonzero((G.dec >= dlo) * (G.dec <= dhi))
        print('Matching Dec range', dlo,dhi, '->', len(K), 'Gaia')
        I,J,d = match_radec(galex_asc.ra, galex_asc.dec, G.ra[K], G.dec[K], 3./3600., nearest=True)
        print(len(I), 'ASC matches')
        rows = G.rows[K][J]
        allrows.append(rows)
        I,J,d = match_radec(galex_msc.ra, galex_msc.dec, G.ra[K], G.dec[K], 3./3600., nearest=True)
        print(len(I), 'MSC matches')
        rows = G.rows[K][J]
        allrows.append(rows)
    del G
    rows = np.hstack(allrows)
    rows = np.unique(rows)
    del allrows
    print('Reading', len(rows))
    G = fits_table('/global/cscratch1/sd/dstn/gaia-all.fits', rows=rows, columns=['file_index', 'index'])
    print('Read', len(G), 'Gaia references')

    Gfiles = fits_table('/global/cscratch1/sd/dstn/gaia-all-filenames.fits')
    file_inds = np.unique(G.file_index)
    TT = []
    for f_ind in file_inds:
        ii = np.flatnonzero(G.file_index == f_ind)
        rows = G.index[ii]
        fn = '/global/project/projectdirs/cosmo/staging/gaia/dr2/fits-flat/' + Gfiles.gaia_filenames[f_ind].strip()
        print('Reading', len(rows), 'from', fn)
        gg = fits_table(fn, rows=rows)
        TT.append(gg)
    del G
    print('Merging...')
    TT = merge_tables(TT)

    print('Re-matching...')
    I,J,d = match_radec(galex_asc.ra, galex_asc.dec, TT.ra, TT.dec, 3./3600., nearest=True)
    print('Matched', len(I), 'ASC')

    GM = galex_asc[I]
    TM = TT[J]
    for c in TM.get_columns():
        C = TM.get(c)
        if c in GM.get_columns():
            c = 'gaia_' + c
        GM.set(c, C)
    GM.writeto('/global/cscratch1/sd/dstn/galex-asc-gaia-match.fits')

    print('Re-matching...')
    I,J,d = match_radec(galex_msc.ra, galex_msc.dec, TT.ra, TT.dec, 3./3600., nearest=True)
    print('Matched', len(I), 'MSC')
    GM = galex_msc[I]
    TM = TT[J]
    for c in TM.get_columns():
        C = TM.get(c)
        if c in GM.get_columns():
            c = 'gaia_' + c
        GM.set(c, C)
    GM.writeto('/global/cscratch1/sd/dstn/galex-msc-gaia-match.fits')
