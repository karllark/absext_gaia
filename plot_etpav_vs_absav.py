
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os.path
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from astropy.table import Table
from astropy.io import fits

from measure_extinction.stardata import StarData

def compute_absext(redplx, compplx,
                   redmag, compmag):
    """
    Compute the absolute extinction given the parallaxes and magnitudes
    for reddened and comparision stars.  Done assuming both have the
    same intrinsic luminosity.
    """
    if (redplx[0] > 0) & (compplx[0] > 0):
        dmagr = math.log10(redplx[0])
        dmagc = math.log10(compplx[0])
        gext = (redmag[0] - compmag[0] + 5.*dmagr - 5.*dmagc)
        gext_unc = (redmag[1]**2 + compmag[1]**2
                    + 5.*redplx[1]/(redplx[0]*math.log(10.))
                    + 5.*compplx[1]/(compplx[0]*math.log(10.)))
        return (gext, gext_unc)
    else:
        return (None, None)

if __name__ == "__main__":

    # filenames
    # FUSE Sample
    gtablename = 'data/mwext_fuse_gaia.dat'
    ptablename = 'data/mwext_fuse_pairs.dat'
    fitspath = '/home/kgordon/Dust/Ext/FUSE/ext1_fuse'
    # IUE Sample
    # gtablename = 'data/mwext_resolv_gaia.dat'
    # ptablename = 'data/mwext_fuse_pairs.dat'
    # fitspath = '/home/kgordon/Dust/Ext/ext0_fuse'

    # gtablename = 'data/mwext_resolv_gaia.dat'

    # get the gaia data
    gt = Table.read(gtablename, format='ascii.commented_header')

    # get the extinction pair names
    pt = Table.read(ptablename, format='ascii.fixed_width')

    # compute the absolute extinction and tabulate the extrapolated value
    gext = []
    gext_unc = []
    etp_av = []
    etp_av_unc = []
    abs_av = []
    abs_av_unc = []
    for rname, cname in zip(pt['rname'].data, pt['cname'].data):
        ri, = np.where(gt['name'].data == rname)
        ci, = np.where(gt['name'].data == cname)
        if len(ri)*len(ci) > 0:
            pr = gt['parallax'][ri[0]]
            pr_unc = gt['parallax_error'][ri[0]]
            pc = gt['parallax'][ci[0]]
            pc_unc = gt['parallax_error'][ci[0]]
            gmagr = gt['G_mag'][ri[0]]
            gmagr_unc = gmagr*gt['G_flux_error'][ri[0]]/gt['G_flux'][ri[0]]
            gmagc = gt['G_mag'][ci[0]]
            gmagc_unc = gmagc*gt['G_flux_error'][ci[0]]/gt['G_flux'][ci[0]]

            redplx = (gt['parallax'][ri[0]], gt['parallax_error'][ri[0]])
            redmag = (gt['G_mag'][ri[0]],
                      gt['G_mag'][ri[0]]*gt['G_flux_error'][ri[0]]
                      / gt['G_flux'][ri[0]])
            compplx = (gt['parallax'][ci[0]], gt['parallax_error'][ci[0]])
            compmag = (gt['G_mag'][ci[0]],
                       gt['G_mag'][ci[0]]*gt['G_flux_error'][ci[0]]
                       / gt['G_flux'][ci[0]])
            absext = compute_absext(redplx, compplx, redmag, compmag)
            if absext[0] is not None:
                gext.append(absext[0])
                gext_unc.append(absext[1])

                # calculate the V band abs extinction
                redname = '/home/kgordon/Python_git/measured_extcurves/data/DAT_files/%s.dat' % (rname)
                compname = '/home/kgordon/Python_git/measured_extcurves/data/DAT_files/%s.dat' % (cname)
                if os.path.isfile(redname) and os.path.isfile(compname):
                    redstar = StarData(redname, photonly=True)
                    compstar = StarData(compname, photonly=True)
                    v_absext = compute_absext(redplx, compplx,
                                              redstar.data['BAND'].bands['V'],
                                              compstar.data['BAND'].bands['V'])
                    abs_av.append(v_absext[0])
                    abs_av_unc.append(v_absext[1])
                else:
                    abs_av.append(-1.0)
                    abs_av_unc.append(-1.0)

                # get the extrapolated derived value from FITS extinction curve
                fname = "%s/%s_%s/%s_%s_ext_bin.fits" \
                    % (fitspath, rname.lower(), cname.lower(),
                       rname.lower(), cname.lower())
                if os.path.isfile(fname):
                    hdu = fits.open(fname)
                    etp_av.append(hdu[0].header['AV'])
                    etp_av_unc.append(hdu[0].header['AV_UNC'])
                else:
                    etp_av.append(-1.0)
                    etp_av_unc.append(-1.0)

    # plot the two distances versus each other
    fontsize = 18

    font = {'size': fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=1)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('xtick.minor', width=2)
    matplotlib.rc('ytick.major', width=2)
    matplotlib.rc('ytick.minor', width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

    # plot 1 to 1 line
    # ax.plot([1.0, 1e5], [1.0, 1e5], 'k-')

    # ax.errorbar(etp_av, gext, xerr=etp_av_unc, yerr=gext_unc,
    #            fmt='bo', capsize=3)
    # ax.errorbar(etp_av, abs_av, xerr=etp_av_unc, yerr=abs_av_unc,
    #            fmt='go', capsize=3)

    ax.plot(gext, gext_unc, 'bo', label='G')
    ax.plot(abs_av, abs_av_unc, 'go', label='V')
    ax.legend()

    ax.set_xlabel(u'$A(\lambda)$')
    ax.set_ylabel(u'$\sigma[A(\lambda)]$')
#    ax.set_xlabel(u'$A(V)_{etp}$')
#    ax.set_ylabel(u'$A(G)_{abs}$')
#    ax.set_xscale('log')
#    ax.set_yscale('log')
    ax.set_xlim(-2., 4.)
    ax.set_ylim(0., 1.)

    plt.show()
