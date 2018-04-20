from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from astropy.table import Table


if __name__ == "__main__":

    # get the gaia data
    gtablename = 'data/mwext_resolv_gaia.dat'
    gt = Table.read(gtablename, format='ascii.commented_header')

    # get the spectroscopic distances
    st = Table.read('data/valencic04_datafile4.txt', format='cds')

    # list of standards
    iue_standards_i = ['hd188209', 'hd167756', 'hd204172', 'hd150898',
                       'hd064760', 'hd150168', 'hd091316', 'hd040111',
                       'hd165024']
    iue_standards_iii = ['hd063922', 'hd119159', 'hd046328', 'hd062747',
                         'hd051283', 'hd063465', 'hd079447', 'hd195986']
    iue_standards_v = ['hd047839', 'hd214680', 'hd036512', 'hd031726',
                       'hd064802', 'hd032630', 'hd065904', 'hd034759']
    standards = iue_standards_i + iue_standards_iii + iue_standards_v

    # get the two distances in matched vectors
    distname = []
    gdist = []
    gdist_min = []
    gdist_max = []
    sdist = []
    sdist_min = []
    sdist_max = []
    for k, cname in enumerate(gt['name']):
        indxs, = np.where(cname.upper() == st['Name'])
        if len(indxs) != 1:
            print('error: ', cname, len(indxs))
        else:
            m = indxs[0]
            distname.append(cname)
            sdist.append(st['Distance'][m])
            sdist_min.append(st['MinDist'][m])
            sdist_max.append(st['MaxDist'][m])
            gdist.append(1.e3/gt['parallax'][k])
            gdist_min.append(1e3/(gt['parallax'][k] + gt['parallax_error'][k]))
            gdist_max.append(1e3/(gt['parallax'][k] - gt['parallax_error'][k]))

    gdist = np.array(gdist)
    gdist_min = np.array(gdist_min)
    gdist_max = np.array(gdist_max)
    sdist = np.array(sdist)
    sdist_min = np.array(sdist_min)
    sdist_max = np.array(sdist_max)

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
    ax.plot([1.0, 1e5], [1.0, 1e5], 'k-')

    xerror_asym = [gdist-gdist_min, gdist_max-gdist]
    yerror_asym = [sdist-sdist_min, sdist_max-sdist]
    ax.errorbar(gdist, sdist, xerr= xerror_asym, yerr=yerror_asym,
                fmt='bo', capsize=3)

    # paint the standard stars a different color
    for cname in standards:
        indxs, = np.where(cname == distname)
        if len(indxs) == 1:
            k = indxs[0]
            ax.errorbar([gdist[k]], [sdist[k]], 'ro')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e2, 2e4)
    ax.set_ylim(1e2, 2e4)

    plt.show()
