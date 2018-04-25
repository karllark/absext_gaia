from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math

import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
import pyvo as vo


def get_gaia_single_star(sname, tap_service):

    # get the coordinates from Simbad
    res = Simbad.query_object(sname)

    # put into a coordinate object to handle the unit conversions
    coordstr = '%s %s' % (res['RA'].data[0], res['DEC'].data[0])
    coord = SkyCoord(coordstr, unit=(u.hourangle, u.degree), frame='icrs')

    query_str = """SELECT * FROM gaiadr2.gaia_source
                WHERE CONTAINS(POINT('ICRS',ra,dec),
                CIRCLE('ICRS',%s,%s,0.00027))=1""" % (coord.ra.degree,
                                                   coord.dec.degree)

    tap_results = tap_service.run_async(query_str)
    #tap_results = tap_service.search(query_str)

    if len(tap_results) > 0:
        return (tap_results[0], len(tap_results))
    else:
        return (None, 0)

if __name__ == "__main__":

    # initialize the GAIA TAP service
    tap_service = vo.dal.TAPService("http://gaia.ari.uni-heidelberg.de/tap")
    print(tap_service)

    # get the list of starnames
    itablename = 'data/mwext_fuse.dat'
    # itablename = 'data/mwext_small.dat'
    snames = Table.read(itablename,
                        format='ascii.fixed_width', guess=False)

    # query GAIA and create a table
    ot = Table(names=('name', 'parallax', 'parallax_error', 'nfound',
                      'G_mag', 'G_flux', 'G_flux_error',
                      'AG', 'AG_lower', 'AG_upper'),
               dtype=('S15', 'f', 'f', 'i', 'f', 'f', 'f', 'f', 'f', 'f'))
    for sname in snames['name']:
        print('trying ', sname)
        # get the GAIA info for one star
        sres, nres = get_gaia_single_star(sname, tap_service)
        if nres == 0:
            print(sname, 'not found')
            # ot.add_row((sname, 0.0, 0.0, 0, 0.0, 0.0, 0.0))
        elif not math.isfinite(sres['parallax']):
            print(sname, 'no parralax')
        else:
            print(sname, sres['parallax'], sres['parallax_error'])
            ot.add_row((sname, sres['parallax'], sres['parallax_error'], nres,
                        sres['phot_g_mean_mag'],
                        sres['phot_g_mean_flux'],
                        sres['phot_g_mean_flux_error'],
                        sres['a_g_val'],
                        sres['a_g_percentile_lower'],
                        sres['a_g_percentile_upper']))

    # output the resulting table
    otablename = itablename.replace('.dat', '_gaia.dat')
    ot.write(otablename, format='ascii.commented_header', overwrite=True)

# possibly useful code for direct GAIA query?
# gaiadr1.gaia_source

# SELECT+*+FROM+fp_psc+WHERE+CONTAINS(POINT('J2000',ra,dec),CIRCLE('J2000',66.76957,26.10453,0.01))=1"


# r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
# r.pprint()
