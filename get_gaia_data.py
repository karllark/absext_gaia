import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import pyvo as vo

res = Simbad.query_object("hd283809")
# print(res.colnames)
# print(res['RA'].data, res['DEC'].data)

coordstr = '%s %s' % (res['RA'].data[0], res['DEC'].data[0])
print(coordstr)
coord = SkyCoord(coordstr, unit=(u.hourangle, u.degree), frame='icrs')
width = u.Quantity(0.1, u.deg)
height = u.Quantity(0.1, u.deg)
print(coord)
print(coord.ra.degree, coord.dec.degree)

tap_service = vo.dal.TAPService("http://gaia.ari.uni-heidelberg.de/tap")
# print(tap_service.tables)

tap_results = tap_service.search("SELECT TOP 10 * FROM  gaiadr1.gaia_source")

query_str = """SELECT * FROM gaiadr1.gaia_source
               WHERE CONTAINS(POINT('ICRS',ra,dec),
               CIRCLE('ICRS',%s,%s,0.01))=1""" % (coord.ra.degree,
                                                  coord.dec.degree)
print(query_str)
tap_results = tap_service.search(query_str)
#                                    CIRCLE('ICRS',70.3530175,25.91344194,0.01))=1

for tr in tap_results:
    print(tr['phot_g_mean_mag'])
    print(tr['parallax'], tr['parallax_error'])
    print(tr['phot_g_mean_flux'], tr['phot_g_mean_flux_error'])
    # print(tr.keys())
# gaiadr1.gaia_source

# SELECT+*+FROM+fp_psc+WHERE+CONTAINS(POINT('J2000',ra,dec),CIRCLE('J2000',66.76957,26.10453,0.01))=1"


# r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
# r.pprint()
