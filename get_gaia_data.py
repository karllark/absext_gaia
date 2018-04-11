import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad

res = Simbad.query_object("hd283809")
# print(res.colnames)
# print(res['RA'].data, res['DEC'].data)

coordstr = '%s %s' % (res['RA'].data[0], res['DEC'].data[0])
print(coordstr)
coord = SkyCoord(coordstr, unit=(u.hourangle, u.degree), frame='icrs')
width = u.Quantity(0.1, u.deg)
height = u.Quantity(0.1, u.deg)
print(coord)

r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
r.pprint()
