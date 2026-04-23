import os
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astroquery.mast import Catalogs
import requests
from astropy.table import Table
from urllib.parse import quote as urlencode
import json
import sys
from astroquery.gaia import Gaia
import numpy as np


c = SkyCoord(l=0.5*u.deg, b=-1.4*u.deg, frame='galactic')
c_right = SkyCoord(l=0.5*u.deg - (23*5*u.arcmin)/2, b=-1.4*u.deg, frame='galactic')
c_left = SkyCoord(l=0.5*u.deg + (23*5*u.arcmin)/2, b=-1.4*u.deg, frame='galactic')
c_up = SkyCoord(l=0.5*u.deg, b=-1.4*u.deg + (45*u.arcmin)/2, frame='galactic')
c_down = SkyCoord(l=0.5*u.deg, b=-1.4*u.deg - (45*u.arcmin)/2, frame='galactic')

# Define Galactic center and box size
l0 = 0.5 * u.deg
b0 = -1.4 * u.deg
half_l = (23 * 5 * u.arcmin) / 2   # same as your left/right extent
half_b = (45 * u.arcmin) / 2       # same as your up/down extent

lmin = (l0 - half_l).to(u.deg)
lmax = (l0 + half_l).to(u.deg)
bmin = (b0 - half_b).to(u.deg)
bmax = (b0 + half_b).to(u.deg)

def smallest_ra_interval_deg(ra_deg):
    """Return smallest RA interval containing all points on a circle."""
    s = np.sort(np.mod(ra_deg, 360.0))
    gaps = np.diff(np.r_[s, s[0] + 360.0])
    i = np.argmax(gaps)  # largest empty gap
    ramin = s[(i + 1) % len(s)]
    ramax = s[i]
    crosses_zero = ramin > ramax
    return ramin, ramax, crosses_zero

# Sample the 4 edges of the Galactic rectangle
N = 400
l_line = np.linspace(lmin.value, lmax.value, N)
b_line = np.linspace(bmin.value, bmax.value, N)

edge_l = np.concatenate([
    l_line,                     # bottom edge (b=bmin)
    l_line,                     # top edge    (b=bmax)
    np.full(N, lmin.value),     # left edge   (l=lmin)
    np.full(N, lmax.value),     # right edge  (l=lmax)
])
edge_b = np.concatenate([
    np.full(N, bmin.value),
    np.full(N, bmax.value),
    b_line,
    b_line,
])

gal_edge = SkyCoord(l=edge_l * u.deg, b=edge_b * u.deg, frame="galactic")
icrs_edge = gal_edge.icrs

ra_vals = icrs_edge.ra.deg
dec_vals = icrs_edge.dec.deg

ramin, ramax, ra_wrap = smallest_ra_interval_deg(ra_vals)
decmin = float(np.min(dec_vals))
decmax = float(np.max(dec_vals))

if not ra_wrap:
    ra_clause = f"ra BETWEEN {ramin:.8f} AND {ramax:.8f}"
else:
    ra_clause = f"(ra >= {ramin:.8f} OR ra <= {ramax:.8f})"

query = f"""
SELECT
    source_id,
    ra, dec,
    phot_rp_mean_mag,
    teff_gspphot,
    mh_gspphot,
    logg_gspphot,
    distance_gspphot
FROM gaiadr3.gaia_source
WHERE
    {ra_clause}
    AND dec BETWEEN {decmin:.8f} AND {decmax:.8f}
    AND phot_rp_mean_mag IS NOT NULL
    AND teff_gspphot IS NOT NULL
    AND mh_gspphot IS NOT NULL
    AND logg_gspphot > 4
    AND distance_gspphot IS NOT NULL
"""

job = Gaia.launch_job_async(query)
tbl = job.get_results()
out_csv = os.path.join(os.getcwd(), "gaia_results.csv")
tbl.write(out_csv, format="csv", overwrite=True)
print(f"Saved CSV: {out_csv}")
# ...existing code...
