from astropy.coordinates import SkyCoord
import astropy.units as units
import numpy as np

try:
    from dustmaps.bayestar import BayestarQuery
    from dustmaps.decaps import DECaPSQuery, DECaPSQueryLite
except ImportError:
    BayestarQuery = None
    DECaPSQuery = None
    DECaPSQueryLite = None

def get_ebv(l, b, distance, catalog = 'Bayestar2019', mode = 'samples'):
    """
    Get the E(B-V) value for a given set of coordinates and distance using the Bayestar dust map.

    Parameters:
    l (array): Galactic longitude in degrees.
    b (array): Galactic latitude in degrees.
    distance (array): Distance in parsecs.

    Returns:
    array: E(B-V) values at the specified coordinates and distances.
    """
    if BayestarQuery is None or DECaPSQueryLite is None:
        raise ImportError("dustmaps is required to query dust maps. Install with: pip install dustmaps")

    # Ensure that the Bayestar dust map is downloaded and available
    if catalog == 'Bayestar2019':
        query = BayestarQuery(max_samples=None)
    if catalog == 'decaps':
        query = DECaPSQueryLite()

    # Create a SkyCoord object for the given coordinates
    coords = SkyCoord(l=l * units.degree, b=b * units.degree, distance=distance * units.pc, frame='galactic')

    # Query the Bayestar dust map for the E(B-V) value at the specified coordinates and distance
    ebv = query(coords, mode=mode)

    return ebv
