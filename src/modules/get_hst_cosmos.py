from astropy.io import fits
from astropy import units as u
import requests
import xmltodict


def get_hst_cosmos(source, opt_view=40*u.arcsec):
    """Retrieve HST Cosmos data.

    :param source: source about which to retrieve HST data
    :type source: astropy object? TODO check this!
    :param opt_view: field of view, defaults to 40*u.arcsec
    :type opt_view: astropy unit, optional
    :return: HST data
    :rtype: FITS file
    """
    api_endpoint = "https://irsa.ipac.caltech.edu/cgi-bin/Cutouts/nph-cutouts"

    params = {'locstr': '{} {}'.format(source['ra'], source['dec']),
              'max_size': '180',
              'sizeX': str(int(opt_view.to(u.arcsec).value)),
              'units': opt_view.unit.name,
              'mode': 'PI',
              'mission': 'COSMOS',
              'ntable_cutouts': '1',
              'cutouttbl1': 'acs_mosaic_2.0'}

    try:
        r = requests.get(url=api_endpoint, params=params)
        data = xmltodict.parse(r.content)
    except requests.exceptions.ConnectionError:
        data = None

    try:
        result = fits.open(data['result']['images']['cutouts']['fits'])
    except:
        result = None

    return result