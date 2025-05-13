import requests
from PIL import Image, ImageEnhance
from io import BytesIO
from urllib.error import HTTPError, URLError

from astropy.io import ascii, fits
from astropy import units as u
from astropy.wcs import WCS
from astroquery.skyview import SkyView
import numpy as np

from src.modules.panstarrs_fcns import getcolorim, geturl
from src.modules.logger import Logger

logger = Logger.get_logger()


def get_skyview(hi_pos, opt_view=6*u.arcmin, survey='DSS2 Blue', cache=True):
    """Retrieve the optical image from a certain pointing.

    :param hi_pos: position in the radio spectral line data
    :type hi_pos: astropy SkyCoord object
    :param opt_view: size of the optical image, defaults to 6*u.arcmin
    :type opt_view: astropy.units.arcmin, optional
    :param survey: survey containing optical data, defaults to 'DSS2 Blue'
    :type survey: str, optional
    :param cache: are the downloaded ancillary files cached
    :type cache: bool
    :return: optical image
    :rtype: astropy HDUList
    """
    # DSS2 Blue images have a 1 arc/pix pixel scale, but retrieving ~the pixel scale doesn't work.
    opt_pixels = (opt_view.to(u.arcsec).value * 2).astype(int)

    try:
        # Get a survey image from SkyView:
        if hi_pos.frame.name == 'galactic':
            path = SkyView.get_images(position=hi_pos.to_string('hmsdms'), coordinates='galactic', width=opt_view[0],
                                        height=opt_view[-1], survey=[survey], pixels=opt_pixels, cache=cache)
        elif (not hi_pos.equinox) or (hi_pos.frame.name == 'icrs'):
            path = SkyView.get_images(position=hi_pos.to_string('hmsdms'), coordinates='ICRS', width=opt_view[0],
                                        height=opt_view[-1], survey=[survey], pixels=opt_pixels, cache=cache)
        # Note that there seems to be a bug in SkyView that it sometimes won't retrieve non-J2000.0.  Keep an eye on this!
        else:
            path = SkyView.get_images(position=hi_pos.to_string('hmsdms'), coordinates=hi_pos.equinox.value,
                                        width=opt_view[0], height=opt_view[-1], survey=[survey], pixels=opt_pixels,
                                        cache=cache)
    except requests.exceptions.ConnectionError:
        path = []
        
    if len(path) != 0:
        logger.info("\tSurvey image retrieved from {}.".format(survey))
        result = path[0]
    else:
        logger.warning("\tNo {} image retrieved.  Bug, or server error?  Try again later?".format(survey))
        result = None

    return result


def get_panstarrs(hi_pos, opt_view=6*u.arcmin):
    """Get PanSTARRS false color image and r-band fits (for the WCS).

    :param hi_pos: position in the radio spectral line data
    :type hi_pos: astropy SkyCoord object
    :param opt_view: size of the optical image, defaults to 6*u.arcmin
    :type opt_view: astropy.units.arcmin, optional
    :return: color image and FITS header
    :rtype: Tuple[color_im, FITS_header] TODO check exact types!
    """
    #  Get PanSTARRS false color image and r-band fits (for the WCS).
    pstar_pixsc = 0.25
    if len(opt_view) > 1:
        logger.warning("\tPanSTARRS only returns square images; taking largest dimension.")
        opt_view = np.max(opt_view)
    try:
        path = geturl(hi_pos.ra.deg, hi_pos.dec.deg, size=int(opt_view.to(u.arcsec).value / pstar_pixsc),
                    filters="r", format="fits")
    except URLError:
        path = []

    if len(path) != 0:
        fits_head = fits.getheader(path[0])
        color_im = getcolorim(hi_pos.ra.deg, hi_pos.dec.deg, size=int(opt_view.to(u.arcsec).value / pstar_pixsc),
                              filters="gri")
        logger.info("\tOptical false color image retrieved from PanSTARRS.")
    else:
        logger.warning("\tNo PanSTARRS false color image retrieved.  Server or connection error," \
              " or no PanSTARRS coverage?")
        fits_head = None
        color_im = None

    return color_im, fits_head


def get_decals(hi_pos, opt_view=6*u.arcmin, decals='decals'):
    """Get DECaLS false color image and g-band fits (for the WCS).

    :param hi_pos: position in the radio spectral line data
    :type hi_pos: astropy SkyCoord object
    :param opt_view: size of the optical image, defaults to 6*u.arcmin
    :type opt_view: astropy.units.arcmin, optional
    :return: color image and FITS header
    :rtype: Tuple[color_im, FITS_header] TODO check exact types!
    """
    # Get DECaLS false color image and fits (for the WCS). Example URL for this script provided by John Wu.
    pixscale = 0.262   # default arcsec/pixel
    dimen = (opt_view.to(u.arcsec).value / pixscale).astype(int)
    if dimen[0] > 3000:
        pixscale = pixscale*dimen[0]/3000.
        dimen = [3000]

    if decals == 'decaps':
        url = 'https://legacysurvey.org/viewer/jpeg-cutout/?ra={}&dec={}&layer=decaps2&' \
              'pixscale={}&width={}&height={}&bands=g'.format(hi_pos.ra.deg, hi_pos.dec.deg, pixscale, dimen[0], dimen[-1])
    elif decals == 'dr9':
        url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra={}&dec={}&layer=ls-dr9&' \
              'pixscale={}&width={}&height={}&bands=g'.format(hi_pos.ra.deg, hi_pos.dec.deg, pixscale, dimen[0], dimen[-1])
    elif decals == 'sdss':
        url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra={}&dec={}&layer=sdss&' \
              'pixscale={}&width={}&height={}&bands=g'.format(hi_pos.ra.deg, hi_pos.dec.deg, pixscale, dimen[0], dimen[-1])
    else:
        url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra={}&dec={}&layer=ls-dr10&' \
              'pixscale={}&width={}&height={}&bands=g'.format(hi_pos.ra.deg, hi_pos.dec.deg, pixscale, dimen[0], dimen[-1])

    try:
        fits_head = fits.getheader(url)
        r = requests.get(url.replace("fits", "jpg").split('bands')[0])
        color_im = Image.open(BytesIO(r.content))
    except HTTPError:
        if decals == 'decaps':
            logger.warning("\tHTTP Error, no DECaPS false color image retrieved. Server error or no DECaPS coverage?")
        elif decals == 'dr9':
            logger.warning("\tHTTP Error, no DECaLS false color image retrieved. Server error or no DECaLS DR9 coverage?")
        else:
            logger.warning("\tHTTP Error, no DECaLS false color image retrieved. Server error or no DECaLS DR10 coverage?")
        fits_head = None
        color_im = None
    except URLError:
        logger.warning("\tURL Error. No false color image retrieved.")
        fits_head = None
        color_im = None

    return color_im, fits_head


def get_wise(hi_pos, opt_view=6*u.arcmin, survey='WISE W1'):
    """Retrieve a WISE image from a certain pointing.

    :param hi_pos: position in the radio spectral line data data
    :type hi_pos: astropy SkyCoord object
    :param opt_view: size of the optical image, defaults to 6*u.arcmin
    :type opt_view: astropy.units.arcmin, optional
    :param survey: survey containing WISE data, defaults to 'WISE W1'
    :type survey: str, optional
    :return: optical image
    :rtype: astropy HDUList
    """

    # IRSA query only accepts things in J2000/ICRS, so:
    # Note 2MASS is available through astroquery! get_irsa --> get_wise
    hi_pos = hi_pos.transform_to('icrs')
    params = {'POS': '{},{}'.format(hi_pos.ra.deg, hi_pos.dec.deg)}

    api_endpoint = "https://irsa.ipac.caltech.edu/ibe/search/wise/allwise/p3am_cdd"
    
    try:
        r = requests.get(url=api_endpoint, params=params)
        tab = ascii.read(r.content.decode(), header_start=44, data_start=48, format='fixed_width')
        params['coadd_id'] = tab['coadd_id'][0]
        params['coaddgrp'] = tab['coadd_id'][0][:2]
        params['coadd_ra'] = tab['coadd_id'][0][:4]
        params['band'] = int(survey[-1])
        params['size'] = str(opt_view[0].value) + ',' + str(opt_view[-1].value) + str(opt_view[0].unit)
        path = str.format("/{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id:s}-w{band:1d}-int-3.fits?"
                        "center={POS:s}&size={size:s}", **params)
    except requests.exceptions.ConnectionError:
        path = []

    try:
        result = fits.open(api_endpoint.replace('search', 'data') + path)
    except:
        result = None

    return result
