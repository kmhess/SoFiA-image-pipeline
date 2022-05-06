import requests
from PIL import Image
from io import BytesIO
from urllib.error import HTTPError

from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astroquery.skyview import SkyView
import numpy as np

from src.modules.panstarrs_fcns import getcolorim, geturl


def get_skyview(hi_pos, opt_view=6*u.arcmin, survey='DSS2 Blue'):
    """Retrieve the optical image from a certain pointing.

    :param hi_pos: position in the HI data
    :type hi_pos: astropy SkyCoord object
    :param opt_view: size of the optical image, defaults to 6*u.arcmin
    :type opt_view: astropy.units.arcmin, optional
    :param survey: survey containing optical data, defaults to 'DSS2 Blue'
    :type survey: str, optional
    :return: optical image
    :rtype: astropy HDUList
    """
    # DSS2 Blue images have a 1 arc/pix pixel scale, but retrieving ~the pixel scale doesn't work.
    opt_pixels = (opt_view.to(u.arcsec).value * 2).astype(int)

    # Get a survey image from SkyView:
    if hi_pos.frame.name == 'galactic':
        path = SkyView.get_images(position=hi_pos.to_string('hmsdms'), coordinates='galactic',
                                  width=opt_view[0], height=opt_view[-1], survey=[survey], pixels=opt_pixels)
    elif (not hi_pos.equinox) or (hi_pos.frame.name == 'icrs'):
        path = SkyView.get_images(position=hi_pos.to_string('hmsdms'), coordinates='ICRS',
                                  width=opt_view[0], height=opt_view[-1], survey=[survey], pixels=opt_pixels)
    # Note that there seems to be a bug in SkyView that it sometimes won't retrieve non-J2000.0.  Keep an eye on this!
    else:
        path = SkyView.get_images(position=hi_pos.to_string('hmsdms'), coordinates=hi_pos.equinox.value,
                                  width=opt_view[0], height=opt_view[-1], survey=[survey], pixels=opt_pixels)
    if len(path) != 0:
        print("\tSurvey image retrieved from {}.".format(survey))
        result = path[0]
    else:
        print("\tWARNING: No {} image retrieved.  Bug, or server error?  Try again later?".format(survey))
        result = None

    return result


def get_panstarrs(hi_pos, opt_view=6*u.arcmin):
    """Get PanSTARRS false color image and r-band fits (for the WCS).

    :param hi_pos: position in the HI data
    :type hi_pos: astropy SkyCoord object
    :param opt_view: size of the optical image, defaults to 6*u.arcmin
    :type opt_view: astropy.units.arcmin, optional
    :return: color image and FITS header
    :rtype: Tuple[color_im, FITS_header] TODO check exact types!
    """
    #  Get PanSTARRS false color image and r-band fits (for the WCS).
    pstar_pixsc = 0.25
    if len(opt_view) > 1:
        print("\tWARNING: PanSTARRS only returns square images; taking largest dimension.")
        opt_view = np.max(opt_view)
    path = geturl(hi_pos.ra.deg, hi_pos.dec.deg, size=int(opt_view.to(u.arcsec).value / pstar_pixsc),
                  filters="r", format="fits")

    if len(path) != 0:
        fits_head = fits.getheader(path[0])
        color_im = getcolorim(hi_pos.ra.deg, hi_pos.dec.deg, size=int(opt_view.to(u.arcsec).value / pstar_pixsc),
                              filters="gri")
        print("\tOptical false color image retrieved from PanSTARRS.")
    else:
        print("\tWARNING: No PanSTARRS false color image retrieved.  Server error or no PanSTARRS coverage?")
        fits_head = None
        color_im = None

    return color_im, fits_head


def get_decals(hi_pos, opt_view=6*u.arcmin):
    """Get DECaLS false color image and g-band fits (for the WCS).

    :param hi_pos: position in the HI data
    :type hi_pos: astropy SkyCoord object
    :param opt_view: size of the optical image, defaults to 6*u.arcmin
    :type opt_view: astropy.units.arcmin, optional
    :return: color image and FITS header
    :rtype: Tuple[color_im, FITS_header] TODO check exact types!
    """
    # Get DECaLS false color image and fits (for the WCS). Example URL for this script provided by John Wu.
    pixscale = 0.262   # default(?) arcsec/pixel
    dimen = (opt_view.to(u.arcsec).value / pixscale).astype(int)
    url = 'https://www.legacysurvey.org/viewer/cutout.fits?ra={}&dec={}&layer=ls-dr9&' \
          'pixscale={}&width={}&height={}&bands=g'.format(hi_pos.ra.deg, hi_pos.dec.deg, pixscale, dimen[0], dimen[-1])

    try:
        fits_head = fits.getheader(url)
        r = requests.get(url.replace("fits", "jpg").split('bands')[0])
        color_im = Image.open(BytesIO(r.content))
    except HTTPError:
        print("\tWARNING: HTTP Error, no DECaLS false color image retrieved. Server error or no DECaLS coverage?")
        fits_head = None
        color_im = None

    return color_im, fits_head
