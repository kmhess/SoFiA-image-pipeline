from modules.panstarrs_fcns import getcolorim, geturl

from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astroquery.skyview import SkyView


def get_skyview(hi_pos, opt_view=6*u.arcmin, survey='DSS2 Blue'):

    # DSS2 Blue images have a 1 arc/pix pixel scale, but retrieving ~the pixel scale doesn't work.
    opt_pixels = int(opt_view.to(u.arcsec).value * 2)

    # Get DSS2 Blue optical image:
    if (not hi_pos.equinox) or (hi_pos.frame.name == 'icrs'):
        path = SkyView.get_images(position=hi_pos.to_string('hmsdms'), coordinates='ICRS',
                                  width=opt_view, height=opt_view, survey=[survey], pixels=opt_pixels,
                                  cache=False)
    # Note that there seems to be a bug in SkyView that it sometimes won't retrieve non-J2000.0.  Keep an eye on this!
    else:
        path = SkyView.get_images(position=hi_pos.to_string('hmsdms'), coordinates=hi_pos.equinox.value,
                                  width=opt_view, height=opt_view, survey=[survey], pixels=opt_pixels,
                                  cache=False)
    if len(path) != 0:
        print("\tOptical image retrieved from {}.".format(survey))
        result = path[0]
    else:
        print("\tWARNING: No {} image retrieved.  Bug, or server error?  Try again later?".format(survey))
        result = None

    return result


def get_panstarrs(hi_pos, opt_view=6*u.arcmin):

    #  Get PanSTARRS false color image and r-band fits (for the WCS).
    pstar_pixsc = 0.25
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


