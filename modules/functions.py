from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
import numpy as np


def chan2freq(channels, fits_name):
    header = fits.getheader(fits_name)
    frequencies = (channels * header['CDELT3'] + header['CRVAL3']) * u.Hz
    return frequencies


def chan2vel(channels, fits_name):
    header = fits.getheader(fits_name)
    # Need to deal with different types of headers and velocity scaling with or without frequency!!! (Also bary vs topo, etc)
    velocities = (channels * header['CDELT3'] + header['CRVAL3']) * u.m / u.s
    return velocities


def get_info(fits_name):

    # For FITS conventions on the equinox, see:
    # https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf

    header = fits.getheader(fits_name)

    try:
        bmaj = header['BMAJ'] * 3600. * u.arcsec
        bmin = header['BMIN'] * 3600. * u.arcsec
        bpa = header['BPA']
        cellsize = header['CDELT2'] * 3600. * u.arcsec
        pix_per_beam = bmaj / cellsize * bmin / cellsize * np.pi / (4 * np.log(2))
    except:
        print("\tWARNING: Couldn't find beam in primary header information; in other extension? " \
              "Assuming beam is 3.5x3.5 pixels")
        cellsize = header['CDELT2'] * 3600. * u.arcsec
        bmaj, bmin, bpa = 3.5 * cellsize, 3.5 * cellsize, 0
        pix_per_beam = bmaj / cellsize * bmin / cellsize * np.pi / (4 * np.log(2))

    chan_width = header['CDELT3']
    if 'FREQ' in header['CTYPE3']:
        units = u.Hz
    else:
        units = u.m / u.s
    chan_width = chan_width * units

    # Add code to deal with reference frame?  AIPS conventions use VELREF: http://parac.eu/AIPSMEM117.pdf
    # spec_sys = header['SPECSYS']

    try:
        equinox = header['EQUINOX']
        if equinox < 1984.0:
            equinox = 'B' + str(equinox)
            frame = 'fk4'
        else:
            equinox = 'J' + str(equinox)
            frame = 'fk5'
    except:
        print("\tWARNING: No equinox information in header; assuming ICRS frame.")
        equinox = None
        frame = 'icrs'

    return {'bmaj': bmaj, 'bmin': bmin, 'bpa': bpa, 'pix_per_beam': pix_per_beam,
            'chan_width': chan_width, 'equinox': equinox, 'frame': frame, 'cellsize': cellsize}


def get_radecfreq(catalog, original):

    header = fits.getheader(original)
    wcs = WCS(header)
    Xc = catalog['x']
    Yc = catalog['y']
    if header['NAXIS'] == 3:
        subcoords = wcs.wcs_pix2world(Xc, Yc, 1, 0)   # origin follows: spatial, spectral, stokes?
    if header['NAXIS'] == 4:
        subcoords = wcs.wcs_pix2world(Xc, Yc, 1, 0, 0)
    ra, dec, freq = subcoords[0], subcoords[1], subcoords[2]

    return ra, dec, freq


def get_subcube(source, original):

    hdu_orig = fits.open(original)

    if hdu_orig[0].header['NAXIS'] == 4:
        stokes_dim, z_dim, y_dim, x_dim = 0, 1, 2, 3
    if hdu_orig[0].header['NAXIS'] == 3:
        z_dim, y_dim, x_dim = 0, 1, 2

    # Some lines stolen from cubelets in  SoFiA:
    # Could consider allowing a user specified range in z.
    cubeDim = hdu_orig[0].data.shape
    Xc = source['x']
    Yc = source['y']
    Xmin = source['x_min']
    Ymin = source['y_min']
    Xmax = source['x_max']
    Ymax = source['y_max']
    cPixXNew = int(Xc)
    cPixYNew = int(Yc)
    maxX = 2 * max(abs(cPixXNew - Xmin), abs(cPixXNew - Xmax))
    maxY = 2 * max(abs(cPixYNew - Ymin), abs(cPixYNew - Ymax))
    XminNew = cPixXNew - maxX
    if XminNew < 0: XminNew = 0
    YminNew = cPixYNew - maxY
    if YminNew < 0: YminNew = 0
    XmaxNew = cPixXNew + maxX
    if XmaxNew > cubeDim[x_dim] - 1: XmaxNew = cubeDim[x_dim] - 1
    YmaxNew = cPixYNew + maxY
    if YmaxNew > cubeDim[y_dim] - 1: YmaxNew = cubeDim[y_dim] - 1

    if len(cubeDim) == 4:
        subcube = hdu_orig[0].data[0, :, int(YminNew):int(YmaxNew) + 1, int(XminNew):int(XmaxNew) + 1]
    elif len(cubeDim) == 3:
        subcube = hdu_orig[0].data[:, int(YminNew):int(YmaxNew) + 1, int(XminNew):int(XmaxNew) + 1]
    else:
        print("WARNING: Original cube does not have 3-4 dimensions.")
        subcube = None

    hdu_orig.close()

    return subcube
