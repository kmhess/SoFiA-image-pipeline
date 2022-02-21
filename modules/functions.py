from astropy import constants as const
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
import numpy as np


def chan2freq(channels, fits_name):
    """Convert channels to frequencies.

    :param channels: which channels to convert
    :type channels: Iterable[int]
    :param fits_name: name of the FITS file
    :type fits_name: str
    :return: frequencies
    :rtype: Iterable[float]
    """
    header = fits.getheader(fits_name)
    frequencies = (header['CDELT3'] * (channels - (header['CRPIX3'] - 1)) + header['CRVAL3']) * u.m / u.s
    return frequencies


def chan2vel(channels, fits_name):
    """Convert channels to velocities.

    N.B.: This assumes the channels have uniform width in velocity space,
          which may not be the case!

    :param channels: the channels to convert
    :type channels: Iterable[int]
    :param fits_name: name of the FITS file
    :type fits_name: str
    :return: calculated velocities
    :rtype: Iterable[float]
    """
    print("\tWARNING: Assuming channels are uniform width in velocity.")
    header = fits.getheader(fits_name)
    # Need to deal with different types of headers and velocity scaling with or without frequency!!! (Also bary vs topo, etc)
    velocities = (header['CDELT3'] * (channels - (header['CRPIX3'] - 1)) + header['CRVAL3']) * u.m / u.s
    return velocities


def felo2vel(channels, fits_name):
    # Formula taken from here: https://www.astro.rug.nl/software/kapteyn/spectralbackground.html#aips-axis-type-felo
    print("\tWARNING: Axis type FELO...this conversion may not be precise (may be off by ~10 km/s).")
    c = const.c.to(u.m/u.s).value
    header = fits.getheader(fits_name)
    fr = header['RESTFREQ'] / (1 + header['CRVAL3'] / c)
    df = -1 * header['RESTFREQ'] * header['CDELT3'] * c / ((c + header['CRVAL3']) * (c + header['CRVAL3']))
    velocities = header['CRVAL3'] + c * header['RESTFREQ'] * (1 / (fr + (channels - header['CRPIX3']) * df) - 1 / fr)
    return(velocities)


def sbr2nhi(sbr, bunit, bmaj, bmin):
    """Get the HI column density from sbr.

    :param sbr: SBR
    :type sbr: float
    :param bunit: unit in which sbr is measured
    :type bunit: str
    :param bmaj: major axis of the beam
    :type bmaj: float
    :param bmin: minor axis of the bea,
    :type bmin: float
    :return: column density
    :rtype: float
    """
    if bunit == 'Jy/beam*m/s':
      nhi = 1.104e+21 * sbr / bmaj / bmin
    elif bunit == 'Jy/beam*Hz':
      nhi = 2.330e+20 * sbr / bmaj / bmin
    else:
      print("\tWARNING: Mom0 imag units are not Jy/beam*m/s or Jy/beam*Hz. Cannot convert to HI column density.")
      nhi = sbr
    return nhi


def get_info(fits_name, beam=None):
    """Get the beam info from a FITS file.

    :param fits_name: name of the FITS file
    :type fits_name: str
    :param beam: beam specifications, defaults to None. Specifications are
        given in arcsec (axes) and degrees (position_angle), and formatted as
        {[major_axis, minor_axis, position_angle]|[major_axis, minor_axis]|
        [position_angle]}
    :type beam: Iterable[float], optional
    :return: The characteristics of the beam
    :rtype: dict
    """

    # For FITS conventions on the equinox, see:
    # https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf

    header = fits.getheader(fits_name)

    cellsize = header['CDELT2'] * 3600. * u.arcsec

    if len(beam) == 3:
        print(f"\tUsing user specified beam: {beam[0]} arcsec by {beam[1]} arcsec; PA: {beam[2]} deg")
        bmaj = beam[0] * u.arcsec
        bmin = beam[1] * u.arcsec
        bpa = beam[2]
    elif len(beam) == 2:
        print(f"\tWARNING: assuming PA = 0. Using user specified beam: {beam[0]} arcsec by {beam[1]} arcsec.")
        bmaj = beam[0] * u.arcsec
        bmin = beam[1] * u.arcsec
        bpa = 0
    elif len(beam) == 1:
        print(f"\tWARNING: using user specified circular beam size of {beam[0]} arcsec.")
        bmaj = bmin = beam[0] * u.arcsec
        bpa = 0
    else:
        try:
            bmaj = header['BMAJ'] * 3600. * u.arcsec
            bmin = header['BMIN'] * 3600. * u.arcsec
            bpa = header['BPA']
        except:
            print("\tWARNING: Couldn't find beam in primary header information; in other extension? " \
                  "Assuming beam is 3.5x3.5 pixels")
            bmaj, bmin, bpa = 3.5 * cellsize, 3.5 * cellsize, 0

    pix_per_beam = bmaj / cellsize * bmin / cellsize * np.pi / (4 * np.log(2))

    chan_width = header['CDELT3']
    if 'FREQ' in header['CTYPE3']:
        units = u.Hz
    else:
        units = u.m / u.s
    chan_width = chan_width * units

    # Try to determine the reference frame.  AIPS conventions use VELREF: http://parac.eu/AIPSMEM117.pdf
    try:
        spec_sys = header['SPECSYS']
        print("\tFound {} reference frame specified in SPECSYS in header.".format(spec_sys))
    except:
        try:
            velref = header['VELREF']
            if velref == 1: spec_sys = 'LSR'
            if velref == 2: spec_sys = 'HELIOCEN'
            if velref == 3: spec_sys = 'TOPOCENT'
            print("\tDerived {} reference frame from VELREF in header using AIPS convention.".format(spec_sys))
        except:
            spec_sys = 'TOPOCENT'
            print("\tNo SPECSYS or VELREF in header, assuming data in TOPOCENT reference frame.")

    # Try to determine the spectral coordinates
    spec_axis = header['CTYPE3']
    print("\tFound spectral axis type {} in header.".format(spec_axis))
    if ("-" in spec_axis) and spec_sys:
        print("\tWARNING: dropping end of spectral axis type. Using SPECSYS/VELREF for reference frame.")
        spec_axis = spec_axis.split ("-")[0]
    elif ("-" in spec_axis) and (not spec_sys):
        print("\tWARNING: attempting to use end of spectral axis type for reference frame.")
        spec_axis = spec_axis.split("-")[0]
        spec_sys = spec_axis.split("-")[1]

    # Try to determine the equinox of the observations
    try:
        equinox = header['EQUINOX']
        if equinox < 1984.0:
            equinox = 'B' + str(equinox)
            frame = 'fk4'
        else:
            equinox = 'J' + str(equinox)
            frame = 'fk5'
        print("\tFound {} equinox in header.".format(equinox))
    except:
        try:
            equinox = header['EPOCH']
            if equinox < 1984.0:
                equinox = 'B' + str (equinox)
                frame = 'fk4'
            else:
                equinox = 'J' + str (equinox)
                frame = 'fk5'
            print("\tWARNING: Using deprecated EPOCH in header for equinox: {}.".format(equinox))
        except:
            print("\tWARNING: No equinox information in header; assuming ICRS frame.")
            equinox = None
            frame = 'icrs'

    return {'bmaj': bmaj, 'bmin': bmin, 'bpa': bpa, 'pix_per_beam': pix_per_beam, 'chan_width': chan_width,
            'equinox': equinox, 'frame': frame, 'cellsize': cellsize, 'spec_sys': spec_sys, 'spec_axis': spec_axis}


def get_radecfreq(catalog, original):
    """Get the right ascension, declination, and frequeny of a catalog object.

    :param catalog: catalog objet header
    :type catalog: astropy.Header? TODO check in function calls
    :param original: name of the original file
    :type original: str
    :return: right ascension, declination, and frequency
    :rtype: tuple
    """

    header = fits.getheader(original)
    wcs = WCS(header)
    # Get the x,y-position of the catalog object
    Xc = catalog['x']
    Yc = catalog['y']
    if header['NAXIS'] == 3:
        subcoords = wcs.wcs_pix2world(Xc, Yc, 1, 0)   # origin follows: spatial, spectral, stokes?
    if header['NAXIS'] == 4:
        subcoords = wcs.wcs_pix2world(Xc, Yc, 1, 0, 0)
    ra, dec, freq = subcoords[0], subcoords[1], subcoords[2]

    return ra, dec, freq


def get_subcube(source, original):
    """Retrieve a subcube from a datacube

    :param source: source object
    :type source: Astropy Header? TODO check!
    :param original: original data file
    :type original: str
    :return: subcube of data
    :rtype: NDArray? Or Astropy dataformat? TODO check
    """

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
