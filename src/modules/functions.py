from astropy.coordinates import SkyCoord
from astropy import constants as const
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
import numpy as np
from pvextractor import extract_pv_slice, PathFromCenter

from src.modules.logger import Logger

logger = Logger.get_logger()

HI_restfreq = 1420405751.77 * u.Hz


###################################################################


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
    # Don't know how to deal with cubelets having diff CRPIX3 from orig data; catalog data is in ref to orig base 0
    # frequencies = (header['CDELT3'] * (channels - (header['CRPIX3'] - 1)) + header['CRVAL3']) * u.Hz  # parameter.offset=True ?
    frequencies = (header['CDELT3'] * channels + header['CRVAL3']) * u.Hz     # parameter.offset=False ?
    # frequencies = (header['CDELT3'] * (channels - (header['CRPIX3']-1)) + header['CRVAL3']) * u.Hz
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
    logger.warning("\tAssuming channels are uniform width in velocity.")
    header = fits.getheader(fits_name)
    # Don't know how to deal with cubelets having diff CRPIX3 from orig data; catalog data is in ref to orig base 0
    velocities = (header['CDELT3'] * (channels - (header['CRPIX3'] - 1)) + header['CRVAL3']) * u.m / u.s
    # velocities = (header['CDELT3'] * channels + header['CRVAL3']) * u.m / u.s
    return velocities


def felo2vel(channels, fits_name):
    """Converts channels to velocities for a cube with non-linear channels.

    N.B.: This conversion differs from the output of SoFiA-2 which uses wcslib and therefore may not be strictly correct.

    :param channels:
    :type channels: Iterable[int]
    :param fits_name: name of the FITS file
    :type fits_name: str
    :return: calculated velocities
    :rtype: Iterable[float]
    """
    # Formula taken from here: https://www.astro.rug.nl/software/kapteyn/spectralbackground.html#aips-axis-type-felo
    logger.warning("\tAxis type FELO...this conversion may not be precise (may be off by ~10 km/s).")
    c = const.c.to(u.m/u.s).value
    header = fits.getheader(fits_name)
    fr = header['RESTFREQ'] / (1 + header['CRVAL3'] / c)
    df = -1 * header['RESTFREQ'] * header['CDELT3'] * c / ((c + header['CRVAL3']) * (c + header['CRVAL3']))
    velocities = header['CRVAL3'] + c * header['RESTFREQ'] * (1 / (fr + (channels - header['CRPIX3']) * df) - 1 / fr)
    return velocities


def sbr2nhi(sbr, bunit, bmaj, bmin, source, spec_line=None, Jykms=False):
    """Get the HI column density from sbr.  See Section 15 of Meyer et al (2017) for equations: 
    https://ui.adsabs.harvard.edu/abs/2017PASA...34...52M/abstract

    :param sbr: SBR
    :type sbr: float
    :param bunit: unit in which sbr is measured
    :type bunit: str
    :param bmaj: major axis of the beam in arcseconds
    :type bmaj: float
    :param bmin: minor axis of the beam
    :type bmin: float
    :param source: source object
    :type source: Astropy table
    :param spec_line: spectral line properties
    :type spec_line: dict
    :return: column density
    :rtype: float
    """

    c = const.c.to(u.m/u.s).value
    spec_line = spec_line['name']
    
    if (bunit == 'Jy/beam*m/s') or (bunit == 'Jy/beam*M/S'):
        logger.warning("\tAssumes velocity axis of cube is in the *observed* frame. If cube is in source rest frame, "
              "the column density is (1+z) times greater than shown.")
        if ('v_rad' in source.colnames): # or (cube_params['spec_axis'] == 'VRAD'): # Taken from make_images.py.
            # First convert to observed frequency, then to redshift following these equations:
            # https://web-archives.iram.fr/ARN/may95/node4.html
            logger.warning("\tHI column density calculation assumes optical velocity convention. Data is in radio convention!")
            vel_sys = source['v_rad']
            freq_sys = HI_restfreq * (1 - vel_sys/c)
            z = (HI_restfreq - freq_sys) / freq_sys
        elif 'v_opt' in source.colnames:
            vel_sys = source['v_opt']
            z = vel_sys / c
        elif 'v_app' in source.colnames:
            vel_sys = source['v_app']
            z = vel_sys / c

        if (spec_line == None) or (spec_line == 'HI'):
            nhi = 1.104e+21 * (1 + z) ** 2 * sbr / bmaj / bmin
        else:
            logger.warning("\tMom0 units not corrected for non-HI lines in units of Jy/beam*m/s.")
            nhi = sbr

    elif (bunit == 'Jy/beam*Hz') or (bunit == 'beam-1 Jy*Hz'):
        if (spec_line == None) or (spec_line == 'HI'):
            freq_sys = source['freq']
            z = (HI_restfreq.value - freq_sys) / freq_sys
            nhi = 2.330e+20 * (1 + z) ** 4 * sbr / bmaj / bmin
        else:
            # Units of Jy/beam Hz or Jy/beam km/s depending on what was done earlier
            nhi = sbr
            # Units of Kelvin Hz
            # nhi = 1.222e+3 * sbr / bmaj / bmin / line['restfreq'].to(u.GHz).value**2
        
    else:
        if (spec_line == None) or (spec_line == 'HI'):
            logger.warning("\tMom0 image units are not Jy/beam*m/s or Jy/beam*Hz. Cannot convert to HI column density.")
        else:
            logger.warning("\tMom0 image units are not Jy/beam*m/s or Jy/beam*Hz.")
        nhi = sbr
    
    if np.isfinite(nhi):
        nhi_ofm = int(np.floor(np.log10(np.abs(nhi))))
    else:
        nhi_ofm = 0
    
    if (spec_line == None) or (spec_line == 'HI'):
        nhi_label = r'$N_\mathrm{{HI}}$ = {0:.1f} x $10^{{ {1:d} }}$ cm$^{{-2}}$'.format(nhi/10**nhi_ofm, nhi_ofm)
        nhi_labels = r'$N_\mathrm{{HI}}$ = $2^n$ x {0:.1f} x $10^{{ {1:d} }}$ cm$^{{-2}}$ ($n$=0,1,...)'.format(nhi/10**nhi_ofm, nhi_ofm)
    else:
        if spec_line == 'Unknown':
            subscript = ''
        else:
            subscript = spec_line
        if Jykms == False:
            nhi_label = r'$S_\mathrm{{{0:s}}}$ = {1:.1f} x $10^{{ {2:d} }}$ Jy/bm Hz'.format(spec_line, nhi/10**nhi_ofm, nhi_ofm)
            nhi_labels = r'$S_\mathrm{{{:s}}}$ = $2^n$ x {:.1f}x$10^{{ {:d} }}$ Jy/bm Hz ($n$=0,1...)'.format(subscript,
                                                                                                              nhi/10**nhi_ofm, nhi_ofm)
        else:
            nhi_label = r'$S_\mathrm{{{0:s}}}$ = {1:.1f} x $10^{{ {2:d} }}$ Jy/bm km/s'.format(spec_line, nhi/10**nhi_ofm, nhi_ofm)
            nhi_labels = r'$S_\mathrm{{{:s}}}$ = $2^n$ x {:.1f}x$10^{{ {:d} }}$ Jy/bm km/s ($n$=0,1...)'.format(subscript,
                                                                                                              nhi/10**nhi_ofm, nhi_ofm)

    return nhi, nhi_label, nhi_labels


def get_info(fits_name, beam=None, source_id=0):
    """Get the beam info from a FITS file.

    :param fits_name: name of the FITS file
    :type fits_name: str
    :param beam: beam specifications, defaults to None. Specifications are
        given in arcsec (axes) and degrees (position_angle), and formatted as
        {[major_axis, minor_axis, position_angle]|[major_axis, minor_axis]|
        [position_angle]}
    :type beam: Iterable[float], optional
    :return: The characteristics of the beam and coordinate system of the image.
    :rtype: dict
    """

    # For FITS conventions on the equinox, see:
    # https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf

    header = fits.getheader(fits_name)

    cellsize = header['CDELT2'] * 3600. * u.arcsec

    default_beam = False
    if len(beam) == 3:
        logger.info(f"\t\tUsing user specified beam: {beam[0]} arcsec by {beam[1]} arcsec; PA: {beam[2]} deg")
        bmaj = beam[0] * u.arcsec
        bmin = beam[1] * u.arcsec
        bpa = beam[2]
    elif len(beam) == 2:
        logger.warning(f"\tAssuming PA = 0. Using user specified beam: {beam[0]} arcsec by {beam[1]} arcsec.")
        bmaj = beam[0] * u.arcsec
        bmin = beam[1] * u.arcsec
        bpa = 0
    elif len(beam) == 1:
        logger.warning(f"\tUsing user specified circular beam size of {beam[0]} arcsec.")
        bmaj = bmin = beam[0] * u.arcsec
        bpa = 0
    else:
        try:
            bmaj = header['BMAJ'] * 3600. * u.arcsec
            bmin = header['BMIN'] * 3600. * u.arcsec
            bpa = header['BPA']
            if bmaj * bmin == 0.0:
                logger.warning("\tBMAJ and/or BMIN in header = 0! Assuming beam is 3.5x3.5 pixels")
                logger.warning("\t\tColumn density and beam plotted as order of magnitude estimate ONLY. ")
                logger.warning("\t\tRerun with -b and provide beam info to remove red strikethroughs on plots.")
                bmaj, bmin, bpa = 3.5 * cellsize, 3.5 * cellsize, 0
                default_beam = True
            else:
                logger.info(f"\t\tFound {bmaj:.1f} by {bmin:.1f} beam with PA={bpa:.1f} deg in primary header.")
        except:
            logger.warning("\tCouldn't find beam in primary header information; in other extension? "
                  "Assuming beam is 3.5x3.5 pixels")
            logger.warning("\t\tColumn density and beam plotted as order of magnitude estimate ONLY. ")
            logger.warning("\t\tRerun with -b and provide beam info to remove red strikethroughs on plots.")
            bmaj, bmin, bpa = 3.5 * cellsize, 3.5 * cellsize, 0
            default_beam = True

    pix_per_beam = bmaj / cellsize * bmin / cellsize * np.pi / (4 * np.log(2))

    # Try catching cubes in Galactic coordinates first
    if 'GLON' in header['CTYPE1']:
        logger.info("\t\tFound data is in Galactic spatial frame.")
        equinox = None
        frame = 'galactic'
    # If not Galacticc, try to determine the equinox of the observations
    else:
        try:
            equinox = header['EQUINOX']
            if equinox < 1984.0:
                equinox = 'B' + str(equinox)
                frame = 'fk4'
            else:
                equinox = 'J' + str(equinox)
                frame = 'fk5'
            logger.info("\t\tFound {} equinox in header.".format(equinox))
        except KeyError:
            try:
                equinox = header['EPOCH']
                if equinox < 1984.0:
                    equinox = 'B' + str(equinox)
                    frame = 'fk4'
                else:
                    equinox = 'J' + str(equinox)
                    frame = 'fk5'
                logger.warning("\tUsing deprecated EPOCH in header for equinox: {}.".format(equinox))
            except KeyError:
                logger.warning("\tNo equinox information in header; assuming ICRS frame.")
                equinox = None
                frame = 'icrs'

    # Try to determine the reference frame.  AIPS conventions use VELREF: http://parac.eu/AIPSMEM117.pdf
    spec_sys = False
    try:
        spec_sys = header['SPECSYS']
        logger.info("\t\tFound {} reference frame specified in SPECSYS in header.".format(spec_sys))
    except:
        try:
            spec_sys = header['SPECSYS3']
            logger.info("\t\tFound {} reference frame specified in SPECSYS3 in header.".format(spec_sys))
        except:
            try:
                velref = header['VELREF']
                if velref == 1: spec_sys = 'LSR'
                if velref == 2: spec_sys = 'HELIOCEN'
                if velref == 3: spec_sys = 'TOPOCENT'
                logger.info("\t\tDerived {} reference frame from VELREF in header using AIPS convention.".format(spec_sys))
            except:
                # Comment this message out for now...program checks later.
                logger.info("\t\tNo SPECSYS, SPECSYS3, or VELREF in header to define reference frame, will check CTYPE3.")
                pass

    # Try to determine the spectral properties
    if ((fits_name[-9:] != 'cube.fits') and (source_id != 0)) or ('mom' in fits_name):
        logger.warning("\tRetrieving info from a moment map or other 2D image?")
        chan_width = None
        spec_axis = None
    else:
        spec_axis = header['CTYPE3']
        chan_width = header['CDELT3']
        if 'FREQ' in spec_axis:
            units = u.Hz
        else:
            units = u.m / u.s
        chan_width = chan_width * units

        logger.info("\t\tFound CTYPE3 spectral axis type {} in header.".format(spec_axis))
        if ("-" in spec_axis) and spec_sys:
            logger.warning("\tDropping end of spectral axis type. Using SPECSYS/VELREF for reference frame.")
            spec_axis = spec_axis.split ("-")[0]
        elif ("-" in spec_axis) and (not spec_sys):
            spec_sys = spec_axis.split("-")[1]
            spec_axis = spec_axis.split("-")[0]
            if spec_sys == 'HEL': spec_sys = 'HELIOCEN'
            logger.warning("\tAttempting to use end of CTYPE3 for reference frame: {}".format(spec_sys))

    if not spec_sys:
        logger.info("\t\tNo SPECSYS, VELREF, or reference frame in CTYPE3, assuming data in TOPOCENT reference frame.")
        spec_sys = 'TOPOCENT'

    return {'bmaj': bmaj, 'bmin': bmin, 'bpa': bpa, 'pix_per_beam': pix_per_beam, 'default_beam': default_beam,
            'chan_width': chan_width, 'equinox': equinox, 'frame': frame, 'cellsize': cellsize, 'spec_sys': spec_sys,
            'spec_axis': spec_axis}


## Get uniform wcs object
def get_wcs_info(fits_name):
    """Get uniform wcs object for plotting

    :param fits_name: basename for the source for data files
    :type fits_name: str
    :return: world coordinate system object
    :rtype: WCS object
    """
    cubeh = fits.getheader(fits_name)
    cubew = WCS(cubeh, fix=True, translate_units='shd').sub(2)

    if cubew.wcs.equinox != 2000.0:
        sky = cubew.sub(2).pixel_to_world(cubew.wcs.crpix[0], cubew.wcs.crpix[1])
        sky_icrs = sky.transform_to("icrs")
        
        cubeh['EPOCH'] = 2000.0
        cubeh['CRVAL1'] = sky_icrs.ra.deg
        cubeh['CRVAL2'] = sky_icrs.dec.deg
        hiwcs = WCS(cubeh).sub(2)
    else:
        hiwcs = cubew

    return hiwcs, cubew


def get_radecfreq(catalog, original):
    """Get the right ascension, declination, and frequeny of a catalog object.

    :param catalog: catalog object header
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
    :type source: Astropy table
    :param original: original data file
    :type original: str
    :return: subcube of data
    :rtype: NDArray
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
        logger.warning("Original cube does not have 3-4 dimensions.")
        subcube = None

    hdu_orig.close()

    return subcube


def create_pv(source, filename, opt_view=6*u.arcmin, min_axis=False):
    """

    :param source: source object
    :type source: Astropy table
    :param filename: name of FITS file
    :type filename: str
    :param opt_view: requested size of the image for regriding
    :type opt_view: quantity
    :param min_axis: flag for extracting major or minor axis
    :type min_axis: boolean
    :return: position-velocity slice of the mask cube
    :rtype: FITS HDU
    """

    pos_angle = source['kin_pa']
    if min_axis == True:
        pos_angle += 90.
    slice = PathFromCenter(center=SkyCoord(ra=source['pos_x'], dec=source['pos_y'], unit='deg'),
                           length=opt_view, angle=pos_angle*u.deg, width=6*u.arcsec)
    mask = fits.open(filename)
    try:
        mask_pv = extract_pv_slice(mask[0].data, slice, wcs=WCS(mask[0].header, fix=True, translate_units='shd'))
    except ValueError:
        logger.warning('\tpvextractor is complaining about non-square pixels, try with assert_square = False')
        try:
            mask_pv = extract_pv_slice(mask[0].data, slice, wcs=WCS(mask[0].header, fix=True, translate_units='shd'),
                                       assert_square=False)
        except:
            logger.error('\tCannot extract pv slice of mask. Try upgrading to latest version of pvextractor (v>=0.4) from github:')
            logger.error('\t\t"python3 -m pip install git+https://github.com/radio-astro-tools/pvextractor"')
            mask_pv = None
    mask.close()

    return mask_pv


def line_lookup(spec_line):
    """Return rest frequency and object to convert to optical velocities for requested line.
    https://www.narrabri.atnf.csiro.au/observing/spectral.html
    https://splatalogue.online/#/basic

    :param spec_line:
    :type spec_line:
    :return:
    :rtype:
    """

    rad_opt = 'Optical'
    if (spec_line == '') or (spec_line == 'HI'):
        spec_line == 'HI'
        restfreq_line = 1420405751.77 * u.Hz
        convention = u.doppler_optical(restfreq_line)
    elif spec_line == 'CO(1-0)':
        restfreq_line = 115.27120180 * u.GHz
        convention = u.doppler_optical(restfreq_line)
    elif spec_line == 'CO(2-1)':
        restfreq_line = 230.538 * u.GHz
        convention = u.doppler_optical(restfreq_line)
    elif spec_line == 'CO(3-2)':
        restfreq_line = 345.79599 * u.GHz
        convention = u.doppler_optical(restfreq_line)
    elif spec_line == 'OH_1612':
        restfreq_line = 1.6122309 * u.GHz
        convention = u.doppler_optical(restfreq_line)
    elif spec_line == 'OH_1665':
        restfreq_line = 1.6654018 * u.GHz
        convention = u.doppler_optical(restfreq_line)
    elif spec_line == 'OH_1667':
        restfreq_line = 1.6673590 * u.GHz
        convention = u.doppler_optical(restfreq_line)
    elif spec_line == 'OH_1720':
        restfreq_line = 1.7205299 * u.GHz
        convention = u.doppler_optical(restfreq_line)
    else:
        logger.error("\t\tUnrecognized spectral line, try 'HI'; 'CO(1-0)', etc; or 'OH_1667', etc. " \
                     "Continuing to make plots without a known redshift.")
        spec_line = 'Unknown'
        restfreq_line = None
        convention = None

    return {'name': spec_line, 'restfreq': restfreq_line, 'convention':convention, 'rad_opt':rad_opt}


def plot_labels(source, ax, default_beam, x_color='k', id_label=''):
    """Plot labels on spatial plots depending on the coordinate frame.

    :param source: source object
    :type source: Astropy table
    :param ax: matplotlib axes instance
    :type ax: axes object
    :param default_beam: whether the synthesized beam is known from data/user or not
    :type default_beam: bool
    :param x_color: color of galaxy position marker
    :type x_color: str
    :return:
    """

    if 'l' in source.colnames:
        x_coord, y_coord = 'glon', 'glat'
        # x_label, y_label = 'Galactic Longitude [deg]', 'Galactic Latitude [deg]'
        x_label, y_label = r'$\it{{l}}$ [deg]', r'$\it{{b}}$ [deg]'
    else:
        x_coord, y_coord = 'ra', 'dec'
        x_label, y_label = 'RA (ICRS)', 'Dec (ICRS)'

    if source['id'] != 0:
        ax.scatter(source['pos_x'], source['pos_y'], marker='x', c=x_color, linewidth=0.75,
                transform=ax.get_transform('world'))
    ax.set_title(source['name'] + id_label, fontsize=24)
    ax.tick_params(axis='both', which='major', labelsize=22, length=6, width=2)
    ax.coords[x_coord].set_axislabel(x_label, fontsize=22)
    ax.coords[y_coord].set_axislabel(y_label, fontsize=22)
    if default_beam:
        ax.scatter(0.92, 0.9, marker='x', c='red', s=500, linewidth=5, transform=ax.transAxes, zorder=99)
        ax.plot([0.1, 0.9], [0.05, 0.05], c='red', linewidth=3, transform=ax.transAxes, zorder=100)
        ax.text(0.5, 0.5, 'Not calculated with correct beam', transform=ax.transAxes, fontsize=40, color='gray',
                alpha=0.5, ha='center', va='center', rotation=30, zorder=101)

    return


def make_hist_arr(xx=None, yy=None, yy_err=None):
    """Create arrays for plotting spectra as with square channel-like response.

    :param x: x-axis parameter
    :type x: array
    :param y: y-axis parameter
    :type y: array
    :param y_err: y-axis error
    :type y_err: array
    :return:
    """

    i = 0
    new_x, new_y, new_y_error = [], [], []
    delta_xx = xx[i] - xx[i+1]
    for x, y, y_err in zip(xx, yy, yy_err):
        new_x.append(x + delta_xx / 2)
        new_x.append(x)
        new_x.append(x - delta_xx / 2)
        new_y.append(y)
        new_y.append(y)
        new_y.append(y)
        new_y_error.append(0)
        new_y_error.append(y_err)
        new_y_error.append(0)
        i += 0
        if x != xx[-1]:
            delta_xx = xx[i] - xx[i+1]

    return new_x, new_y, new_y_error


def make_header(source, opt_view=6*u.arcmin):
    """Return a dummy header when an ancillary image can't be retrieved.

    :param source: source object
    :type source: Astropy table
    :param opt_view: requested size of the image for regriding
    :type opt_view: quantity
    :return:
    :rtype: FITS HDU
    """

    npix = 16
    hdu = fits.ImageHDU()
    hdu.data = np.ones([npix, npix])
    hdu.header['CTYPE1']  = 'RA---SIN'
    hdu.header['CRPIX1'] = npix / 2 + 1
    hdu.header['CRVAL1'] = source['pos_x']
    hdu.header['CDELT1'] = -1 * opt_view[0].to(u.deg).value / npix
    hdu.header['CTYPE2']  = 'DEC--SIN'
    hdu.header['CRPIX2'] = npix / 2 + 1
    hdu.header['CRVAL2'] = source['pos_y']
    hdu.header['CDELT2'] = opt_view[0].to(u.deg).value / npix

    return hdu.header


def add_source(catalog, fits_name):
    """Return a dummy source when a summary image is requested.

    :param catalog: source catalog
    :type catalog: Astropy table
    :param fits_name: name of the FITS file
    :type fits_name: str
    :return:
    :rtype: Astropy table
    """

    catalog.add_row(catalog[0])
    header = fits.getheader(fits_name)
    wcs = WCS(header, fix=True, translate_units='shd')

    # Change the relevant catalog parameters ... prob need to deal with kin_pa and rms at some point.
    catalog[-1]['name'] = fits_name.split('/')[-1][:-10]
    catalog[-1]['id'] = 0
    catalog[-1]['x'], catalog[-1]['y'], catalog[-1]['z'] = header['NAXIS1']/2, header['NAXIS2']/2, 0
    catalog[-1]['x_min'], catalog[-1]['x_max'] = int(np.min(catalog['x_min'])), int(np.max(catalog['x_max']))
    catalog[-1]['y_min'], catalog[-1]['y_max'] = int(np.min(catalog['y_min'])), int(np.max(catalog['y_max']))
    catalog[-1]['z_min'], catalog[-1]['z_max'] = int(np.min(catalog['z_min'])), int(np.max(catalog['z_max']))
    catalog[-1]['rms'] = np.nanmin(catalog['rms'])
    catalog[-1]['ra'], catalog[-1]['dec'] = wcs.celestial.wcs_pix2world(catalog[-1]['x'], catalog[-1]['y'], 0)
    catalog[-1]['pos_x'], catalog[-1]['pos_y'] = wcs.celestial.wcs_pix2world(catalog[-1]['x'], catalog[-1]['y'], 0)
    if 'freq' in catalog.colnames:
        catalog[-1]['freq'] = (np.min(catalog['freq']) + np.max(catalog['freq'])) / 2
    elif 'vrad' in catalog.colnames:
        catalog[-1]['v_rad'] = (np.min(catalog['v_rad']) + np.max(catalog['v_rad'])) / 2
    elif 'v_opt' in catalog.colnames:
        catalog[-1]['v_opt'] = (np.min(catalog['v_opt']) + np.max(catalog['v_opt'])) / 2
    elif 'v_app' in catalog.colnames:
        catalog[-1]['v_app'] = (np.min(catalog['v_app']) + np.max(catalog['v_app'])) / 2

    return catalog
