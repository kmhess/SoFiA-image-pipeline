#!/usr/bin/env python3

# Import default Python libraries
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import random
import string
import traceback

# Import installed Python libraries
from astropy.table import Table
from astropy import units as u
import importlib.metadata
import numpy as np

# Other imports from src are in code below so logger is initialized correctly
from src.modules.logger import Logger

version = importlib.metadata.version('SoFiA-image-pipeline')

###################################################################

def main():
    parser = ArgumentParser(description="Welcome to the SoFiA Image Pipeline, version {}.\n"
                            "Create images from a SoFiA catalog, and cubelets or fits file. Only works with SoFiA-2 and wcs=True (for now).".format(version),
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('-c', '--catalog', required=True,
                        help='Required: Specify the input XML or ascii catalog name. No default.')

    parser.add_argument('-id', '--source-id', default=[], nargs='*', #type=int,
                        help='Optional: Space-separated list, or range of sources to include in the plotting. If set to 0, do\n'
                             'summary figure. If set to -1, do summary figure and all sources. Default all sources.')

    parser.add_argument('-s', '--surveys', default=[], nargs='*', type=str,
                        help='Optional: Specify SkyView surveys to retrieve from astroquery on which to overlay HI contours.\n'
                             'These additional non-SkyView options are also available: \'decals\',\'decals-dr9\',\'decaps\',\n'
                             '\'sdss\', \'panstarrs\',\'hst\'. \'hst\' only refers to COSMOS HST. Default is "DSS2 Blue"\n' 
                             'if no user provided image. If \'none\' is requested, work in offline mode.')

    parser.add_argument('-ui', '--user-image', default=None,
                        help='Optional: Full path to the FITS image on which to overlay HI contours.')

    parser.add_argument('-ur', '--user-range', default=[10., 99.], nargs=2, type=float,
                        help='Optional: Percentile range used when displaying the user image (see "-ui"). Default is [10,99].')

    parser.add_argument('-line', '--spectral-line', default="HI", type=str,
                        help='Optional: Provide name of spectral line such as "HI"; "CO(1-0)" up to (3-2); or "OH_1667" or\n'
                             'other L-band accessible OH lines. Default is "HI". See github for more details. Work in progress.')

    parser.add_argument('-i', '--image_size', default=[6], nargs=1, type=float,
                        help='Optional: specify the minimum survey image size to retrieve in arcmin.  It will be adjusted if\n'
                             'the mom0 mask is larger. Note max panstarrs image size is 8 arcmin (default: %(default)s).')

    parser.add_argument('-snr', '--snr-range', default=[2., 3.], nargs=2, type=float,
                        help='Optional: specify which SNRmin and SNRmax values should be used to set the lowest reliable \n'
                             'contour in the figures. The contour level is calculated as the median value in the mom0 image\n'
                             'of all pixels whose SNR value is within the given range. Default is [2,3].')

    parser.add_argument('-o', '--original', default=None,
                        help='Optional: specify the original fits data: used for plotting aperture spectra over the full \n'
                             'frequency range of the cube. Otherwise, plot with noise over frequency range in the cubelet. \n'
                             'Uses 2D mask to integrate. (No default).')

    parser.add_argument('-b', '--beam', default=None,
                        help='Optional: specify the beam dimensions as "(bmaj,bmin,bpa)" in arcsec, arcsec, deg. If only \n'
                             '1 value is given, assume a circular beam. If 2 values are given, assume PA = 0. (No default).')

    parser.add_argument('-cw', '--chan_width', default=[None], nargs=1, type=float,
                        help='Optional: specify the channel width in native units of the original data (e.g. Hz or m/s). \n'
                              '(No default).')

    parser.add_argument('-x', '--suffix', default='png',
                        help='Optional: specify the output image file type: png, pdf, eps, jpeg, tiff, etc (default: %(default)s).')

    parser.add_argument('-m', '--imagemagick', nargs='?', type=str, default='', const='magick',
                        help='Optional: combine the main plots into single large image file using IMAGEMAGICK. If this option \n'
                             'is given with no argument we simply assume that image combination is executed by the "magick" \n'
                             'command. Otherwise, the argument of this option gives the full path to the CONVERT executable \n'
                             '(e.g. for older versions of imagemagick). Only the first multiwavelength image specified in \n'
                             '"surveys" argument is plotted next to the spectral line data.')
    
    parser.add_argument('-log', '--logfile-name', type=str, default=None,
                        help='Optional: Set name of the output log file. If not provided, log is named with catalog name and \n'
                             'date time stamp. Output will always have the ".log" suffix.  (No default).')

    parser.add_argument('-noid', '--no-source-id',
                        help='Turn off printing of the catalog source id number in the title of individual source plots.',
                        action='store_true')
    
    parser.add_argument('-cm', '--chan-maps',
                        help='Make a (multipage) pdf file for each source containing channel maps in 4x5 portrait layout. \n'
                             'Contours are plotted at +/- powers of 2 x RMS. No channel maps made if "-spec" flag is set.',
                        action='store_true')

    parser.add_argument('-spec', '--spec-only',
                        help='Skip making spatial images and only make spectra profiles.',
                        action='store_true')

    ###################################################################

    # Parse the arguments above
    args = parser.parse_args()
    catalog_file = args.catalog
    suffix = args.suffix
    original = args.original
    imagemagick = args.imagemagick

    if len(args.source_id) >= 1:
        if args.source_id[0] == '-1':
            args.source_id = [int(s) for s in args.source_id]
        elif '-' in args.source_id[0]:
            s_range = args.source_id[0].split('-')
            args.source_id = np.array(range(int(s_range[1]) - int(s_range[0]) + 1)) + int(s_range[0])
        else:
            args.source_id = [int(s) for s in args.source_id]

    try:
        beam = [float(b) for b in args.beam.split(',')]
    except:
        beam = []

    opt_view = args.image_size * u.arcmin

    # Set up logger
    log_basename = catalog_file.split('_cat')[0]
    auto_logname = True
    if args.logfile_name:
        log_basename = (args.logfile_name).split('.log')[0]
        auto_logname = False
    logger = Logger.get_logger(log_path=log_basename+'.log', auto_logname=auto_logname)#, clear_logs=False))

    print("\n")
    logger.info("*****************************************************************")
    logger.info("\tBeginning SoFiA-image-pipeline (SIP) {}.".format(version))

    if (len(args.surveys) == 0) and (not args.user_image):
        args.surveys = ['DSS2 Blue']
        logger.info("\tNo user specified image and no survey specified: will default to {}.".format(args.surveys[0]))
    elif ('none' in args.surveys):
        args.surveys = []
        if not args.user_image:
            logger.info("\tOffline mode requested: will not make ancillary data overlays.")
        else:
            logger.info("\tOffline mode requested: will only try user specified image.")
    surveys = tuple(args.surveys)

    if (suffix == 'eps') | (suffix == 'ps'):
        logger.warning("\t{} may have issues with transparency or making spectra.".format(suffix))

    if len(args.source_id):
        logger.info("\tWill only process selected sources: {}".format(args.source_id))

    # Read in the catalog file:
    if catalog_file.split(".")[-1] == "xml":
        logger.info("\tReading catalog in XML format.")
        logger.info("\tAlways assumes an xml file comes from SoFiA-2.")
        try:
            catalog = Table.read(catalog_file)
        except FileNotFoundError:
            logger.error("\t{} catalog not found. Typo or wrong directory?\n".format(catalog_file))
            exit()
        except:
            no_cat = True
    elif (catalog_file.split(".")[-1] == "txt") | (catalog_file.split(".")[-1] == "ascii"):
        logger.info("\tReading catalog in ascii format.")
        try:
            catalog = Table.read(catalog_file, format='ascii', header_start=18)  # Depends on SoFiA version!!! Do a brute force tes?
            logger.info("\tCatalog generated by SoFiA-2?")
            no_cat = False
        except FileNotFoundError:
            logger.error("\t{} catalog not found. Typo or wrong directory?\n".format(catalog_file))
            exit()
        except:
            no_cat = True
        if no_cat == True:
            try:
                catalog = Table.read(catalog_file, format='ascii', header_start=0)  # Depends on SoFiA version!!! Do a brute force tes?
                logger.info("\tCatalog generated by a source other than SoFiA-2?")
                no_cat = False
            except:
                no_cat = True
        if no_cat == True:
            logger.error("\tTrouble reading ascii format catalog.  A bug or generated by a different version of SoFiA?")
    else:
        logger.error("\tCatalog must be in xml or ascii format.\n")
        exit()

    # Check what's in the catalog; calculate ra, dec if necessary:
    if (('ra' not in catalog.colnames) and ('l' not in catalog.colnames)) and (not original):
        logger.error("\tLooks like catalog doesn't contain 'ra' and 'dec' or 'l' and 'b' columns. Re-run SoFiA with \n" \
            "\t\t'parameter.wcs = True' or you must include the path to the original fits file to derive \n" \
            "\t\tra, dec from the pixel values in the catalog.")
        logger.info("*****************************************************************\n")
        exit()
    elif ('ra' not in catalog.colnames) and (original):
        # Try to calculate an ra dec if not in the catalog. Deprecated?
        from src.modules.functions import get_radecfreq

        # This may be deprecated given that we insistence on running parameter.wcs, physical, offset = True
        logger.warning("\tLooks like catalog doesn't contain 'ra' and 'dec' columns. But can derive them with \n" \
            "\t\tthe pixel positions in the catalog provided.")
        logger.warning("\tThis probably means you ran SoFiA with 'parameter.wcs = False' which means the units \n" \
            "\t\t in your maps may be completely wacky! (Channel width knowledge is not maintained.)")
        ra, dec, freq = get_radecfreq(catalog, original)
        catalog['ra'] = ra
        catalog['dec'] = dec
        catalog['freq'] = freq
    # Add dummy columns to catalog for plotting source centers generalized to ICRS or Galactic later:
    catalog['pos_x'] = None
    catalog['pos_y'] = None

    # Set up some directories
    cubelet_dir = catalog_file.split("_cat.")[0] + '_cubelets/'
    if not os.path.isdir(cubelet_dir):
        logger.error("\tCubelet directory does not exist. Expecting {}. Need to run SoFiA-2 or restructure your"
            " directories.\n".format(cubelet_dir))
        exit()

    figure_dir = catalog_file.split("_cat.")[0] + '_figures/'
    if not os.path.isdir(figure_dir):
        logger.info("\tMaking figure directory.")
        os.system('mkdir {}'.format(figure_dir))

    logger.info("\tAssuming all requested sources are associated with {} line transition".format(args.spectral_line))    

    src_basename = cubelet_dir + catalog_file.split("/")[-1].split("_cat.")[0]

    # Make all the images on a source-by-source basis.  In future, could parallelize this.
    n_src = 0
    n_fail = 0
    failed_srcs = []

    # Do the hard work
    from src import make_images, make_spectra
    from src.modules import make_chan_maps
    from src.combine_images import combine_images

    from src.modules.functions import line_lookup
    spectral_line = line_lookup(args.spectral_line)

    for source in catalog:

        source['id'] = int(source['id'])  # For SoFiA-1 xml files--this doesn't work bc column type is float.

        if (not len(args.source_id)) or (source['id'] in args.source_id) or (args.source_id[0] == -1):
            logger.info(" ")
            logger.info("\t-Source {}: {}.".format(source['id'], source['name']))
            try:
                if args.spec_only:
                    logger.info("\tSkipping spatial images: only spectral profiles requested.")
                else:
                    x, p = make_images.main(source, src_basename, original, opt_view=opt_view, suffix=suffix, beam=beam,
                                    chan_width=args.chan_width[0], surveys=list(surveys), snr_range=args.snr_range,
                                    user_image=args.user_image, user_range=args.user_range, spec_line=spectral_line,
                                    noid=args.no_source_id)
                if args.chan_maps and not args.spec_only:
                    make_chan_maps.main(source, src_basename, suffix=suffix, beam=beam, noid=args.no_source_id, 
                                        opt_head=x, patch=p)
                make_spectra.main(source, src_basename, original, spec_line=spectral_line, suffix=suffix, 
                                  beam=beam, noid=args.no_source_id)
                n_src += 1
            except:
                failed_srcs.append(int(source['id']))
                n_fail += 1
                traceback.print_exc()
            try:
                if imagemagick:
                    code = ''.join(random.choices(string.ascii_letters + string.digits, k=6))
                    combine_images(source, src_basename, imagemagick, suffix=suffix, surveys=list(surveys), 
                                   user_image=args.user_image, code=code)
            except:
                os.system('rm -rf *_{1}.{0}'.format(suffix, code))
                pass

    if ((0 in args.source_id) or (-1 in args.source_id)) & (not args.spec_only):
        # Make a source object for the overview image
        from src.modules.functions import add_source

        logger.info(" ")
        logger.info("\tMaking summary images of full field.")
        catalog = add_source(catalog=catalog, fits_name=catalog_file.split("_cat.")[0] + '_mom0.fits')
        src_basename = src_basename.split('_cubelets')[0]
        make_images.main(catalog[-1], src_basename, original, opt_view=opt_view, suffix=suffix, beam=beam,
                    chan_width=args.chan_width[0], surveys=list(surveys), snr_range=args.snr_range,
                    user_image=args.user_image, user_range=args.user_range, spec_line=spectral_line,
                    catalog=catalog, noid=False)

    logger.info(" ")
    logger.info("\tDONE! Made images for {} sources.".format(n_src))
    if n_fail > 0:
        failed_msg = ' '.join([str(f) for f in failed_srcs])
        logger.warning(f"\Failed for {n_fail} sources with id number: {failed_msg}")
    logger.info("\tCreated log file: {}".format(Logger._log_path))
    logger.info("*****************************************************************\n")

if __name__ == '__main__':
    main()
