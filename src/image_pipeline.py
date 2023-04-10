#!/usr/bin/env python3

# Import default Python libraries
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import sys

# Import installed Python libraries
from astropy.io import fits  # change to table.read?
from astropy.table import Table
from astropy import units as u
import pkg_resources  # part of setuptools
import numpy as np

from src import make_images, make_spectra
from src.modules.functions import get_radecfreq
from src.combine_images import combine_images

version = pkg_resources.require("SoFiA-image-pipeline")[0].version

###################################################################

def main():
    parser = ArgumentParser(description="Create images from a SoFiA catalog and cubelets, or fits file. \n"
                                        "Only works with SoFiA-2 and wcs=True (for now).",
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('-c', '--catalog', required=True,
                        help='Required: Specify the input XML or ascii catalog name. No default.')

    parser.add_argument('-id', '--source-id', default=[], nargs='*', type=int,
                        help='Space-separated list of sources to include in the plotting. Default all sources')

    parser.add_argument('-x', '--suffix', default='png',
                        help='Optional: specify the output image file type: png, pdf, eps, jpeg, tiff, etc (default: %(default)s).')

    parser.add_argument('-o', '--original', default=None,
                        help='Optional: specify the original fits data: used for plotting HI spectra *with* noise over \n'
                            ' the full frequency range of the cube. Otherwise, plot with noise over frequency range\n'
                            ' in the cubelet.  Uses 2D mask to integrate. (No default).')

    parser.add_argument('-b', '--beam', default=None,
                        help='Optional: specify the beam dimensions (bmaj,bmin,bpa) in arcsec, arcsec, deg. If only 1 value\n'
                            ' is given, assume a circular beam. If 2 values are given, assume PA = 0. (No default).')

    parser.add_argument('-cw', '--chan_width', default=[None], nargs=1, type=float,
                        help='Optional: specify the channel width in native units of the original data (e.g. Hz or m/s).'
                             ' (No default).')

    parser.add_argument('-i', '--image_size', default=[6], nargs=1, type=float,
                        help='Optional: specify the minimum survey image size to retrieve in arcmin.  It will be adjusted if\n'
                            ' the HI mask is larger. Note max panstarrs image size is 8 arcmin (default: %(default)s).')

    parser.add_argument('-snr', '--snr-range', default=[2., 3.], nargs=2, type=float,
                        help='Optional: specify which SNRmin and SNRmax values should be used to set the lowest reliable HI \n'
                            ' contour in the figures. The contour level is calculated as the median value in the mom0 image\n'
                            ' of all pixels whose SNR value is within the given range. Default is [2,3].')

    parser.add_argument('-s', '--surveys', default=[], nargs='*', type=str,
                        help='Specify SkyView surveys to retrieve from astroquery on which to overlay HI contours.\n'
                            ' These additional non-SkyView options are also available: \'decals\',\'panstarrs\',\'hst\'.\n'
                            ' \'hst\' only refers to COSMOS HST (e.g. for CHILES). Default is "DSS2 Blue" if no user\n' 
                            ' provided image.')

    parser.add_argument('-m', '--imagemagick', nargs='?', type=str, default='', const='convert',
                        help='Optional: combine the main plots into single large image file using the IMAGEMAGICK CONVERT task.\n'
                            ' If this option is given with no argument we simply assume that CONVERT is executed by the "convert"\n'
                            ' command. Otherwise, the argument of this option gives the full path to the CONVERT executable.\n'
                            ' Only the first multiwavelength image specified in "surveys" argument is plotted next to the\n'
                            ' spectral line data.')

    parser.add_argument('-ui', '--user-image', default=None,
                        help='Optional: Full path to the FITS image on which to overlay HI contours.')

    parser.add_argument('-ur', '--user-range', default=[10., 99.], nargs=2, type=float,
                        help='Optional: Percentile range used when displaying the user image (see "-ui"). Default is [10,99].')

    ###################################################################

    # Parse the arguments above
    args = parser.parse_args()
    suffix = args.suffix
    original = args.original
    imagemagick = args.imagemagick

    try:
        beam = [float(b) for b in args.beam.split(',')]
    except:
        beam = []
    opt_view = args.image_size * u.arcmin

    print("\n*****************************************************************")
    print("\tBeginning SoFiA-image-pipeline (SIP) {}.".format(version))

    if (len(args.surveys) == 0) and (not args.user_image):
        args.surveys = ['DSS2 Blue']
        print("\tNo user specified image and no survey specified: will default to {}.".format(args.surveys[0]))
    surveys = tuple(args.surveys)

    if (suffix == 'eps') | (suffix == 'ps'):
        print("\tWARNING: {} may have issues with transparency or making spectra.".format(suffix))

    if len(args.source_id):
        print("\tWill only process selected sources: {}".format(args.source_id))

    # Read in the catalog file:
    catalog_file = args.catalog

    if catalog_file.split(".")[-1] == "xml":
        print("\tReading catalog in XML format.")
        print("\tWARNING: always assumes an xml file comes from SoFiA-2.")
        try:
            catalog = Table.read(catalog_file)
            sofia = 2
        except FileNotFoundError:
            print("\tERROR: {} catalog not found. Typo or wrong directory?\n".format(catalog_file))
            exit()
        except:
            no_cat = True
    elif (catalog_file.split(".")[-1] == "txt") | (catalog_file.split(".")[-1] == "ascii"):
        print("\tReading catalog in ascii format.")
        try:
            catalog = Table.read(catalog_file, format='ascii', header_start=18)  # Depends on SoFiA version!!! Do a brute force tes?
            print("\tCatalog generated by SoFiA-2?")
            sofia = 2
            no_cat = False
        except FileNotFoundError:
            print("\tERROR: {} catalog not found. Typo or wrong directory?\n".format(catalog_file))
            exit()
        except:
            no_cat = True
        if no_cat == True:
            try:
                catalog = Table.read(catalog_file, format='ascii', header_start=1)  # Depends on SoFiA version!!! Do a brute force tes?
                print("\tCatalog generated by SoFiA-1 or another source?")
                sofia = 1
                no_cat = False
            except:
                no_cat = True
        if no_cat == True:
            print("\tERROR: Trouble reading ascii format catalog.  A bug or generated by a different version of SoFiA?")
    else:
        print("\tERROR: Catalog must be in xml or ascii format.\n")
        exit()

    # Check what's in the catalog; calculate ra, dec if necessary:
    if (('ra' not in catalog.colnames) and ('l' not in catalog.colnames)) and (not original):
        print("\tERROR: Looks like catalog doesn't contain 'ra' and 'dec' or 'l' and 'b' columns. Re-run SoFiA with \n" \
            "\t\t'parameter.wcs = True' or you must include the path to the original fits file to derive \n" \
            "\t\tra, dec from the pixel values in the catalog.")
        print("*****************************************************************\n")
        exit()
    elif ('ra' not in catalog.colnames) and (original):
        # This may be deprecated given that we insistence on running parameter.wcs, physical, offset = True
        print("\tWARNING: Looks like catalog doesn't contain 'ra' and 'dec' columns. But can derive them with \n" \
            "\t\tthe pixel positions in the catalog provided.")
        print("\tWARNING: This probably means you ran SoFiA with 'parameter.wcs = False' which means the units \n" \
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
        print("\tERROR: Cubelet directory does not exist. Expecting {}. Need to run SoFiA-2 or restructure your"
            " directories.\n".format(cubelet_dir))
        exit()

    figure_dir = catalog_file.split("_cat.")[0] + '_figures/'
    if not os.path.isdir(figure_dir):
        print("\tMaking figure directory.")
        os.system('mkdir {}'.format(figure_dir))

    src_basename = cubelet_dir + catalog_file.split("/")[-1].split("_cat.")[0]

    # Make all the images on a source-by-source basis.  In future, could parallelize this.
    n_src = 0

    for source in catalog:

        source['id'] = int(source['id'])  # For SoFiA-1 xml files--this doesn't work bc column type is float.

        if not len(args.source_id) or source['id'] in args.source_id:
            print("\n\t-Source {}: {}.".format(source['id'], source['name']))
            make_images.main(source, src_basename, opt_view=opt_view, suffix=suffix, sofia=sofia, beam=beam,
                             chan_width=args.chan_width[0], surveys=list(surveys), snr_range=args.snr_range,
                             user_image=args.user_image, user_range=args.user_range)
            make_spectra.main(source, src_basename, original, suffix=suffix, beam=beam)

            if imagemagick:
                combine_images(source, src_basename, imagemagick, suffix=suffix, surveys=list(surveys), user_image=args.user_image)

            n_src += 1

    print("\n\tDONE! Made images for {} sources.".format(n_src))
    print("*****************************************************************\n")


if __name__ == '__main__':
    main()
