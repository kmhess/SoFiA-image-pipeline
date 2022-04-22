# Import Python libraries
from argparse import ArgumentParser, RawTextHelpFormatter
import os
from urllib.error import HTTPError

from astropy.coordinates import SkyCoord
from astropy import units as u

from modules.get_ancillary import *

###################################################################

parser = ArgumentParser(description="Save fits files from requested surveys to disk for later plotting (for example).",
                        formatter_class=RawTextHelpFormatter)

parser.add_argument('-ra', '--right_ascension', default=0.0, required=True,
                    help='Optional: specify the central RA of the image to retrieve in DEGREES. '
                         ' (default: %(default)s).')

parser.add_argument('-dec', '--declination', default=45.0, required=True,
                    help='Optional: specify the central DEC of the image to retrieve in DEGREES. '
                         ' (default: %(default)s).')

parser.add_argument('-i', '--image_size', default=[2.8, 2.3], nargs='+', type=float,#[2.8, 2.3],
                    help='Optional: specify the minimum survey image size to retrieve in DEGREES.'
                         ' (default: %(default)s).')

parser.add_argument('-s', '--surveys', default=['DSS2 Blue'], nargs='*', type=str, required=False,
                    help='Specify SkyView surveys to retrieve from astroquery on which to overlay HI contours.\n'
                         ' These additional non-SkyView options are also available: \'decals\',\'panstarrs\',\'hst\'.\n'
                         ' \'hst\' only refers to COSMOS HST (e.g. for CHILES).')

parser.add_argument('-o', '--outname', default=None, required=False,
                    help='Optional: specify the prefix for the resulting file. Will be appended by survey name.')

###################################################################

# Parse the arguments above
args = parser.parse_args()

surveys = set(args.surveys)
survey_string = ''
if args.outname:
    outname = args.outname + '_'
else:
    outname = ''

opt_view = args.image_size * u.deg
if len(opt_view) > 2:
    print("ERROR: -i image_size expects one or two arguments. Exiting.")
    exit()

# Can make this more flexible to include hmsdms entries.
hi_pos = SkyCoord(ra=args.right_ascension, dec=args.declination, unit='deg')

# If requested retrieve PanSTARRS false color imaging
if ('panstarrs' in surveys):
    if not os.path.isfile(outname + 'panstarrs.jpg'):
        pstar_im, pstar_head = get_panstarrs(hi_pos, opt_view=opt_view)
        if pstar_im:
            pstar_im.save(outname + 'panstarrs.jpg')
            pstar_fits_hdr = fits.PrimaryHDU(header=pstar_head)
            pstar_fits_hdr.writeto(outname + 'panstarrs_hdr.fits')
            survey_string += ' PanSTARRS'
        surveys.remove('panstarrs')
    else:
        print("\tERROR: {} already exists; will not overwrite. Choose a different outname prefix"
              " with `-o` flag. Continuing to next requested survey".format(outname + 'panstarrs.jpg'))

# If requested retrieve DECaLS false color imaging
if 'decals' in surveys:
    if not os.path.isfile(outname + 'decals.jpg'):
        decals_im, decals_head = get_decals(hi_pos, opt_view=opt_view)
        if decals_im:
            decals_im.save(outname + 'decals.jpg')
            decals_fits_hdr = fits.PrimaryHDU(header=decals_head)
            decals_fits_hdr.writeto(outname + 'decals_hdr.fits')
            survey_string += ' DECaLS'
        surveys.remove('decals')
    else:
        print("\tERROR: {} already exists; will not overwrite. Choose a different outname prefix"
              " with `-o` flag. Continuing to next requested survey".format(outname + 'decals.jpg'))

# If requested retrieve any number of survey images available through SkyView.
if len(surveys) > 0:
    for survey in surveys:
        if not os.path.isfile(outname + survey + '.fits'):
            try:
                overlay_image = get_skyview(hi_pos, opt_view=opt_view, survey=survey)
            except ValueError:
                print("\tERROR: \"{}\" may not among the survey hosted at skyview or survey names recognized by "
                      "astroquery. \n\t\tSee SkyView.list_surveys or SkyView.survey_dict from astroquery for valid "
                      "surveys.".format(survey))
            except HTTPError:
                print("\tERROR: http error 404 returned from SkyView query.  Skipping {}.".format(survey))
            overlay_image.writeto(outname + survey + '.fits')
            survey_string += ' {}'.format(survey)
        else:
            print("\tERROR: {} already exists; will not overwrite. Choose a different outname prefix"
                  " with `-o` flag. Continuing to next requested survey".format(outname + survey + '.fits'))

if len(survey_string) > 0:
    print("\n\tDONE! Saved survey images to disk for{}.".format(survey_string))
else:
    print("\n\tDONE! No survey images saved to disk.  Adjust input and try again?")
print("*****************************************************************\n")
