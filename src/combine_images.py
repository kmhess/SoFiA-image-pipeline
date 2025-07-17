import os

from src.modules.logger import Logger

logger = Logger.get_logger()


# Note, although this has been generalized, it seems to work best with png's!

def combine_images(source, src_basename, imgck, suffix='png', surveys='DSS2 Blue', user_image=None, file_size_limit=8e+5,
                   code=None):
    """_summary_

    :param source: source for which to combine images
    :type source: Astropy source object?
    :param src_basename: basename for the source image files
    :type src_basename: str
    :param user_image: path to the FITS image on which HI contours were overlaid
    :type user_image: str
    :param suffix: filetype, defaults to 'png'
    :type suffix: str, optional
    :param code: unique random code to identify temp files
    :type code: str
    """

    # Specify the command to use imagemagick's convert (karma has a convert which may conflict)
    # convert_im = "/usr/local/Cellar/imagemagick/7.1.0-13/bin/convert"

    # Configure expected file name:
    infile = src_basename.replace('cubelets', 'figures') + '_{}_'.format(source['id'])

    # Use terminal commands to assemble figures with imagemagick: https://imagemagick.org/index.php
    logger.info("\tAssembling figures with imagemagick")
    new_file = "{}combo.{}".format(infile, suffix)
    # Remove redundant y-axies for the 2D images:
    for im in ['mom0', 'snr', 'mom1', 'mom2', 'specfull']:
        os.system('{0} {1}{2}.{3} -gravity west -chop 40x0 {2}_{4}.{3}'.format(imgck, infile, im, suffix, code))

    # Use imagemagick to append images together:
    if user_image and os.path.exists('{0}mom0_usr.{1}'.format(infile, suffix)):
        os.system("{0} {1}mom0_usr.{2} mom0_{3}.{2} snr_{3}.{2} mom1_{3}.{2} mom2_{3}.{2} +append"
                  " -gravity south -splice 0x18 temp_{3}.{2}".format(imgck, infile, suffix, code))
    elif surveys and os.path.exists('{0}mom0_{2}.{1}'.format(infile, suffix, 
                                        surveys[0].replace(" ", "").lower().replace('decals-dr9', 'decals'))):
        os.system("{0} {1}mom0_{3}.{2} mom0_{4}.{2} snr_{4}.{2} mom1_{4}.{2} mom2_{4}.{2} +append"
                  " -gravity south -splice 0x18 temp_{4}.{2}".format(imgck, infile, suffix,
                                        surveys[0].replace(" ", "").lower().replace('decals-dr9', 'decals'), code))
    else:
        logger.warning("\tNo ancillary data image available for source {}.".format(source['id']))
        os.system("{0} {1}mom0.{2} snr_{3}.{2} mom1_{3}.{2} mom2_{3}.{2} +append"
                  " -gravity south -splice 0x18 temp_{3}.{2}".format(imgck, infile, suffix, code))
    os.system("{0} {1}spec.{2} -resize 133% temp2_{3}.{2}".format(imgck, infile, suffix, code))
    os.system("{0} specfull_{3}.{2} -resize 133% temp3_{3}.{2}".format(imgck, infile, suffix, code))

    # Remove redundant y-axes for pv plots and create a little space between pv and spectra:
    if os.path.isfile('{0}pv_min.{1}'.format(infile, suffix)):
        if 'freq' in source.colnames:
            os.system('{0} {1}pv_min.{2} -gravity west -chop 132x0 pv_min_{3}.{2}'.format(imgck, infile, suffix, code))
            os.system('{0} {1}pv.{2} -gravity east -chop 128x0 -splice 40x0 pv_{3}.{2}'.format(imgck, infile, 
                                                                                               suffix, code))
            os.system("{0} temp2_{3}.{2} temp3_{3}.{2} pv_{3}.{2} -gravity west -splice 20x0 pv_min_{3}.{2}"
                      " +append temp4_{3}.{2}".format(imgck, infile, suffix, code))
        else:
            os.system('{0} {1}pv_min.{2} -gravity west -chop 40x0 pv_min_{3}.{2}'.format(imgck, infile, suffix, code))
            os.system("{0} temp2_{3}.{2} temp3_{3}.{2} {1}pv.{2} -gravity west -splice 20x0 pv_min_{3}.{2}"
                      " +append temp4_{3}.{2}".format(imgck, infile, suffix, code))
    else:
        os.system("{0} temp2_{3}.{2} temp3_{3}.{2} {1}pv.{2} -gravity west -splice 20x0"
                  " +append temp4_{3}.{2}".format(imgck, infile, suffix, code))
    os.system("{0} temp_{3}.{2} temp4_{3}.{2} -append {1}".format(imgck, new_file, suffix, code))

    new_file_size = os.path.getsize(new_file)

    if new_file_size > file_size_limit:
        logger.info('\tReducing size of combined image to {0:.0f}% of original (it was {1:.1e}B)'.format(100*file_size_limit/new_file_size, 
                                                                                                   new_file_size))
        os.system("{0} {1} -resize {2:.0f}% {1}".format(imgck, new_file, 100*file_size_limit/new_file_size))
    os.system('rm -rf *_{1}.{0}'.format(suffix, code))

    return
