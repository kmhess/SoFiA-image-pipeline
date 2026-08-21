import os
import subprocess

from sip.modules.logger import Logger

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
    for im in ['mom0', 'snr', 'mom1', 'mom2', 'specboth']:
        subprocess.run([imgck, '{0}{1}.{2}'.format(infile, im, suffix), '-gravity', 'west', '-chop', '40x0',
                        '{0}_{1}.{2}'.format(im, code, suffix)])
        
    # Use imagemagick to append images together:
    if user_image and os.path.exists('{0}mom0_usr.{1}'.format(infile, suffix)):
        subprocess.run([imgck, '{}mom0_usr.{}'.format(infile, suffix), 'mom0_{}.{}'.format(code, suffix), 'snr_{}.{}'.format(code, suffix),
                        'mom1_{}.{}'.format(code, suffix), 'mom2_{}.{}'.format(code, suffix), '+append', '-gravity', 'south', '-splice',
                        '0x18', 'temp_{}.{}'.format(code, suffix)])
    elif surveys and os.path.exists('{0}mom0_{2}.{1}'.format(infile, suffix, 
                                        surveys[0].replace(" ", "").lower().replace('decals-dr9', 'decals'))):
        subprocess.run([imgck, '{}mom0_{}.{}'.format(infile, surveys[0].replace(" ", "").lower().replace('decals-dr9', 'decals'),suffix),
                        'mom0_{}.{}'.format(code, suffix), 'snr_{}.{}'.format(code, suffix), 'mom1_{}.{}'.format(code, suffix),
                        'mom2_{}.{}'.format(code, suffix), '+append', '-gravity', 'south', '-splice', '0x18',
                        'temp_{}.{}'.format(code, suffix)])
    else:
        logger.warning("\tNo ancillary data image available for source {}.".format(source['id']))
        subprocess.run([imgck, '{}mom0.{}'.format(infile, suffix), 'snr_{}.{}'.format(code, suffix), 'mom1_{}.{}'.format(code, suffix), 
                        'mom2_{}.{}'.format(code, suffix), '+append', '-gravity', 'south', '-splice', '0x18', 
                        'temp_{}.{}'.format(code, suffix)])
    subprocess.run([imgck, '{}spec.{}'.format(infile, suffix), '-resize', '133%', 'temp2_{}.{}'.format(code, suffix)])
    subprocess.run([imgck, 'specboth_{}.{}'.format(code, suffix), '-resize', '133%', 'temp3_{}.{}'.format(code, suffix)])

    # Remove redundant y-axes for pv plots and create a little space between pv and spectra:
    if os.path.isfile('{0}pv.{1}'.format(infile, suffix)):
        if os.path.isfile('{0}pv_min.{1}'.format(infile, suffix)):
            if 'freq' in source.colnames:
                subprocess.run([imgck, '{}pv_min.{}'.format(infile, suffix), '-gravity', 'west', '-chop', '132x0', 
                                'pv_min_{}.{}'.format(code, suffix)])
                subprocess.run([imgck, '{}pv.{}'.format(infile, suffix), '-gravity', 'east', '-chop', '128x0', '-splice', '40x0', 
                                'pv_{}.{}'.format(code, suffix)])
                subprocess.run([imgck, 'temp2_{}.{}'.format(code, suffix), 'temp3_{}.{}'.format(code, suffix),
                                'pv_{}.{}'.format(code, suffix), '-gravity', 'west', '-splice', '20x0', 
                                'pv_min_{}.{}'.format(code, suffix), '+append', 'temp4_{}.{}'.format(code, suffix)])
            else:
                subprocess.run([imgck, '{}pv_min.{}'.format(infile, suffix), '-gravity', 'west', '-chop', '40x0', 'pv_min_{}.{}'.format(code, suffix)])
                subprocess.run([imgck, 'temp2_{}.{}'.format(code, suffix), 'temp3_{}.{}'.format(code, suffix), '{}pv.{}'.format(infile, suffix), 
                                '-gravity', 'west', '-splice', '20x0', 'pv_min_{}.{}'.format(code, suffix), '+append', 'temp4_{}.{}'.format(code, suffix)])
        else:
            subprocess.run([imgck, 'temp2_{}.{}'.format(code, suffix), 'temp3_{}.{}'.format(code, suffix), '{}pv.{}'.format(infile, suffix), 
                            '-gravity', 'west', '-splice', '20x0', '+append', 'temp4_{}.{}'.format(code, suffix)])
        subprocess.run([imgck, 'temp_{}.{}'.format(code, suffix), 'temp4_{}.{}'.format(code, suffix), '-append', new_file])
    else:
        subprocess.run([imgck, 'temp2_{}.{}'.format(code, suffix), 'temp3_{}.{}'.format(code, suffix), '-gravity', 'west', '-splice', '20x0',
                        '+append', 'temp4_{}.{}'.format(code, suffix)])
        subprocess.run([imgck, 'temp_{}.{}'.format(code, suffix), 'temp4_{}.{}'.format(code, suffix), '-append', new_file])

    new_file_size = os.path.getsize(new_file)

    if new_file_size > file_size_limit:
        logger.info('\tReducing size of combined image to {0:.0f}% of original (it was {1:.1e}B)'.format(100*file_size_limit/new_file_size, 
                                                                                                         new_file_size))
        subprocess.run([imgck, new_file, '-resize', '{0:.0f}%'.format(100*file_size_limit/new_file_size), new_file])
    subprocess.run(['rm -rf *_{}.{}'.format(code, suffix)], shell=True)

    return
