import os
import random
import string


# Note, although this has been generalized, it seems to work best with png's!

def combine_images(source, src_basename, imgck, suffix='png', surveys='DSS2 Blue', user_image=None, file_size_limit = 8e+5):
    """_summary_

    :param source: source for which to combine images
    :type source: Astropy source object?
    :param src_basename: basename for the source image files
    :type src_basename: str
    :param suffix: filetype, defaults to 'png'
    :type suffix: str, optional
    """

    # Specify the command to use imagemagick's convert (karma has a convert which may conflict)
    # convert_im = "/usr/local/Cellar/imagemagick/7.1.0-13/bin/convert"

    # Configure expected file name:
    infile = src_basename.replace('cubelets', 'figures') + '_{}_'.format(source['id'])
    code = ''.join(random.choices(string.ascii_letters + string.digits, k=6))

    # Use terminal commands to assemble figures with imagemagick: https://imagemagick.org/index.php
    print("\tAssembling figures with imagemagick")
    new_file = "{}combo.{}".format(infile, suffix)
    if user_image and os.path.exists('{0}mom0_usr.{1}'.format(infile, suffix)):
        os.system("{0} {1}mom0_usr.{2} {1}mom0.{2} {1}snr.{2} {1}mom1.{2} {1}mom2.{2}"
                  " +append temp{3}.{2}".format(imgck, infile, suffix, code))
    elif surveys:
        os.system("{0} {1}mom0_{3}.{2} {1}mom0.{2} {1}snr.{2} {1}mom1.{2} {1}mom2.{2} +append"
                  " temp{4}.{2}".format(imgck, infile, suffix,
                                        surveys[0].replace(" ", "").lower().replace('decals-dev', 'decals'), code))
    else:
        print("\tWARNING: No ancillary data image available for source {}.".format(source['id']))
        os.system("{0} {1}mom0.{2} {1}snr.{2} {1}mom1.{2} {1}mom2.{2} +append temp{3}.{2}".format(imgck, infile, suffix, code))
    os.system("{0} {1}spec.{2} -resize 125% temp2{3}.{2}".format(imgck, infile, suffix, code))
    os.system("{0} {1}specfull.{2} -resize 125% temp3{3}.{2}".format(imgck, infile, suffix, code))
    os.system("{0} temp2{3}.{2} temp3{3}.{2} {1}pv.{2} {1}pv_min.{2} +append temp4{3}.{2}".format(imgck, infile,
                                                                                                  suffix, code))
    os.system("{0} temp{3}.{2} temp4{3}.{2} -append {1}".format(imgck, new_file, suffix, code))
    new_file_size = os.path.getsize(new_file)
    if new_file_size > file_size_limit:
        print('\tReducing size of combined image to {0:.0f}% of original (it was {1:.1e}B)'.format(100*file_size_limit/new_file_size, new_file_size))
        os.system("{0} {1} -resize {2:.0f}% {1}".format(imgck, new_file, 100*file_size_limit/new_file_size))
    os.system('rm -rf temp*{1}.{0}'.format(suffix, code))

    return
