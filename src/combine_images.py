import os


# Note, although this has been generalized, it seems to work best with png's!

def combine_images(source, src_basename, imgck, suffix='png', surveys='DSS2 Blue', user_image=None):
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
    convert_im = "convert"

    # Configure expected file name:
    infile = src_basename.replace('cubelets', 'figures') + '_{}_'.format(source['id'])

    # Use terminal commands to assemble figures with imagemagick: https://imagemagick.org/index.php
    print("\tAssembling figures with imagemagick")
    new_file = "{}combo.{}".format(infile, suffix)
    if user_image and os.path.exists('{0}mom0_usr.{1}'.format(infile, suffix)):
        os.system("{0} {1}mom0_usr.{2} {1}mom0hi.{2} {1}snr.{2} {1}mom1.{2}"
                  " +append temp.{2}".format(imgck, infile, suffix))
    else:
        os.system("{0} {1}mom0{3}.{2} {1}mom0hi.{2} {1}snr.{2} {1}mom1.{2}"
                  " +append temp.{2}".format(imgck, infile, suffix, surveys[0].replace(" ", "").lower()))
    os.system("{0} {1}spec.{2} -resize 125% temp2.{2}".format(imgck, infile, suffix))
    os.system("{0} {1}specfull.{2} -resize 125% temp3.{2}".format(imgck, infile, suffix))
    os.system("{0} temp2.{2} temp3.{2} {1}pv.{2} +append temp4.{2}".format(imgck, infile, suffix))
    os.system("{0} temp.{2} temp4.{2} -append {1}".format(imgck, new_file, suffix))
    os.system('rm -rf temp*.{}'.format(suffix))

    return
