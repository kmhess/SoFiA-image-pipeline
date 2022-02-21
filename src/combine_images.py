import os


# Note, although this has been generalized, it seems to work best with png's!
def combine_images(source, src_basename, imgck, suffix='png'):

    # Configure expected file name:
    infile = src_basename.replace('cubelets', 'figures') + '_{}_'.format(source['id'])

    # Use terminal commands to assemble figures with imagemagick: https://imagemagick.org/index.php
    print("\tAssembling figures with imagemagick")
    new_file = "{}combo.{}".format(infile, suffix)
    os.system("{0} {1}mom0dss2blue.{2} {1}mom0hi.{2} {1}snr.{2} {1}mom1.{2}"
              " +append temp.{2}".format(imgck, infile, suffix))
    os.system("{0} {1}spec.{2} -resize 125% temp2.{2}".format(imgck, infile, suffix))
    os.system("{0} {1}specfull.{2} -resize 125% temp3.{2}".format(imgck, infile, suffix))
    os.system("{0} temp2.{2} temp3.{2} {1}pv.{2} +append temp4.{2}".format(imgck, infile, suffix))
    os.system("{0} temp.{2} temp4.{2} -append {1}".format(imgck, new_file, suffix))
    os.system('rm -rf temp*.{}'.format(suffix))

    return
