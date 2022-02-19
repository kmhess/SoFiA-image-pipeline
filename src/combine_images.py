import os


def combine_images(source, src_basename):

    # Specify the command to use imagemagick's convert (karma has a convert which may conflict)
    # convert_im = "/usr/local/Cellar/imagemagick/7.1.0-13/bin/convert"
    convert_im = "convert"

    # Configure expected file name:
    infile = src_basename.replace('cubelets', 'figures') + '_{}_'.format(source['id'])

    # Use terminal commands to assemble figures with imagemagick: https://imagemagick.org/index.php
    print("\n\tAssembling figures with imagemagck for source {}: {}.".format(source['id'], source['name']))
    new_png = "{}combo.png".format(infile)
    os.system(
        "{} {}mom0dss2blue.png {}mom0hi.png {}snr.png {}mom1.png +append temp.png".format(convert_im, infile, infile,
                                                                                          infile, infile, infile))
    os.system("{} {}spec.png -resize 125% temp2.png".format(convert_im, infile))
    os.system("{} {}specfull.png -resize 125% temp3.png".format(convert_im, infile))
    os.system("{} temp2.png temp3.png {}pv.png +append temp4.png".format(convert_im, infile, infile, infile))
    os.system("{} temp.png temp4.png -append {}".format(convert_im, new_png))
    os.system('rm -rf temp.png temp2.png temp3.png temp4.png')

    return
