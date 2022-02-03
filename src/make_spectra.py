from datetime import datetime
import os

from astropy.io import ascii, fits
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np

from modules.functions import chan2freq
from modules.functions import get_info
from modules.functions import get_subcube


HI_restfreq = 1420405751.77 * u.Hz
optical_HI = u.doppler_optical(HI_restfreq)


###################################################################


# Save full spectrum to a txt file:
def get_noise_spec(source, src_basename, original=None):

    outfile = src_basename.replace('cubelets', 'figures') + '_{}_specfull.txt'.format(source['id'])

    if not os.path.isfile(outfile):

        if not original:
            print("\tWARNING: Original data cube not provided: making spectrum of subcube with noise.")
            fits_file = src_basename + '_{}_cube.fits'.format(source['id'])
            cube = fits.getdata(fits_file)
            mask = fits.getdata(src_basename + '_{}_mask.fits'.format (source['id']))
            channels = np.asarray(range(cube.shape[0])) + source['z_min']
        else:
            print("\tOriginal data cube provided: making full spectrum image with noise.")
            fits_file = original
            cube = get_subcube(source, original)
            mask = get_subcube(source, original[:-5] + '_mask.fits')
            channels = np.asarray(range(cube.shape[0]))

        mask2d = np.sum(mask, axis=0)
        frequency = chan2freq(channels, fits_file)
        spectrum = np.nansum(cube[:, mask2d != 0], axis=1)
        n_pix = 0 * channels + np.sum(mask2d != 0)

        with open('temp.txt', 'w') as f:
            f.write("# Integrated source spectrum with noise\n")
            f.write("# Creator: hessFDP_pipeline.py\n") #%s\n#\n" % sofia_version_full)
            f.write("# \n")
            f.write("# The source spectrum, with noise, is calculated by integrating over\n")
            f.write("# the 2D mask of the source in every channel.  This means every row \n")
            f.write("# has the same number of contributing pixels and the noise is the same\n")
            f.write("# for every point.\n")
            f.write("# \n")

        ascii.write([channels, frequency, spectrum, n_pix], 'temp2.txt', format='fixed_width_two_line',
                    names=['Channel', 'Frequency', 'Flux density', 'Pixels'])
        os.system("cat temp.txt temp2.txt > {}".format(outfile))
        os.system("rm temp.txt temp2.txt")


# Make full spectrum plot:
def make_specfull(source, src_basename, cube_params, suffix='png', full=False):

    outfile = src_basename.replace('cubelets', 'figures') + '_{}_specfull.{}'.format(source['id'], suffix)

    if not os.path.isfile(outfile):

        print("\tMaking HI spectrum plot covering the full frequency range.")
        spec = ascii.read(outfile[:-1*len(suffix)] + 'txt')  #, names=['Channel', 'Frequency', '"Flux density"', 'Pixels'])
        optical_velocity = (spec['Frequency'] * u.Hz).to(u.km / u.s, equivalencies=optical_HI)

        if full == True:
            fig = plt.figure(figsize=(15, 4))
        else:
            fig = plt.figure(figsize=(8, 4))

        ax_spec = fig.add_subplot(111)
        ax_spec.plot([optical_velocity[-1].value-10, optical_velocity[0].value+10], [0, 0], '--', color='gray')
        ax_spec.errorbar(optical_velocity[:].value, spec['Flux density'] / cube_params['pix_per_beam'], elinewidth=0.75,
                         yerr=source['rms'] * np.sqrt(spec['Pixels'] / cube_params['pix_per_beam']), capsize=1)
        ax_spec.set_title(source['name'])
        ax_spec.set_xlim(optical_velocity[-1].value-5, optical_velocity[0].value+5)
        ax_spec.set_ylabel("Integrated Flux [Jy]")
        ax_spec.set_xlabel("Optical Velocity [km/s]")

        spectrumJy = spec["Flux density"] / cube_params['pix_per_beam']
        if full == True:
            maskmin = chan2freq(source['z_min'], hdu=hdu_pb).to(u.km / u.s, equivalencies=optical_HI).value
            maskmax = chan2freq(source['z_max'], hdu=hdu_pb).to(u.km / u.s, equivalencies=optical_HI).value
            ax_spec.plot([maskmin, maskmin], [np.nanmin(spectrumJy), np.nanmax(spectrumJy)], ':', color='gray')
            ax_spec.plot([maskmax, maskmax], [np.nanmin(spectrumJy), np.nanmax(spectrumJy)], ':', color='gray')

        # Condition from Apertif experience that if the RFI is *really* bad, plot based on strength of HI profile
        if (np.max(spectrumJy) > 2.) | (np.min(spectrumJy) < -1.):
            ax_spec.set_ylim(np.max(spectrumJy[source['z_min']:source['z_max']+1]) * -2,
                             np.max(spectrumJy[source['z_min']:source['z_max']+1]) * 2)

        fig.savefig(outfile, bbox_inches='tight')

    return


# Make SoFiA masked spectrum plot (no noise):
def make_spec(source, src_basename, cube_params, suffix='png'):

    outfile = src_basename.replace('cubelets', 'figures') + '_{}_spec.{}'.format(source['id'], suffix)

    if not os.path.isfile(outfile):

        print("\tMaking HI SoFiA masked spectrum plot.")
        spec = ascii.read(src_basename + '_{}_spec.txt'.format(source['id']),
                          names=['Channel', 'Frequency', 'Flux density', 'Pixels'])
        optical_velocity = (spec['Frequency'] * u.Hz).to(u.km / u.s, equivalencies=optical_HI)

        fig = plt.figure(figsize=(8, 4))
        ax_spec = fig.add_subplot(111)
        ax_spec.plot([optical_velocity[-1].value-10, optical_velocity[0].value+10], [0, 0], '--', color='gray')
        ax_spec.errorbar(optical_velocity[:].value, spec['Flux density'] / cube_params['pix_per_beam'], elinewidth=0.75,
                         yerr=source['rms'] * np.sqrt(spec['Pixels'] / cube_params['pix_per_beam']), capsize=1)
        ax_spec.set_title(source['name'])
        ax_spec.set_xlim(optical_velocity[-1].value-5, optical_velocity[0].value+5)
        ax_spec.set_ylabel("Integrated Flux [Jy]")
        ax_spec.set_xlabel("Optical Velocity [km/s]")
        fig.savefig(outfile, bbox_inches='tight')

    return


def main(source, src_basename, original=None, suffix='png'):

    print("\n\tStart making spectral profiles of the spectral line source {}: {}.".format(source['id'], source['name']))

    # Get beam information from the source cubelet
    cube_params = get_info(src_basename + '_{}_cube.fits'.format(source['id']))

    # Make plot of SoFiA masked spectrum
    make_spec(source, src_basename, cube_params, suffix=suffix)

    # If desired,
    outfile = src_basename.replace('cubelets', 'figures') + '_{}_specfull.txt'.format(source['id'])
    if original or (not os.path.isfile(outfile)):
        get_noise_spec(source, src_basename, original)

    # Make full spectrum plot
    make_specfull(source, src_basename, cube_params, suffix=suffix, full=False)

    plt.close('all')

    print("\tDone making spectral profiles of the spectral line source {}: {}.".format(source['id'], source['name']))

    return


if __name__ == '__main__':

    main(source, src_basename, original=None, suffix='png')
