# from datetime import datetime
import os

from astropy.io import ascii, fits
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np

from modules.functions import chan2freq, chan2vel, get_info, get_subcube, felo2vel


HI_restfreq = 1420405751.77 * u.Hz
optical_HI = u.doppler_optical(HI_restfreq)


###################################################################


# Save full spectrum to a txt file:
def get_noise_spec(source, src_basename, cube_params, original=None):

    outfile = src_basename.replace('cubelets', 'figures') + '_{}_specfull.txt'.format(source['id'])

    if not os.path.isfile(outfile):

        if not original:
            print("\tWARNING: Original data cube not provided: making spectrum of subcube with noise.")
            fits_file = src_basename + '_{}_cube.fits'.format(source['id'])
            cube = fits.getdata(fits_file)
            mask = fits.getdata(src_basename + '_{}_mask.fits'.format(source['id']))
            spec_template = ascii.read(src_basename + '_{}_spec.txt'.format(source['id']),
                                       names=['chan', 'col2', 'f_sum', 'n_pix'])
            channels = spec_template['chan']
        else:
            print("\tOriginal data cube provided: making full spectrum image with noise.")
            fits_file = original
            cube = get_subcube(source, original)
            mask = get_subcube(source, original[:-5] + '_mask.fits')
            spec_template = None
            channels = np.asarray(range(cube.shape[0]))

        mask2d = np.sum(mask, axis=0)
        spectrum = np.nansum(cube[:, mask2d != 0], axis=1)
        n_pix = 0 * channels + np.sum(mask2d != 0)

        with open('temp.txt', 'w') as f:
            f.write("# Integrated source spectrum with noise\n")
            f.write("# Creator: SoFiA-image-pipeline.py\n") #%s\n#\n" % sofia_version_full)
            f.write("# \n")
            f.write("# The source spectrum, with noise, is calculated by integrating over\n")
            f.write("# the 2D mask of the source in every channel.  This means every row \n")
            f.write("# has the same number of contributing pixels and the noise is the same\n")
            f.write("# for every point.\n")
            f.write("# \n")

            if 'freq' in source.colnames:
                frequency = spec_template['col2'] if spec_template else chan2freq(channels, fits_file)
                ascii.write([channels, frequency, spectrum, n_pix], 'temp2.txt', format='fixed_width_two_line',
                            names=['chan', 'freq', 'f_sum', 'n_pix'])
            elif 'FELO' in cube_params['spec_axis']:
                velocities = spec_template['col2'] if spec_template else felo2vel(channels, fits_file)
                ascii.write([channels, velocities, spectrum, n_pix], 'temp2.txt', format='fixed_width_two_line',
                            names=['chan', 'velo', 'f_sum', 'n_pix'])
            else:
                velocities = spec_template['col2'] if spec_template else chan2vel(channels, fits_file)
                ascii.write([channels, velocities, spectrum, n_pix], 'temp2.txt', format='fixed_width_two_line',
                            names=['chan', 'velo', 'f_sum', 'n_pix'])
        os.system("cat temp.txt temp2.txt > {}".format(outfile))
        os.system("rm temp.txt temp2.txt")

        return


# Make full spectrum plot:
def make_specfull(source, src_basename, cube_params, suffix='png', full=False):

    outfile = src_basename.replace('cubelets', 'figures') + '_{}_specfull.{}'.format(source['id'], suffix)

    if not os.path.isfile(outfile):

        print("\tMaking HI spectrum plot covering the full frequency range.")
        convention = 'Optical'
        if 'freq' in source.colnames:
            spec = ascii.read(outfile[:-1*len(suffix)] + 'txt')
            optical_velocity = (spec['freq'] * u.Hz).to(u.km / u.s, equivalencies=optical_HI).value
            maskmin = (spec['freq'][spec['chan'] == source['z_min']] * u.Hz).to(u.km / u.s, equivalencies=optical_HI).value
            maskmax = (spec['freq'][spec['chan'] == source['z_max']] * u.Hz).to(u.km / u.s, equivalencies=optical_HI).value
        else:
            if 'vrad' in source.colnames: convention = 'Radio'
            spec = ascii.read(outfile[:-1 * len(suffix)] + 'txt', names=['chan', 'velo', 'f_sum', 'n_pix'])
            optical_velocity = (spec['velo'] * u.m / u.s).to(u.km / u.s).value
            maskmin = (spec['velo'][spec['chan'] == source['z_min']] * u.m / u.s).to(u.km / u.s).value
            maskmax = (spec['velo'][spec['chan'] == source['z_max']] * u.m / u.s).to(u.km / u.s).value

        if full == True:
            fig = plt.figure(figsize=(15, 4))
        else:
            fig = plt.figure(figsize=(8, 4))

        ax_spec = fig.add_subplot(111)
        ax_spec.plot([np.min(optical_velocity) - 10, np.max(optical_velocity) + 10], [0, 0], '--', color='gray')
        ax_spec.errorbar(optical_velocity, spec['f_sum'] / cube_params['pix_per_beam'], elinewidth=0.75,
                         yerr=source['rms'] * np.sqrt(spec['n_pix'] / cube_params['pix_per_beam']), capsize=1)
        ax_spec.set_title(source['name'])
        ax_spec.set_xlim(np.min(optical_velocity) - 5, np.max(optical_velocity) + 5)
        ax_spec.set_ylabel("Integrated Flux [Jy]")
        ax_spec.set_xlabel("{} {} Velocity [km/s]".format(cube_params['spec_sys'].capitalize(), convention))

        spectrumJy = spec["f_sum"] / cube_params['pix_per_beam']

        # Plot limit of SoFiA mask
        ymin, ymax = ax_spec.get_ylim()
        ax_spec.plot([maskmin, maskmin], [0.95*ymin, 0.95*ymax], ':', color='gray')
        ax_spec.plot([maskmax, maskmax], [0.95*ymin, 0.95*ymax], ':', color='gray')

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
        convention = 'Optical'
        if 'freq' in source.colnames:
            spec = ascii.read(src_basename + '_{}_spec.txt'.format(source['id']),
                              names=['chan', 'freq', 'f_sum', 'n_pix'])
            optical_velocity = (spec['freq'] * u.Hz).to(u.km / u.s, equivalencies=optical_HI).value
        else:
            if 'vrad' in source.colnames: convention = 'Radio'
            spec = ascii.read(src_basename + '_{}_spec.txt'.format(source['id']),
                              names=['chan', 'velo', 'f_sum', 'n_pix'])
            optical_velocity = (spec['velo'] * u.m / u.s).to(u.km / u.s).value

        # Get spec units (Jy or Jy/beam) to check whether the division by the beam area has been made
        ll = 0
        while not ('f_sum' in spec.meta['comments'][ll] and 'chan' in spec.meta['comments'][ll] and
                   'n_pix' in spec.meta['comments'][ll]):
            ll += 1
        specunits = (spec.meta['comments'][ll+1].split()[spec.meta['comments'][ll].split().index('f_sum')])

        fig = plt.figure(figsize=(8, 4))
        ax_spec = fig.add_subplot(111)
        ax_spec.plot([np.min(optical_velocity) - 10, np.max(optical_velocity) + 10], [0, 0], '--', color='gray')
        if specunits == 'Jy/beam':
            ax_spec.errorbar(optical_velocity, spec['f_sum'] / cube_params['pix_per_beam'], elinewidth=0.75,
                             yerr=source['rms'] * np.sqrt(spec['n_pix'] / cube_params['pix_per_beam']), capsize=1)
        elif specunits == 'Jy':
            ax_spec.errorbar(optical_velocity, spec['f_sum'], elinewidth=0.75,
                             yerr=source['rms'] * np.sqrt(spec['n_pix'] / cube_params['pix_per_beam']), capsize=1)
        ax_spec.set_title(source['name'])
        ax_spec.set_xlim(np.min(optical_velocity) - 5, np.max(optical_velocity) + 5)
        ax_spec.set_ylabel("Integrated Flux [Jy]")
        ax_spec.set_xlabel("{} {} Velocity [km/s]".format(cube_params['spec_sys'].capitalize(), convention))
        fig.savefig(outfile, bbox_inches='tight')

    return


def main(source, src_basename, original=None, suffix='png', beam=None):

    print("\tStart making spectral profiles")

    # Get beam information from the source cubelet
    cube_params = get_info(src_basename + '_{}_cube.fits'.format(source['id']), beam)

    # Make plot of SoFiA masked spectrum
    make_spec(source, src_basename, cube_params, suffix=suffix)

    # Make text file of spectrum with noise; use full frequency range of original cube if provided:
    # Can be a bit more precise here in the output options/specification.
    outfile = src_basename.replace('cubelets', 'figures') + '_{}_specfull.txt'.format(source['id'])
    if original or (not os.path.isfile(outfile)):
        get_noise_spec(source, src_basename, cube_params, original)

    # Make plot of spectrum with noise
    make_specfull(source, src_basename, cube_params, suffix=suffix, full=False)

    plt.close('all')

    print("\tDone making spectral profiles of the spectral line source {}: {}.".format(source['id'], source['name']))

    return


if __name__ == '__main__':
    pass
    # main(source, src_basename, original=None, suffix='png')
