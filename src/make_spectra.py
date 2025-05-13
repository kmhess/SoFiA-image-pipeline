# from datetime import datetime
import os

from astropy import constants as const
from astropy.io import ascii, fits
from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pkg_resources  # part of setuptools

from src.modules.functions import chan2freq, chan2vel, get_info, get_subcube, felo2vel, line_lookup, make_hist_arr
from src.modules.logger import Logger

logger = Logger.get_logger()

version = pkg_resources.require("SoFiA-image-pipeline")[0].version

###################################################################


# Save full spectrum to a txt file:
def get_noise_spec(source, src_basename, cube_params, original=None):

    outfile = src_basename.replace('cubelets', 'figures') + '_{}_specfull.txt'.format(source['id'])

    if not os.path.isfile(outfile):

        try:
            if not original:
                logger.warning("\tOriginal data cube not provided: making spectrum of subcube with noise.")
                fits_file = src_basename + '_{}_cube.fits'.format(source['id'])
                cube = fits.getdata(fits_file)
                mask = fits.getdata(src_basename + '_{}_mask.fits'.format(source['id']))
                spec_template = ascii.read(src_basename + '_{}_spec.txt'.format(source['id']),
                                            names=['chan', 'col2', 'f_sum', 'n_pix'])
                channels = spec_template['chan']
            else:
                fits_file = original
                cube = get_subcube(source, original)
                if os.path.isfile(original) and os.path.isfile(original[:-5] + '_mask.fits'):
                    logger.info("\tOriginal data cube provided: making full spectrum image with noise.")
                    mask = get_subcube(source, original[:-5] + '_mask.fits')
                elif os.path.isfile(original) and os.path.isfile(src_basename.split('_cubelets')[0] + '_mask.fits'):
                    logger.info("\tOriginal data cube provided, original file name differs from catalog file: Making full spectrum image with noise.")
                    mask = get_subcube(source, src_basename.split('_cubelets')[0] + '_mask.fits')
                spec_template = None
                channels = np.asarray(range(cube.shape[0]))
        except:
            logger.warning("\tWrong name provided for original file, or original mask file doesn't exist, so can't generate a *_specfull.txt with noise.")
            return

        mask2d = np.sum(mask, axis=0)
        spectrum = np.nansum(cube[:, mask2d != 0], axis=1)
        n_pix = 0 * channels + np.sum(mask2d != 0)

        with open('temp.txt', 'w') as f:
            f.write("# Integrated source spectrum with noise\n")
            f.write("# Creator: SoFiA-image-pipeline (SIP) version {}\n".format(version))
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
            elif 'VRAD' in cube_params['spec_axis']:
                velocities = spec_template['col2'] if spec_template else chan2vel(channels, fits_file)
                ascii.write([channels, velocities, spectrum, n_pix], 'temp2.txt', format='fixed_width_two_line',
                            names=['chan', 'v_rad', 'f_sum', 'n_pix'])
            else:
                velocities = spec_template['col2'] if spec_template else chan2vel(channels, fits_file)
                ascii.write([channels, velocities, spectrum, n_pix], 'temp2.txt', format='fixed_width_two_line',
                            names=['chan', 'velo', 'f_sum', 'n_pix'])
        os.system("cat temp.txt temp2.txt > {}".format(outfile))
        os.system("rm temp.txt temp2.txt")

        return


# Make full spectrum plot:
def make_specfull(source, src_basename, cube_params, original, spec_line=None, suffix='png'):

    outfile2 = src_basename.replace('cubelets', 'figures') + '_{}_specfull.{}'.format(source['id'], suffix)
    specfile = src_basename + '_{}_spec_aperture.txt'.format(source['id'])

    if original or (not os.path.isfile(specfile)):
        specfile = src_basename.replace('cubelets', 'figures') + '_{}_specfull.txt'.format(source['id'])

    logger.info('\tUsing {} to make aperture spectrum plot.'.format(specfile))

    long_format = 200

    if not os.path.isfile(outfile2):

        # Get frequency information for spectral line in question:
        line = line_lookup(spec_line)

        try:
            logger.info("\tMaking {} spectrum plot including noise.".format(spec_line))
            if 'freq' in source.colnames:
                # Calculate source quantities for labels
                v_sys = (source['freq'] * u.Hz).to(u.km/u.s, equivalencies=line['convention']).value
                spec_z = (line['restfreq'].to(u.Hz) - source['freq'] * u.Hz) / (source['freq'] * u.Hz).decompose()
                z_label = r"$z_\mathrm{{{0:s}}}$ = {1:.5f}".format(line['name'], spec_z.value)
                # SoFiA-2 puts out frequency w20/w50 in Hz units
                w50_vel = (const.c * source['w50'] / (source['freq'])).to(u.km/u.s).value
                w20_vel = (const.c * source['w20'] / (source['freq'])).to(u.km/u.s).value
                # Calculate spectral axes quantities for plotting
                spec = ascii.read(specfile, names=['chan', 'freq', 'f_sum', 'n_pix'])
                optical_velocity = (source['freq'] - spec['freq'])/spec['freq'] * const.c.to(u.km/u.s).value
                maskmin = (source['freq'] - spec['freq'][spec['chan'] == source['z_min']]) / source['freq'] * \
                                                                                            const.c.to(u.km / u.s).value
                maskmax = (source['freq'] - spec['freq'][spec['chan'] == source['z_max']]) / source['freq'] * \
                                                                                            const.c.to(u.km / u.s).value
                v_sys = 0
            else:
                # Calculate source quantities for labels
                if 'v_rad' in source.colnames:
                    line['rad_opt'] = 'Radio'
                    v_sys = (source['v_rad'] * u.m / u.s).to(u.km / u.s).value
                    spec_z = ''
                    z_label = ''
                elif 'v_opt' in source.colnames:
                    v_sys = (source['v_opt'] * u.m / u.s).to(u.km / u.s).value
                    spec_z = (source['v_opt'] * u.m / u.s / const.c).decompose()
                    z_label = r"$z_\mathrm{{{0:s}}}$ = {1:.5f}".format(line['name'], spec_z.value)
                else:
                    v_sys = (source['v_app'] * u.m / u.s).to(u.km / u.s).value
                    spec_z = (source['v_app'] * u.m / u.s / const.c).decompose()
                    z_label = r"$z_\mathrm{{app}}$ = {:.5f}".format(spec_z.value)
                # SoFiA-2 puts out velocity w20/w50 in pixel units. https://github.com/SoFiA-Admin/SoFiA-2/issues/63
                w50_vel = (source['w50'] * u.m / u.s).to(u.km / u.s).value
                w20_vel = (source['w20'] * u.m / u.s).to(u.km / u.s).value
                # Calculate spectral axes quantities for plotting. Force velocity column to common name.
                spec = ascii.read(specfile, names=['chan', 'velo', 'f_sum', 'n_pix'])
                optical_velocity = (spec['velo'] * u.m / u.s).to(u.km / u.s).value
                maskmin = (spec['velo'][spec['chan'] == source['z_min']] * u.m / u.s).to(u.km / u.s,
                                                                                         equivalencies=line['convention']).value
                maskmax = (spec['velo'][spec['chan'] == source['z_max']] * u.m / u.s).to(u.km / u.s,
                                                                                         equivalencies=line['convention']).value
            v_sys_label = "$W_{{20}}$ = {} km/s".format(int(w20_vel))
            # if original or len(spec) >= long_format:
            #     v_sys_label += "  $W_{{20}}$ = {} km/s".format(int(w20_vel))
            # else:
            #     v_sys_label += " km/s"
            if 'snr' in source.colnames:
                v_sys_label += ",  SNR = {:.1f}".format(source['snr'])

        except FileNotFoundError:
            logger.warning("\tNo existing _specfull.txt file. Perhaps there is no cube to generate one, or need to specify original.")
            fig2, ax2_spec, outfile2 = None, None, None
            return fig2, ax2_spec, outfile2

        # Get spec units (Jy or Jy/beam) to check whether the division by the beam area has been made
        ll = 0
        try:
            while not ('f_sum' in spec.meta['comments'][ll] and 'chan' in spec.meta['comments'][ll] and
                    'n_pix' in spec.meta['comments'][ll]):
                ll += 1
            specunits = (spec.meta['comments'][ll+1].split()[spec.meta['comments'][ll].split().index('f_sum')])
        except IndexError:
            # If spectral plot is not made by SoFiA or does not have units for f_sum, assume the beam not accounted for
            specunits = 'Jy/beam'

        if specunits == 'Jy/beam':
            flux_dens = spec['f_sum'] / cube_params['pix_per_beam']
        elif specunits == 'Jy':
            flux_dens = spec['f_sum']
        y_error = source['rms'] * np.sqrt(spec['n_pix'] / cube_params['pix_per_beam'])

        # Could be more clever about picking the size of the figure when there are a lot of channels. Leave for later.
        if original or len(spec) >= long_format:
            fig2 = plt.figure(figsize=(14, 4))
        else:
            fig2 = plt.figure(figsize=(9, 4))

        ax2_spec = fig2.add_subplot(111)
        ax2_spec.plot([np.min(optical_velocity) - 10, np.max(optical_velocity) + 10], [0, 0], '--', color='gray')
        # If there are lots of channels, don't plot errors (too crowded and can tell from noise.) Cut off currently arbitrary.
        if len(spec) <= 100:
            opt_vel, f_sum, y_err = make_hist_arr(xx=optical_velocity, yy=flux_dens, yy_err=y_error)
            ax2_spec.errorbar(opt_vel, f_sum, elinewidth=0.75, yerr=y_err, capsize=1)
        elif len(spec) <= 200:
            opt_vel, f_sum, y_err = make_hist_arr(xx=optical_velocity, yy=flux_dens, yy_err=y_error * 0)
            ax2_spec.errorbar(opt_vel, f_sum, elinewidth=0.75, yerr=y_err, capsize=0)
        else:
            logger.info("\tInput *_specfull.txt is >=200 channels; expanding figure, not including error bars (noise should be indicative).")
            ax2_spec.plot(optical_velocity, flux_dens)
        ax2_spec.text(0.05, 0.90, z_label, ha='left', va='center', transform=ax2_spec.transAxes, color='black', fontsize=17)
        ax2_spec.text(0.5, 0.06, v_sys_label, ha='center', va='center', transform=ax2_spec.transAxes, color='black', fontsize=17)
        ax2_spec.set_title(source['name'], fontsize=20)
        ax2_spec.set_xlim(np.min(optical_velocity) - 5, np.max(optical_velocity) + 5)
        ax2_spec.set_ylabel("Integrated Flux [Jy]", fontsize=17)
        if 'freq' in source.colnames:
            ax2_spec.set_xlabel("Rest frame velocity [km/s]", fontsize=17)
        elif line['rad_opt'] == 'Optical':
            ax2_spec.set_xlabel("{} cz [km/s]".format(cube_params['spec_sys'].capitalize()), fontsize=17)
        else:
            ax2_spec.set_xlabel("{} {} Recessional Velocity [km/s]".format(cube_params['spec_sys'].capitalize(), 
                                                                           line['rad_opt']), fontsize=17)
        ax2_spec.tick_params(axis='both', which='major', labelsize=16, length=5, width=1.8)
        ax2_spec.autoscale(False)
        if not original or len(spec) < long_format:
            ax2_spec.xaxis.set_major_locator(plt.MaxNLocator(7))
        if 'freq' in source.colnames:
            ax2b_spec = ax2_spec.twiny()
            freq1 = (spec['freq'][-1] * u.Hz).to(u.MHz)
            freq2 = (spec['freq'][0] * u.Hz).to(u.MHz)
            ax2b_spec.set_xlabel('Frequency [MHz]', fontsize=17)
            if freq1 >= 2*u.GHz:
                freq1 = freq1.to(u.GHz)
                freq2 = freq2.to(u.GHz)
                ax2b_spec.set_xlabel('Frequency [GHz]', fontsize=17)
            ax2b_spec.set_xlim(freq1.value, freq2.value)
            ax2b_spec.tick_params(labelsize=16, length=5, width=1.8)
            ax2b_spec.ticklabel_format(style='plain', useOffset=False)
            if not original or len(spec) < long_format:
                ax2b_spec.xaxis.set_major_locator(plt.MaxNLocator(7))
        spectrumJy = spec["f_sum"] / cube_params['pix_per_beam']
        galspec_max = np.nanmax(spectrumJy[np.where(spec['chan'] == source['z_min'])[0][0]:
                                           np.where(spec['chan'] == source['z_max'])[0][0]+1])
        # Minimum within mask could be positive or negative, but to consider values around 0.
        galspec_min = -1 * np.abs(np.nanmin(spectrumJy[np.where(spec['chan'] == source['z_min'])[0][0]:
                                                       np.where(spec['chan'] == source['z_max'])[0][0]+1])).value

        # Plot limit of SoFiA mask
        ymin, ymax = ax2_spec.get_ylim()
        ax2_spec.plot([maskmin, maskmin], [0.95*ymin, 0.95*ymax], ':', color='gray')
        ax2_spec.plot([maskmax, maskmax], [0.95*ymin, 0.95*ymax], ':', color='gray')

        ax2_spec.plot([v_sys, v_sys], np.array([-0.05, 0.05])*(ymax-ymin), color='gray')

        # # Condition from Apertif experience that if the RFI is *really* bad, plot based on strength of HI profile
        # ax_margin_percent = 0.15
        # if (ymax > 5.*galspec_max) | (ymin < np.nanmin([-2.*np.abs(galspec_min), -5.*np.nanstd(spectrumJy).value])):
        #     print("\tWARNING: Suspect there is a lot of noise in the full spectrum?  Trying to adjust"
        #           " y-axis to source.")
        #     ax_margin = (1. + ax_margin_percent) * np.array([np.nanmin([-2.*np.abs(galspec_min),
        #                                                                 -3*np.nanstd(spectrumJy).value]), 3.*galspec_max])
        #     print(ax_margin)
        #     ax2_spec.set_ylim(ax_margin)

    else:
        logger.warning('\t{} already exists. Will not overwrite.'.format(outfile2))
        fig2, ax2_spec, outfile2 = None, None, None

    return fig2, ax2_spec, outfile2


# Make SoFiA masked spectrum plot (no noise):
def make_spec(source, src_basename, cube_params, spec_line=None, suffix='png'):

    outfile1 = src_basename.replace('cubelets', 'figures') + '_{}_spec.{}'.format(source['id'], suffix)
    fits_file = src_basename + '_{}_cube.fits'.format(source['id'])

    # For estimating position of z_w20, z_w50, z_wm50 which are given in pixel space:
    fits_file = src_basename + '_{}_cube.fits'.format(source['id'])

    if not os.path.isfile(outfile1):

        # Get frequency information for spectral line in question:
        line = line_lookup(spec_line)

        try:
            logger.info("\tMaking {} SoFiA masked spectrum plot.".format(spec_line))
            if 'freq' in source.colnames:
                # Calculate source quantities for labels
                v_sys = (source['freq'] * u.Hz).to(u.km/u.s, equivalencies=line['convention']).value
                spec_z = (line['restfreq'].to(u.Hz) - source['freq'] * u.Hz) / (source['freq'] * u.Hz).decompose()
                z_label = r"$z_\mathrm{{{0:s}}}$ = {1:.5f}".format(line['name'], spec_z.value)
                # SoFiA-2 puts out frequency w20/w50 in Hz units
                w50_vel = (const.c * source['w50'] / (source['freq'])).to(u.km/u.s).value
                w20_vel = (const.c * source['w20'] / (source['freq'])).to(u.km/u.s).value
                # Calculate spectral axes quantities for plotting
                spec = ascii.read(src_basename + '_{}_spec.txt'.format(source['id']),
                                  names=['chan', 'freq', 'f_sum', 'n_pix'])
                optical_velocity = (source['freq'] - spec['freq'])/spec['freq'] * const.c.to(u.km/u.s).value
                v_sys = 0
                if 'z_w20' in source.colnames:
                    z_w20 = chan2freq(source['z_w20'], fits_file)
                    z_w20_vel = ((source['freq'] * u.Hz - z_w20) / source['freq'] * const.c.to(u.km / u.s)).value
                    w20_min_vel = z_w20_vel - ((source['w20'] * u.Hz / 2) / (source['freq'] * u.Hz) * const.c.to(u.km / u.s)).value
                    w20_max_vel = z_w20_vel + ((source['w20'] * u.Hz / 2) / (source['freq'] * u.Hz) * const.c.to(u.km / u.s)).value
            else:
                # Calculate source quantities for labels
                if 'v_rad' in source.colnames:
                    line['rad_opt'] = 'Radio'
                    v_sys = (source['v_rad'] * u.m / u.s).to(u.km / u.s).value
                    spec_z = ''
                    z_label = ''
                elif 'v_opt' in source.colnames:
                    v_sys = (source['v_opt'] * u.m / u.s).to(u.km / u.s).value
                    spec_z = (source['v_opt'] * u.m / u.s / const.c).decompose()
                    z_label = r"$z_\mathrm{{{0:s}}}$ = {1:.5f}".format(line['name'], spec_z.value)
                else:
                    v_sys = (source['v_app'] * u.m / u.s).to(u.km / u.s).value
                    spec_z = (source['v_app'] * u.m / u.s / const.c).decompose()
                    z_label = r"$z_\mathrm{{app}}$ = {:.5f}".format(spec_z.value)
                # SoFiA-2 puts out velocity w20/w50 in pixel units. https://github.com/SoFiA-Admin/SoFiA-2/issues/63
                w50_vel = (source['w50'] * u.m / u.s).to(u.km / u.s).value
                w20_vel = (source['w20'] * u.m / u.s).to(u.km / u.s).value
                # Calculate spectral axes quantities for plotting. Force velocity column to common name.
                spec = ascii.read(src_basename + '_{}_spec.txt'.format(source['id']),
                                  names=['chan', 'velo', 'f_sum', 'n_pix'])
                optical_velocity = (spec['velo'] * u.m / u.s).to(u.km / u.s).value
                if 'z_w20' in source.colnames:
                    z_w20 = chan2vel(source['z_w20'], fits_file)
                    w20_min_vel = (z_w20 - source['w20'] * u.m / u.s / 2).to(u.km / u.s).value
                    w20_max_vel = (z_w20 + source['w20'] * u.m / u.s / 2).to(u.km / u.s).value
            if 'snr' in source.colnames:
                v_sys_label = "$W_{{20}}$ = {} km/s,  SNR = {:.1f}".format(int(w20_vel), source['snr'])
            else:
                v_sys_label = "$W_{{20}}$ = {} km/s".format(int(w20_vel))

        except FileNotFoundError:
            logger.warning("\tNo *_spec.txt file.")
            fig1, ax1_spec, outfile1 = None, None, None
            return fig1, ax1_spec, outfile1

        # Get spec units (Jy or Jy/beam) to check whether the division by the beam area has been made
        ll = 0
        while not ('f_sum' in spec.meta['comments'][ll] and 'chan' in spec.meta['comments'][ll] and
                   'n_pix' in spec.meta['comments'][ll]):
            ll += 1
        specunits = (spec.meta['comments'][ll+1].split()[spec.meta['comments'][ll].split().index('f_sum')])

        fig1 = plt.figure(figsize=(8, 4))
        ax1_spec = fig1.add_subplot(111)
        ax1_spec.plot([np.min(optical_velocity) - 10, np.max(optical_velocity) + 10], [0, 0], '--', color='gray')
        y_error = source['rms'] * np.sqrt(spec['n_pix'] / cube_params['pix_per_beam'])
        if specunits == 'Jy/beam':
            opt_vel, f_sum, y_err = make_hist_arr(xx=optical_velocity, yy=spec['f_sum'] / cube_params['pix_per_beam'], 
                                                  yy_err=y_error)
            ax1_spec.errorbar(opt_vel, f_sum, elinewidth=0.75, yerr=y_err, capsize=1)
        elif specunits == 'Jy':
            opt_vel, f_sum, y_err = make_hist_arr(xx=optical_velocity, yy=spec['f_sum'], yy_err=y_error)
            ax1_spec.errorbar(opt_vel, f_sum, elinewidth=0.75, yerr=y_err, capsize=1)
        ax1_spec.text(0.05, 0.90, z_label, ha='left', va='center', transform=ax1_spec.transAxes, color='black', fontsize=17)
        ax1_spec.text(0.5, 0.06, v_sys_label, ha='center', va='center', transform=ax1_spec.transAxes, color='black', fontsize=17)
        ax1_spec.set_title(source['name'], fontsize=20)
        ax1_spec.set_xlim(np.min(optical_velocity) - 5, np.max(optical_velocity) + 5)
        ax1_spec.set_ylabel("Integrated Flux [Jy]", fontsize=17)
        if 'freq' in source.colnames:
            ax1_spec.set_xlabel("Rest frame velocity [km/s]", fontsize=17)
        elif line['rad_opt'] == 'Optical':
            ax1_spec.set_xlabel("{} cz [km/s]".format(cube_params['spec_sys'].capitalize()), fontsize=17)
        else:
            ax1_spec.set_xlabel("{} {} Recessional Velocity [km/s]".format(cube_params['spec_sys'].capitalize(), 
                                                                           line['rad_opt']), fontsize=17)
        ax1_spec.tick_params(axis='both', which='major', labelsize=16, length=5, width=1.8)
        ax1_spec.autoscale(False)
        ax1_spec.xaxis.set_major_locator(plt.MaxNLocator(7))
        if 'freq' in source.colnames:
            ax1b_spec = ax1_spec.twiny()
            freq1 = (spec['freq'][-1] * u.Hz).to(u.MHz)
            freq2 = (spec['freq'][0] * u.Hz).to(u.MHz)
            ax1b_spec.set_xlabel('Frequency [MHz]', fontsize=17)
            if freq1 >= 2*u.GHz:
                freq1 = freq1.to(u.GHz)
                freq2 = freq2.to(u.GHz)
                ax1b_spec.set_xlabel('Frequency [GHz]', fontsize=17)
            ax1b_spec.set_xlim(freq1.value, freq2.value)
            ax1b_spec.ticklabel_format(style='plain', useOffset=False)
            ax1b_spec.tick_params(labelsize=16, length=5, width=1.8)
            ax1b_spec.xaxis.set_major_locator(plt.MaxNLocator(7))
        # Plot limit of w20 width
        ymin, ymax = ax1_spec.get_ylim()
        if 'z_w20' in source.colnames:
            ax1_spec.plot([w20_min_vel, w20_min_vel], [0.95*ymin, 0.95*ymax], ':', color='red')
            ax1_spec.plot([w20_max_vel, w20_max_vel], [0.95*ymin, 0.95*ymax], ':', color='red')
        ax1_spec.plot([v_sys, v_sys], np.array([-0.05, 0.05])*(ymax-ymin), color='gray')
    else:
        logger.warning('\t{} already exists. Will not overwrite.'.format(outfile1))
        fig1, ax1_spec, outfile1 = None, None, None

    return fig1, ax1_spec, outfile1


def main(source, src_basename, original=None, spec_line=None, suffix='png', beam=None):

    logger.info("\tStart making spectral profiles")

    # Get beam information from the source cubelet
    try:
        cube_params = get_info(src_basename + '_{}_cube.fits'.format(source['id']), beam)
    except FileNotFoundError:
        try:
            cube_params = get_info(src_basename + '_{}_mom0.fits'.format(source['id']), beam)
        except FileNotFoundError:
            logger.error("\tNo cubelet or mom0 to match source {}."
                  " Can't determine coordinate system to plot spectrum.\n".format(source['id']))
            return

    # Make plot of SoFiA masked spectrum
    fig1, ax1_spec, outfile1 = make_spec(source, src_basename, cube_params, spec_line=spec_line, suffix=suffix)

    # Make text file of spectrum with noise; use full frequency range of original cube if provided:
    # Can be a bit more precise here in the output options/specification.
    sofia_aper_spec = src_basename + '_{}_spec_aperture.txt'.format(source['id'])
    specfull_file = src_basename.replace('cubelets', 'figures') + '_{}_specfull.txt'.format(source['id'])
    if original or ((not os.path.isfile(sofia_aper_spec)) and (not os.path.isfile(specfull_file))):
        get_noise_spec(source, src_basename, cube_params, original)

    # Make plot of spectrum with noise
    fig2, ax2_spec, outfile2 = make_specfull(source, src_basename, cube_params, original, spec_line=spec_line,
                                             suffix=suffix)
    if outfile1 and outfile2:
        ymin = min([ax1_spec.get_ylim()[0], ax2_spec.get_ylim()[0]])
        ymax = max([ax1_spec.get_ylim()[1], ax2_spec.get_ylim()[1]])
        ax1_spec.set_ylim([ymin, ymax])
        ax2_spec.set_ylim([ymin, ymax])
    if outfile1:
        fig1.savefig(outfile1, bbox_inches='tight')
    if outfile2:
        fig2.savefig(outfile2, bbox_inches='tight')
    plt.close('all')

    logger.info("\tDone making spectral profiles.")

    return


if __name__ == '__main__':
    pass
    # main(source, src_basename, original=None, suffix='png')
