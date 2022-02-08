import os

from astropy import constants as const
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from matplotlib import colors
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np
from reproject import reproject_interp

from modules.functions import get_info
from modules.functions import chan2freq
from modules.get_ancillary import *

HI_restfreq = 1420405751.77 * u.Hz
optical_HI = u.doppler_optical(HI_restfreq)


###################################################################

# Overlay HI contours on optical image
def make_mom0dss2(source, src_basename, cube_params, patch, opt, suffix='png'):

    outfile = src_basename.replace('cubelets', 'figures') + '_{}_mom0.{}'.format(source['id'], suffix)

    if not os.path.isfile(outfile):
        try:
            print("\tMaking DSS2 Blue optical overlaid with HI contours.")
            hdulist_hi = fits.open(src_basename + '_{}_mom0.fits'.format(str(source['id'])))
        except FileNotFoundError:
            print("\tNo mom0 fits file. Perhaps you ran SoFiA without generating moments?")
            return

        hi_reprojected, footprint = reproject_interp(hdulist_hi, opt[0].header)

        base_contour = 3 * source['rms'] * np.abs(cube_params['chan_width'].value)
        nhi19 = 2.33e20 * base_contour / (cube_params['bmaj'].value * cube_params['bmin'].value) / 1e19
        nhi_label = "N_HI = {:.1f}, {:.1f}, {:.0f}, {:.0f}e+19".format(nhi19 * 1, nhi19 * 2, nhi19 * 4, nhi19 * 8)

        fig = plt.figure(figsize=(8, 8))
        ax1 = fig.add_subplot(111, projection=WCS(opt[0].header))
        ax1.imshow(opt[0].data, cmap='viridis', vmin=np.percentile(opt[0].data, 10),
                   vmax=np.percentile(opt[0].data, 99.8), origin='lower')
        ax1.contour(hi_reprojected, cmap='Oranges', linewidths=1, levels=base_contour * 2 ** np.arange(10))
        ax1.scatter(source['ra'], source['dec'], marker='x', c='black', linewidth=0.75,
                    transform=ax1.get_transform('fk5'))
        ax1.set_title(source['name'], fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=18)
        ax1.coords['ra'].set_axislabel('RA (ICRS)', fontsize=20)
        ax1.coords['dec'].set_axislabel('Dec (ICRS)', fontsize=20)
        ax1.text (0.5, 0.05, nhi_label, ha='center', va='center', transform=ax1.transAxes,
                  color='white', fontsize=18)
        ax1.add_patch(Ellipse((0.92, 0.9), height=patch['height'], width=patch['width'], angle=cube_params['bpa'],
                              transform=ax1.transAxes, edgecolor='white', linewidth=1))

        fig.savefig(outfile, bbox_inches='tight')

        hdulist_hi.close()

    else:
        print('\t{} already exists. Will not overwrite.'.format(outfile))

    return


# Make HI grey scale image
def make_mom0(source, src_basename, cube_params, patch, opt_head, suffix='png'):

    outfile = src_basename.replace('cubelets', 'figures') + '_{}_mom0hi.{}'.format(source['id'], suffix)

    if not os.path.isfile(outfile):
        try:
            print("\tMaking HI grey scale map.")
            hdulist_hi = fits.open(src_basename + '_{}_mom0.fits'.format(str(source['id'])))
        except FileNotFoundError:
            print("\tNo mom0 fits file. Perhaps you ran SoFiA without generating moments?")
            return

        hi_reprojected, footprint = reproject_interp(hdulist_hi, opt_head)

        base_contour = 3 * source['rms'] * np.abs(cube_params['chan_width'].value)
        nhi19 = 2.33e20 * base_contour / (cube_params['bmaj'].value * cube_params['bmin'].value) / 1e19
        nhi_label = "N_HI = {:.1f}, {:.1f}, {:.0f}, {:.0f}e+19".format(nhi19 * 1, nhi19 * 2, nhi19 * 4, nhi19 * 8)

        fig = plt.figure(figsize=(8, 8))
        ax1 = fig.add_subplot(111, projection=WCS(opt_head))
        im = ax1.imshow(hi_reprojected, cmap='gray_r', origin='lower')
        ax1.set(facecolor="white")  # Doesn't work with the color im
        ax1.contour(hi_reprojected, cmap='Oranges_r', linewidths=1.2, levels=base_contour * 2 ** np.arange(10))
        ax1.scatter(source['ra'], source['dec'], marker='x', c='white', linewidth=0.75,
                    transform=ax1.get_transform('fk5'))
        ax1.set_title(source['name'], fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=18)
        ax1.coords['ra'].set_axislabel('RA (ICRS)', fontsize=20)
        ax1.coords['dec'].set_axislabel('Dec (ICRS)', fontsize=20)
        ax1.text(0.5, 0.05, nhi_label, ha='center', va='center', transform=ax1.transAxes, fontsize=18)
        ax1.add_patch(Ellipse((0.92, 0.9), height=patch['height'], width=patch['width'], angle=cube_params['bpa'],
                              transform=ax1.transAxes, facecolor='darkorange', edgecolor='black', linewidth=1))
        cb_ax = fig.add_axes([0.91, 0.11, 0.02, 0.76])
        cbar = fig.colorbar(im, cax=cb_ax)
        cbar.set_label("HI Intensity [Jy/beam*Hz]", fontsize=18)

        fig.savefig(outfile, bbox_inches='tight')

        hdulist_hi.close()

    else:
        print('\t{} already exists. Will not overwrite.'.format(outfile))

    return


# Make HI significance image
def make_snr(source, src_basename, cube_params, patch, opt_head, suffix='png'):

    outfile = src_basename.replace('cubelets', 'figures') + '_{}_snr.{}'.format(source['id'], suffix)

    if not os.path.isfile(outfile):
        try:
            print("\tMaking pixel SNR map.")
            hdulist_snr = fits.open(src_basename + '_{}_snr.fits'.format(str(source['id'])))
        except FileNotFoundError:
            print("\tNo SNR fits file. Perhaps you ran SoFiA without generating moments?")
            return

        hdulist_hi = fits.open(src_basename + '_{}_mom0.fits'.format(str(source['id'])))
        snr_reprojected, footprint = reproject_interp(hdulist_snr, opt_head)
        hi_reprojected, footprint = reproject_interp(hdulist_hi, opt_head)

        base_contour = 3 * source['rms'] * np.abs(cube_params['chan_width'].value)
        nhi19 = 2.33e20 * base_contour / (cube_params['bmaj'].value * cube_params['bmin'].value) / 1e19

        wa_cmap = colors.ListedColormap(['w', 'royalblue', 'limegreen', 'yellow', 'orange', 'r'])
        boundaries = [0, 1, 2, 3, 4, 5, 6]
        norm = colors.BoundaryNorm(boundaries, wa_cmap.N, clip=True)
        fig = plt.figure(figsize=(8, 8))
        ax1 = fig.add_subplot(111, projection=WCS(opt_head))
        ax1.set(facecolor="white")  # Doesn't work with the color im
        im = ax1.imshow(snr_reprojected, cmap=wa_cmap, origin='lower', norm=norm)
        ax1.contour(hi_reprojected, linewidths=2, levels=[base_contour, ], colors=['k', ])
        ax1.scatter(source['ra'], source['dec'], marker='x', c='black', linewidth=0.75,
                    transform=ax1.get_transform('fk5'))
        ax1.set_title(source['name'], fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=18)
        ax1.coords['ra'].set_axislabel('RA (ICRS)', fontsize=20)
        ax1.coords['dec'].set_axislabel('Dec (ICRS)', fontsize=20)
        ax1.text(0.5, 0.05, "N_HI = {:.1f}e+19".format(nhi19), ha='center', va='center',
                 transform=ax1.transAxes, fontsize=18)
        ax1.add_patch(Ellipse((0.92, 0.9), height=patch['height'], width=patch['width'], angle=cube_params['bpa'],
                              transform=ax1.transAxes, facecolor='gold', edgecolor='indigo', linewidth=1))
        cb_ax = fig.add_axes([0.91, 0.11, 0.02, 0.76])
        cbar = fig.colorbar(im, cax=cb_ax)
        cbar.set_label("Pixel SNR", fontsize=18)
        fig.savefig(outfile, bbox_inches='tight')
        hdulist_hi.close()

    else:
        print('\t{} already exists. Will not overwrite.'.format(outfile))

    return


# Make velocity map for object
def make_mom1(source, src_basename, cube_params, patch, opt_head, opt_view=6*u.arcmin, suffix='png', sofia=2):

    outfile = src_basename.replace('cubelets', 'figures') + '_{}_mom1.{}'.format(source['id'], suffix)

    if not os.path.isfile(outfile):

        try:
            print("\tMaking velocity map.")
            mom1 = fits.open(src_basename + '_{}_mom1.fits'.format(source['id']))
        except FileNotFoundError:
            print("\tNo mom1 fits file. Perhaps you ran SoFiA without generating moments?")
            return

        for i in range(mom1[0].data.shape[0]):
            for j in range(mom1[0].data.shape[1]):
                mom1[0].data[i][j] = (mom1[0].data[i][j] * u.Hz).to(u.km/u.s, equivalencies=optical_HI).value
                # Set crazy mom1 values to nan:
                # if (mom1[0].data[i][j] > velmax) | (mom1[0].data[i][j] < velmin):
                #     mom1[0].data[i][j] = np.nan
        mom1_reprojected, footprint = reproject_interp(mom1, opt_head)
        # mom1_reprojected[significance<2.0] = np.nan

        kinpa = source['kin_pa'] * u.deg
        v_sys = (source['freq'] * u.Hz).to(u.km/u.s, equivalencies=optical_HI).value
        w50 = (const.c * source['w50'] * u.Hz / (source['freq'] * u.Hz)).to(u.km/u.s,
                                                                            equivalencies=optical_HI).value
        w20 = (const.c * source['w20'] * u.Hz / (source['freq'] * u.Hz)).to(u.km/u.s,
                                                                            equivalencies=optical_HI).value
        if sofia == 2:
            freqmin = chan2freq(source['z_min'], src_basename + '_{}_cube.fits'.format(source['id']))
            freqmax = chan2freq(source['z_max'], src_basename + '_{}_cube.fits'.format(source['id']))
        elif sofia == 1:
            freqmin = chan2freq(source['z_min'], src_basename + '_{}.fits'.format(source['id']))
            freqmax = chan2freq(source['z_max'], src_basename + '_{}.fits'.format(source['id']))
        velmax = freqmin.to(u.km / u.s, equivalencies=optical_HI).value + 5
        velmin = freqmax.to(u.km / u.s, equivalencies=optical_HI).value - 5

        v_sys_label = "v_sys = {}   W_50 = {}  W_20 = {}".format(int(v_sys), int(w50), int(w20))
        hi_pos = SkyCoord(source['ra'], source['dec'], unit='deg')

        fig = plt.figure(figsize=(8, 8))
        ax1 = fig.add_subplot(111, projection=WCS(opt_head))
        im = ax1.imshow(mom1_reprojected, cmap='RdBu_r', origin='lower')  #vmin=velmin, vmax=velmax, origin='lower')
        # ax1.contour(hi_reprojected, linewidths=1, levels=[sensitivity, ], colors=['k', ])
        if velmax - velmin > 200:
            levels = [v_sys - 100, v_sys - 50, v_sys, v_sys + 50, v_sys + 100]
            clevels = ['white', 'gray', 'black', 'gray', 'white']
        else:
            levels = [v_sys - 50, v_sys, v_sys + 50]
            clevels = ['lightgray', 'black', 'lightgray']
        ax1.contour(mom1_reprojected, colors=clevels, levels=levels, linewidths=0.6)
        # Plot HI center of galaxy
        ax1.scatter(source['ra'], source['dec'], marker='x', c='black', linewidth=0.75,
                    transform=ax1.get_transform('fk5'))
        ax1.plot([(hi_pos.ra + 0.5 * opt_view * np.sin(kinpa) / np.cos(hi_pos.dec)).deg,
                  (hi_pos.ra - 0.5 * opt_view * np.sin(kinpa) / np.cos(hi_pos.dec)).deg],
                 [(hi_pos.dec + 0.5 * opt_view * np.cos(kinpa)).deg,
                  (hi_pos.dec - 0.5 * opt_view * np.cos(kinpa)).deg],
                 c='black', linestyle='--', linewidth=0.75, transform=ax1.get_transform('icrs'))
        ax1.set_title(source['name'], fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=18)
        ax1.coords['ra'].set_axislabel('RA (ICRS)', fontsize=20)
        ax1.coords['dec'].set_axislabel('Dec (ICRS)', fontsize=20)
        ax1.text(0.5, 0.05, v_sys_label, ha='center', va='center', transform=ax1.transAxes,
                 color='black', fontsize=18)
        ax1.add_patch(Ellipse((0.92, 0.9), height=patch['height'], width=patch['width'], angle=cube_params['bpa'],
                              transform=ax1.transAxes, edgecolor='darkred', linewidth=1))
        cb_ax = fig.add_axes([0.91, 0.11, 0.02, 0.76])
        cbar = fig.colorbar(im, cax=cb_ax)
        # cbar.set_label("Barycentric Optical Velocity [km/s]", fontsize=18)
        cbar.set_label("Optical velocity [km/s]", fontsize=18)

        fig.savefig(outfile, bbox_inches='tight')

        mom1.close()

    else:
        print('\t{} already exists. Will not overwrite.'.format(outfile))

    return


# Overlay HI contours on false color optical image
def make_panstarrs(source, src_basename, cube_params, patch, color_im, opt_head, suffix='png'):

    outfile = src_basename.replace ('cubelets', 'figures') + '_{}_mom0pstr.{}'.format (source['id'], suffix)

    if not os.path.isfile(outfile):
        print("\tMaking PanSTARRS image overlaid with HI contours.")
        hdulist_hi = fits.open(src_basename + '_{}_mom0.fits'.format(str(source['id'])))
        hi_reprojected, footprint = reproject_interp(hdulist_hi, opt_head)

        base_contour = 3 * source['rms'] * np.abs(cube_params['chan_width'].value)
        nhi19 = 2.33e20 * base_contour / (cube_params['bmaj'].value * cube_params['bmin'].value) / 1e19
        nhi_label = "N_HI = {:.1f}, {:.1f}, {:.0f}, {:.0f}e+19".format (nhi19 * 1, nhi19 * 2, nhi19 * 4, nhi19 * 8)

        fig = plt.figure(figsize=(8, 8))
        ax1 = fig.add_subplot(111, projection=WCS(opt_head))
        # ax1.set_facecolor("darkgray")   # Doesn't work with the color im
        ax1.imshow(color_im, origin='lower')
        ax1.contour(hi_reprojected, cmap='Oranges', linewidths=1, levels=base_contour * 2 ** np.arange(10))
        ax1.scatter(source['ra'], source['dec'], marker='x', c='white', linewidth=0.75,
                    transform=ax1.get_transform('fk5'))
        ax1.set_title(source['name'], fontsize=20)
        ax1.tick_params(axis='both', which='major', labelsize=18)
        ax1.coords['ra'].set_axislabel('RA (ICRS)', fontsize=20)
        ax1.coords['dec'].set_axislabel('Dec (ICRS)', fontsize=20)
        ax1.text(0.5, 0.05, nhi_label, ha='center', va='center', transform=ax1.transAxes,
                 color='white', fontsize=18)
        ax1.add_patch(Ellipse((0.92, 0.9), height=patch['height'], width=patch['width'], angle=cube_params['bpa'],
                              transform=ax1.transAxes, edgecolor='lightgray', linewidth=1))
        fig.savefig(outfile, bbox_inches='tight')
    else:
        print('\t{} already exists. Will not overwrite.'.format (outfile))

    return


# Make pv plot for object
def make_pv(source, src_basename, cube_params, suffix='png'):

    outfile = src_basename.replace ('cubelets', 'figures') + '_{}_pv.{}'.format (source['id'], suffix)

    if not os.path.isfile(outfile):
        try:
            print("\tMaking pv slice.")
            pv = fits.open(src_basename + '_{}_pv.fits'.format(str(source['id'])))
        except FileNotFoundError:
            print("\tNo pv fits file.  Perhaps you ran SoFiA-2 which doesn't produce this yet?")
            return

        wcs_pv = WCS(pv[0].header)
        ang1, freq1 = wcs_pv.wcs_pix2world(0, 0, 0)
        ang2, freq2 = wcs_pv.wcs_pix2world(pv[0].header['NAXIS1'] - 1, pv[0].header['NAXIS2'] - 1, 0)
        pv_rms = np.nanstd(pv[0].data)

        fig = plt.figure(figsize=(8, 8))
        ax1 = fig.add_subplot(111, projection=WCS(pv[0].header))
        ax1.imshow(pv[0].data, cmap='gray', aspect='auto')
        # if np.all (np.isnan (pv[0].data)): continue
        ax1.contour(pv[0].data, colors='black', levels=[-2 * pv_rms, 2 * pv_rms, 4 * pv_rms])
        ax1.autoscale(False)
        ax1.plot([0.0, 0.0], [freq1, freq2], c='orange', linestyle='--', linewidth=0.75,
                 transform=ax1.get_transform('world'))
        ax1.plot([ang1, ang2], [HI_restfreq.value, HI_restfreq.value], c='orange', linestyle='--',
                 linewidth=0.75, transform=ax1.get_transform('world'))
        ax1.set_title(source['name'], fontsize=16)
        ax1.tick_params(axis='both', which='major', labelsize=18)
        ax1.set_xlabel('Angular Offset [deg]', fontsize=16)
        ax1.set_ylabel('Frequency [Hz]', fontsize=16)
        ax1.coords[1].set_ticks_position('l')
        freq_yticks = ax1.get_yticks()  # freq auto yticks from matplotlib
        ax2 = ax1.twinx()
        vel1 = const.c.to(u.km / u.s).value * (HI_restfreq.value / freq1 - 1)
        vel2 = const.c.to(u.km / u.s).value * (HI_restfreq.value / freq2 - 1)
        ax2.set_ylim(vel2, vel1)
        ax2.set_ylabel('Optical Velocity [km/s]')
        ax1.text(0.5, 0.05, 'Kinematic PA = {:5.1f} deg'.format(source['kin_pa']), ha='center', va='center',
                 transform=ax1.transAxes, color='orange', fontsize=18)
        fig.savefig(outfile, bbox_inches='tight')
        pv.close()

    else:
        print('\t{} already exists. Will not overwrite.'.format (outfile))

    return


def main(source, src_basename, opt_view=6*u.arcmin, opt_pixels=900, suffix='png', sofia=2):

    print("\n\tStart making spatial images of the spectral line source {}: {}.".format(source['id'], source['name']))

    # Get beam information from the source cubelet
    if sofia == 2:
        cube_params = get_info(src_basename + '_{}_cube.fits'.format(source['id']))
    elif sofia == 1:
        cube_params = get_info (src_basename + '_{}.fits'.format (source['id']))

    # Get the position of the source to retrieve an optical image
    hi_pos = SkyCoord(ra=source['ra'], dec=source['dec'], unit='deg',
                      equinox=cube_params['equinox'], frame=cube_params['frame'])

    # SkyView (and maybe other queries??) won't retrieve non-ICRS or non-J2000 coordinates every time,
    # so let's transform everything to ICRS ....Need to keep an eye on this (may have been a server issue).
    hi_pos_icrs = hi_pos.transform_to('icrs')

    # Calculate the size of the optical image for the moment maps
    Xc = source['x']
    Yc = source['y']
    Xmin = source['x_min']
    Ymin = source['y_min']
    Xmax = source['x_max']
    Ymax = source['y_max']
    Xsize = np.array([((Xmax - Xc) * cube_params['cellsize']).to(u.arcmin).value,
                      ((Xc - Xmin) * cube_params['cellsize']).to(u.arcmin).value])
    Ysize = np.array([((Ymax - Yc) * cube_params['cellsize']).to(u.arcmin).value,
                      ((Yc - Ymin) * cube_params['cellsize']).to(u.arcmin).value])
    if np.any(Xsize > opt_view.value / 2) | np.any(Ysize > opt_view.value / 2):
        opt_view = np.max([Xsize, Ysize]) * 2 * 1.05 * u.arcmin
        print("\tImage size bigger than default. Now {:.2f} arcmin".format(opt_view.value))

    # Get optical images, based on the HI position and given image size.
    dss2 = get_dss2(hi_pos_icrs, opt_view, opt_pixels)
    pstar_view = opt_view
    if opt_view > 8*u.arcmin:
        pstar_view = 8*u.arcmin
        print("\tAdjusted viewing size greater than 8 arcmin.  This is the SIP imposed limit on\n" \
              "\t\tPanSTARRS images (they are much bigger than DSS2).")
    pstar_im, pstar_head = get_panstarrs(hi_pos_icrs, opt_view=pstar_view)

    # Temporarily replace with ICRS ra/dec for plotting purposes in the rest (won't change catalog file.):
    source['ra'] = hi_pos_icrs.ra.deg
    source['dec'] = hi_pos_icrs.dec.deg

    # Calculate the size of the beam (plotted as a fraction of the image size)
    patch_height = (cube_params['bmaj'] / opt_view).decompose()
    patch_width = (cube_params['bmin'] / opt_view).decompose()
    patch = {'width': patch_width, 'height': patch_height}

    # Create the first optical overlay figure:
    if dss2:
        make_mom0dss2(source, src_basename, cube_params, patch, dss2, suffix=suffix)
        opt_head = dss2[0].header
        dss2.close()

    # In theory can have different sizes for the panstarrs and dss2 images, but then need to recalculate the patch
    # to represent the beam size.  This isn't done rigorously yet.
    if pstar_im:
        patch_height = (cube_params['bmaj'] / pstar_view).decompose()
        patch_width = (cube_params['bmin'] / pstar_view).decompose()
        patch_pstar = {'width': patch_width, 'height': patch_height}
        make_panstarrs(source, src_basename, cube_params, patch, pstar_im, pstar_head, suffix=suffix)

    # Use dss2 image as the base for regridding the HI since it is relatively small (although set by the number of pixels...
    # need to change this to take into account pixel scale to be rigorous.
    # If dss2 doesn't exist, use panstarrs image as the base for regridding.
    if not dss2:
        opt_head = pstar_head
        patch = patch_pstar

    # Make the rest of the images if there is an optical image to regrid to.
    # Could change this to make images no matter what...
    if dss2 or pstar_im:
        make_mom0(source, src_basename, cube_params, patch, opt_head, suffix=suffix)
        make_snr(source, src_basename, cube_params, patch, opt_head, suffix=suffix)
        make_mom1(source, src_basename, cube_params, patch, opt_head, opt_view=opt_view, suffix=suffix, sofia=2)

    # Make pv if it was created (only in SoFiA-1); not dependent on optical image.
    make_pv(source, src_basename, cube_params, suffix=suffix)

    plt.close('all')

    print("\tDone making spatial images of the spectral line source {}: {}.".format(source['id'], source['name']))

    return True


if __name__ == '__main__':

    main(source, src_basename, opt_view=6*u.arcmin, opt_pixels=900, suffix='png')
