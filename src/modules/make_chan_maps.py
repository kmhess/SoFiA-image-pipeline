import os

from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np

from src.modules.functions import get_info, get_wcs_info
from src.modules.functions import chan2freq

from src.modules.logger import Logger

logger = Logger.get_logger()


###################################################################


def main(source, src_basename, suffix='png', beam=None, noid=False, opt_head=None, patch=None):
        #  original, opt_view=6*u.arcmin, snr_range=[2, 3], spec_line=None, 

    logger.info("\tMaking CHANNEL maps.")
    id_label = ''
    if noid == False:
        id_label = ', #{}'.format(source['id'])

    base_contour = source['rms'] * 2

    # pos_x and pos_y were set in make_images.py  To make code stand-alone, need to repeat that calculation here
    # print(source['pos_x'], source['pos_y'])

    # opt_head was calculated in make_images.py  Keep consistent my just taking it in from there
    # print(opt_head)

    infile = src_basename + '_{}_cube.fits'.format(source['id'])
    outfile = src_basename.replace('cubelets', 'figures') + '_{}_chan_maps.pdf'.format(source['id'])

    if not os.path.isfile(outfile):
        # Get beam information from the source cubelet
        try:
            cube_params = get_info(infile, beam)
        except FileNotFoundError:
            logger.error("\tNo cubelet to match source {}. Can't make channel maps.\n".format(source['id']))
            return
        hdulist_hi = fits.open(infile)

        hiwcs, cubew = get_wcs_info(infile)
        owcs = WCS(opt_head)

        colors_noise = plt.cm.gray(np.linspace(0, 1, 256))
        colors_galaxy = plt.cm.afmhot(np.linspace(1, 0.4, 256))
        all_colors = np.vstack((colors_noise, colors_galaxy))
        pvd_map = colors.LinearSegmentedColormap.from_list('pvd_map', all_colors)
        divnorm = colors.TwoSlopeNorm(vmin=-3*source['rms'], vcenter=+3*source['rms'], vmax=15*source['rms'])

        if 'l' in source.colnames:
            x_coord, y_coord = 'glon', 'glat'
            x_label, y_label = r'$\it{{l}}$ [deg]', r'$\it{{b}}$ [deg]'
        else:
            x_coord, y_coord = 'ra', 'dec'
            x_label, y_label = 'Right Ascension (ICRS)', 'Declination (ICRS)'

        xsize = 4
        ysize = 5
        remove_margin = 10
        end = hdulist_hi[0].data.shape[0] - remove_margin
        print('FREQUENCIES ARE WRONG SO FAR: ABSOLUTE VS RELATIVE CHANNEL')
        chans_ghz = chan2freq(np.array(range(remove_margin, end, 1)), infile).to(u.GHz)

        chan_per_panel = xsize * ysize

        # Save channel maps to a multi-page pdf
        with PdfPages(outfile) as pdf:
            for c, c_ghz in zip(range(hdulist_hi[0].data.shape[0])[remove_margin:end], chans_ghz):
                chan = hdulist_hi[0].data[c, :, :]
                panel = int((c - remove_margin) % chan_per_panel)
                column = panel % xsize
                row = int(((panel - column) / xsize) % chan_per_panel)

                if panel == 0:
                    if end - c < chan_per_panel:
                        if (end - c) % xsize == 0:
                            ysize = int((end - c) / xsize)
                        else:
                            ysize = int((end - c) / xsize) + 1
                    fig = plt.figure(figsize=(xsize*3, ysize*3))
                    gs = fig.add_gridspec(ysize, xsize, hspace=0, wspace=0)
                    axs = gs.subplots(sharex='col', sharey='row', subplot_kw=dict(projection=owcs))
                    # Some brute force relocating of labels because tight_layout screws up panel spacing:
                    fig.suptitle(source['name'] + id_label, fontsize=24, y=0.91 + (0.09-0.09*(ysize/5))*1.3)
                    fig.supxlabel(x_label, fontsize=22, y=0.05 - (0.05 * (5/ysize)-0.05)*1.4)
                    fig.supylabel(y_label, fontsize=22, x=-0.01 - (0.01 * (4/xsize)-0.01)*12)

                if axs.ndim == 1:
                    axs = axs[None,:]
                axs[row,column].imshow(chan, cmap=pvd_map, norm=divnorm, origin='lower', 
                                       transform=axs[row,column].get_transform(cubew))
                # axs[row,column].text(0.5, 0.5, '{:.2f}'.format(c_ghz), fontsize=14)
                axs[row,column].text(0.5, 0.5, 'XXX GHz'.format(c_ghz), fontsize=14)
                axs[row,column].tick_params(axis='both', which='major', labelsize=14, length=6, width=2)
                axs[row,column].set(facecolor="white")  # Doesn't work with the color im
                # Plot positive contours
                try:
                    axs[row,column].contour(chan, cmap='Oranges_r', linewidths=1.2, 
                                            levels=base_contour * 2 ** np.arange(10),
                                            transform=axs[row,column].get_transform(cubew))
                except:
                    pass
                # Plot negative contours
                try:
                    axs[row,column].contour(chan, cmap='YlOrBr_r', linewidths=1.2, linestyles='dashed',
                                            levels=-base_contour * 2 ** np.arange(10, -1, -1),
                                            transform=axs[row,column].get_transform(cubew))
                except:
                    pass
                # Plot beam in the first panel
                if panel == 0:
                    axs[row,column].add_patch(Ellipse((0.92, 0.9), height=patch['height'], width=patch['width'],
                                                      angle=cube_params['bpa'], transform=axs[row,column].transAxes,
                                                      edgecolor='black', linewidth=1, zorder=99))

                for ax in axs[:-1].flatten():
                    ax.tick_params(axis='x', which='both', labelbottom=False)
                for ax in axs[:,1:].flatten():
                    ax.tick_params(axis='y', which='both', labelleft=False)
                if column == 0:
                    axs[row,column].coords[y_coord].set_axislabel('')
                if row == ysize - 1:
                    axs[row,column].coords[x_coord].set_axislabel('')
                    axs[row,column].coords[x_coord].set_ticklabel(rotation=30, ha='right', rotation_mode='anchor')

                axs[row,column].set_xlim(0, opt_head['NAXIS1'])
                axs[row,column].set_ylim(0, opt_head['NAXIS2'])

                if ((row+1) * (column+1) == chan_per_panel) | (c == end - 1): #hdulist_hi[0].data.shape[0] - 1):
                    column += 1
                    if column != xsize: #and row != 0:
                        while column < xsize:
                            axs[row,column].coords[x_coord].set_axislabel('')
                            axs[row,column].coords[x_coord].set_ticklabel_visible(False)
                            # axs[row,column].coords[x_coord].set_ticklabel(rotation=30, ha='right', 
                            #                                               rotation_mode='anchor')
                            axs[row,column].set_xlim(0, opt_head['NAXIS1'])
                            axs[row,column].set_ylim(0, opt_head['NAXIS2'])
                            column += 1
                    pdf.savefig(fig, bbox_inches='tight')
                    plt.close()

        hdulist_hi.close()

    else:
        logger.warning('\t{} already exists. Will not overwrite.'.format(outfile))

    return
