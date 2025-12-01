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
        boundary = 1   # Number of source-free channels in channel map
 
        # SoFiA indexed from 1; FITS headers indexed from 1; python arrays indexed from 0
        new_py_zmin = int(source['z_min']-1 + hdulist_hi[0].header['CRPIX3']-1 - boundary)
        new_py_zmax = int(source['z_max']-1 + hdulist_hi[0].header['CRPIX3']-1 + boundary)
        new_py_zmin = np.max([0, new_py_zmin])
        new_py_zmax = np.min([new_py_zmax, len(hdulist_hi[0].data[:,0,0])])
        n_panels = new_py_zmax - new_py_zmin
        chans_ghz = chan2freq(np.array(range(source['z_min']-1-boundary, source['z_max']-1+boundary, 1)), 
                              infile).to(u.GHz)
        chan_per_page = xsize * ysize
        page = 0

        # Save channel maps to a multi-page pdf
        with PdfPages(outfile) as pdf:
            for c, c_ghz in zip(range(new_py_zmin,new_py_zmax,1), chans_ghz):
                chan_im = hdulist_hi[0].data[c, :, :]
                mask_im = hdulist_mask[0].data[c, :, :]
                chan_num = c - new_py_zmin    # indexed from 0
                panel = int((chan_num) % chan_per_page)
                column = panel % xsize
                row = int(((panel - column) / xsize) % chan_per_page)
                if panel == 0:
                    page += 1
                    logger.warning('\tPlotting page {} of channel maps'.format(page))
                    if n_panels - chan_num < chan_per_page:
                        if (n_panels - chan_num) % xsize == 0:
                            ysize = int((n_panels - (chan_num)) / xsize)
                        else:
                            ysize = int((n_panels - (chan_num)) / xsize) + 1
                    fig = plt.figure(figsize=(xsize*3, ysize*3))
                    gs = fig.add_gridspec(ysize, xsize, hspace=0, wspace=0)
                    axs = gs.subplots(sharex='col', sharey='row', subplot_kw=dict(projection=owcs))
                    # Some brute force relocating of labels because tight_layout screws up panel spacing:
                    fig.suptitle(source['name'] + id_label, fontsize=24, y=0.91 + (0.09-0.09*(ysize/5))*1.3)
                    fig.supxlabel(x_label, fontsize=22, y=0.05 - (0.05 * (5/ysize)-0.05)*1.4)
                    fig.supylabel(y_label, fontsize=22, x=-0.01 - (0.01 * (4/xsize)-0.01)*12)

                if axs.ndim == 1:
                    axs = axs[None,:]
                axs[row,column].imshow(chan_im, cmap=pvd_map, norm=divnorm, origin='lower', 
                                       transform=axs[row,column].get_transform(cubew))
                axs[row,column].text(0.5, 0.5, '{:.4f} GHz'.format(c_ghz.value), fontsize=14)
                axs[row,column].tick_params(axis='both', which='major', labelsize=14, length=6, width=2)
                axs[row,column].set(facecolor="white")  # Doesn't work with the color im
                # Plot positive contours
                try:
                    axs[row,column].contour(chan_im, cmap='Oranges_r', linewidths=1.2, 
                                            levels=base_contour * 2 ** np.arange(10),
                                            transform=axs[row,column].get_transform(cubew))
                except:
                    pass
                # Plot negative contours
                try:
                    axs[row,column].contour(chan_im, cmap='YlOrBr_r', linewidths=1.2, linestyles='dashed',
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

                if ((row+1) * (column+1) == chan_per_page) | (chan_num == n_panels - 1):
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
