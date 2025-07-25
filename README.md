SoFiA Image Pipeline (SIP)
=====
[![DOI](https://zenodo.org/badge/455147174.svg)](https://zenodo.org/badge/latestdoi/455147174)
[![Docker build latest](https://github.com/axshen/SoFiA-image-pipeline/actions/workflows/docker-build-latest.yml/badge.svg)](https://github.com/axshen/SoFiA-image-pipeline/actions/workflows/docker-build-latest.yml)

Introduction
------------
SIP takes a SoFiA generated source catalog and produce images for publication or quick inspection.  Images include spectral line total intensity contours overlaid on multiwavelength images, spectral line moment maps, pixel-by-pixel SNR maps, pv-diagrams with source mask, masked spectra, and aperture spectra using the 2D projection of the source mask.  SIP also generates an image of the aperture and masked spectrum overplotted.

![](docs/sofia_test_output_3_combo.png)
<!-- ![](docs/sofia_test_output_3_specboth.png) -->
<img src="docs/sofia_test_output_3_specboth.png" alt="drawing" width="35%"/>

Requirements
------------
This code has been developed and tested with Python 3.9.7, and appears to work with up to Python 3.13.

Combining individual images with the `-m` option requires [ImageMagick](https://imagemagick.org) be installed.

SIP was written with the outputs from SoFiA 2.4.1 (bba8c43) and later. Through troubleshooting, we have improved the output of both SIP and [SoFiA-2](https://github.com/SoFiA-Admin/SoFiA-2), so for best performance, please use the latest version of each. 

For the best SIP unit conversion performance, SoFiA-2 needs to have been run with the following parameters: 

```
parameter.wcs = True
parameter.physical = True
parameter.offset = True
```

Installation
------------

### Setting up an environment

We recommend using a [virtual environment](https://docs.python.org/3/library/venv.html) to install and rup SIP so as to not conflict with other packages.  For example:
```
python3 -m venv .venv
source .venv/bin/activate
```

### PyPI

You can install the latest SIP released on [PyPI](https://pypi.org/project/sofia-image-pipeline/) by running pip install:

```
pip install sofia-image-pipeline
```

### Development

You can install the latest GitHub version of SIP locally by cloning the repository and running:

```
python3 setup.py develop
```

### Docker

We provide a Docker image for use in containerised applications. This image can be found [here](https://hub.docker.com/r/sofiapipeline/image_pipeline).
To use the latest official Docker image:

```
docker run -it -v <cwd>/<folder>:/app/<folder> sofiapipeline/image_pipeline:latest -c <folder>/<catalog.xml>
```

Usage
------------
SIP works under the assumption that the user has run [SoFiA-2](https://github.com/SoFiA-Admin/SoFiA-2) which generated an xml and/or ascii catalog file, and fits moment maps and cubelets for each source.  SIP assumes that these files are in the same directory structure as created by SoFiA-2 where the catalog file and `*_cubelets/` folder are in the same directory.  The output from SIP will be in a newly created folder next to `*_cubelets/` called `*_figures/`.

```
$ sofia_image_pipeline

usage: sofia_image_pipeline [-h] -c CATALOG [-id [SOURCE_ID ...]] [-s [SURVEYS ...]] [-ui USER_IMAGE] [-ur USER_RANGE USER_RANGE] [-line SPECTRAL_LINE] [-i IMAGE_SIZE] [-snr SNR_RANGE SNR_RANGE] [-o ORIGINAL] [-b BEAM] [-cw CHAN_WIDTH] [-x SUFFIX] [-m [IMAGEMAGICK]]
sofia_image_pipeline: error: the following arguments are required: -c/--catalog
```

### Help

To get help on the parameters: 

```
sofia_image_pipeline -h
```

### Options

```
REQUIRED:
    -c     Catalog file. Can be the ascii file ending in .txt or the XML file from SoFiA-2.
    
OPTIONAL:
    -id    Select certain sources, or range of sources from catalog based on the source id number. Default is all sources. Include `0` to make summary image of all sources.
    -s     List of surveys on which to overlay HI contours. Only the first entry will be used in the combined image if `-m` option is used. If 'none', work in offline mode. Default is 'DSS2 Blue'.
    -ui    User supplied image for overlaying HI contours.  Can use this in combination with `-s` and a list of surveys.
    -ur    Percentile range when plotting the user supplied image.  Requires two values. Default is [10., 99.].
    -line  Specify a spectral line for all sources in catalog. Default is 'HI'.  Also possible is 'CO(1-0)' up to (3-2) and 'OH_1667' and the 3 other OH lines which may fall within L-band observations.
    -i     Minimum image size (ARCMIN). Images will be square. If an HI source exceeds the requested size, a larger image to fit the HI contours will be generated. Default is 6 arcmin.
    -snr   Specify the SNR range within which to plot the lowest HI contour. Requires 2 values. Default is [2.0, 3.0].
    -o     Path to the original data file on which source finding was conducted. This allows the spectrum with noise to be plotted over the full spectral range of the original cube.  
    -b     Synthesized beam dimensions. If the primary header of the FITS files do not contain the beam information, this can be provided by the user. Accepts 1 to 3 values in order (bmaj,bmin,bpa). Format is comma separated, with no spaces.
    -cw    Channel width. This is only necessary if a source cubelet is not available, for example if user only has a mom0.  Channel width must then be provided in the native units of the original data cube (typically Hz or m/s.)
    -x     Output image file type. Any file type accepted by plt.savefig() is in theory valid.  Default is 'png'.
    -m     Make a combined image using ImageMagick.  If a path is provided after this option, it is assumed to be the path to the `convert` executable of ImageMagick. 
```

### Examples

* Use the xml catalog file to output images in default png format:

```
sofia_image_pipeline -c <path/to/catalog.xml>
```

* Use ascii catalog file with output images in pdf format and specify original data set to plot full noise spectrum: 

```
sofia_image_pipeline -c <path/to/catalog.txt> -x pdf -o <path/to/original_cube.fits>
```

* Request spectral line contours on multiple survey images, separated by a space, and make a combined image for each source:

```
sofia_image_pipeline -c <path/to/catalog.txt> -s panstarrs 'GALEX Far UV' -m
```

Test data cube
--------
The SoFiA test data cube can be found through the SoFiA-2 wiki [here](https://github.com/SoFiA-Admin/SoFiA-2/wiki#test-data-cube) (14.2 MB).
This cube contains HI emission from several galaxies around NGC 4036.

The current version of this repo also has a test data set (without instructions) in the `tests/data/` folder.  The cube contains HI emission from several galaxies around UGC 7012.

Advanced tips
--------
* Some troubleshooting is available on the [wiki](https://github.com/kmhess/SoFiA-image-pipeline/wiki).

* If you have a large catalog of sources, start by testing SIP with the `-id N` option, where `N` is a source id number.  Make sure the image and text outputs from SIP for that source are as you expect.  Adjust optional variables as necessary.  Run on your larger catalog.

* SIP is now capable of doing spectral lines other than HI.  So far `HI`, `CO(1-0)`, `CO(2-1)`, `CO(3-2)`, and `OH_1667`, `OH_1665`, `OH_1720`, `OH_1612` are the only allowed options.  This is work in progress.

* SIP can be run in 'offline' mode, by setting `-s none`.  In this case no survey archive data will be downloaded, and SIP will only generate the HI images.  Any surveys in a list in which `none` appears will be ignored.

* SIP now generates a plot called `*specboth.png` which overlays the masked and aperture spectra on the same plot, although it is not in the combo plot.

* When saving files as postscript, use `-x eps` to maintain the figure dimensions.

* Available surveys from `astroquery` can be found by running:
```
from astroquery.skyview import SkyView
SkyView.list_surveys()
# or 
SkyView.survey_dict
```

* In addition to overlaying HI contours on survey images available through `astroquery`, a user can request WISE images (`'WISE W#'` where `#` is the band number); false color images from `'decals'` (DR10), `'decals-dr9'` (DR9), `'panstarrs'`, `'sdss'`, or `'decaps'`; or gray-scale HST-ACS Mosaic images for sources within the COSMOS field with `'hst'`.  The HST image size is currently hardcoded to 40 arcsec on a side. 

* Downloading survey images from `astroquery.SkyView` or other online sources is the greatest limiting factor in the speed of SIP.  To avoid this, for catalogs with a high source density, you may consider downloading one large image to disk before running SIP.  For this purpose, we have included the command line tool `download_usr_fig`. For example:
```
download_usr_fig -h
download_usr_fig -ra 174.465 -dec 21.9743 -s 'Survey Name' -i 0.5 -o my_image
sofia_image_pipeline -c <path/to/catalog.xml> -ui <my_image_SurveyName.fits>
```

* SIP always outputs individual figures for each SoFiA FITS file.  If you did not produce a combined summary image with the `-m` option, you can still create it without re-running SIP if you have ImageMagick installed.  Here is example python code using a for loop over the sources to execute terminal commands.
```
convert_im = "/usr/local/Cellar/imagemagick/7.1.0-13/bin/convert"
for src in sources:
    print("  {}".format(src))
    new_png = "{}_combo.png".format(src)
    survey = "dss2blue"
    os.system("{0} {1}_mom0_{2}.png {1}_mom0.png {1}_snr.png {1}_mom1.png {1}_mom2.png +append temp.png".format(convert_im, src, survey))
    os.system("{} {}_spec.png -resize 125% temp2.png".format(convert_im, src))
    os.system("{} {}_specfull.png -resize 125% temp3.png".format(convert_im, src))
    os.system("{0} temp2.png temp3.png {1}_pv.png {1}_pv_min.png +append temp4.png".format(convert_im, src))
    os.system("{} temp.png temp4.png -append {}".format(convert_im, new_png))
    os.system('rm -rf temp.png temp2.png temp3.png temp4.png')
```
On your system, you will likely be able to replace `/usr/local/Cellar/imagemagick/7.1.0-13/bin/convert` with simply `convert`.  Note `+` or `-` in front of `append` controls if the images are combined horizontally or vertically.
This example places all images on two rows, with the spatial plots on the top row and the spectral plots on the bottom row. 


Known Issues
--------
See the github repo for known bugs and desired enhancements.  We aim to fix serious bugs as quickly as possible.

In addition we are aware of the following issues:
* `download_usr_fig` can download full color images from PanSTARRS and DECaLS, but these can not yet be read as user supplied input to `sofia_image_pipeline`.
* WISE images, PanSTARRS, DECaLS, and DECaPS cannot (yet?) be plotted with Galactic coordinates.
* For data with channels that are not uniform in width (e.g. `SPECSYS = FELO-OPT`), SIP's conversion to km/s is off compared to SoFiA-2's: the programs use formula from [here](https://www.astro.rug.nl/software/kapteyn/spectralbackground.html#aips-axis-type-felo) or use wcslib to do the conversion, respectively.  We haven't tracked down the discrepancy.  To the best of our knowledge, only relatively old radio data observing nearby galaxies, might be in this `FELO` format. 
* No exceptions are caught for `socket.timeout` during downloads (seen for `-s panstarrs` when image requested was 20.5 arcmin). We've noticed that simply rerunning the request at a later time has solved the issue.

### Other caveats
* SIP is not guaranteed to work with output from SoFiA-1. It may produce images, but these may not have the correct units, reference frame, or SIP may crash entirely. Because SoFiA-1 is no longer being developed, there are no plans to make SIP work around these issues.
* The SIP code and documentation is still a work in progress.  Thanks for your patience and suggestions for improvement.

Version history
---------------
* SIP 1.3.5
    * Released 31 October 2024
* SIP 1.3.1
    * Released 13 June 2024
* SIP 1.3.0
    * Released 7 June 2024 
* SIP 1.2.0
    * Released 12 July 2022 
* SIP 1.1.0
    * Released 9 May 2022
* SIP 1.0.0
    * Released 18 April 2022
    
Copyright and licence
---------------------
SIP was created by Kelley M. Hess

© 2022 Kelley M. Hess

This programme is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free 
Software Foundation, either version 3 of the License, or (at your option) any 
later version.

This programme is distributed in the hope that it will be useful, but **without 
any warranty**; without even the implied warranty of **merchantability** or **fitness 
for a particular purpose**. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
this programme. If not, see http://www.gnu.org/licenses/.
