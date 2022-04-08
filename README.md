SoFiA Image Pipeline (SIP)
=====

Introduction
------------
SIP takes a SoFiA generated source catalog and produce images for publication or quick inspection.  Images include moments, pv-diagrams, optical overlay with HI contours and spectra, with and without noise.

Usage
-----
Works under the assumption that the catalog file and cubelets are in the same directory structure as created by SoFiA. Output from this program will be in a folder next to `*_cubelets` called `*_figures`.

Get help on the parameters: 
* `python3 image_pipeline.py -h`

Examples:
* Use the xml file with output images in pdf format: \
`python3 image_pipeline.py -c <path/to/catalog.xml> -x pdf`

* Use ascii file with output images in png format and specify original data set: \
`python3 image_pipeline.py -c <path/to/catalog.txt> -x png -o <path/to/original_cube.fits>`

The program outputs individual *pngs but it may be useful to combine them all into one.  This can be accomplished with something like `imagemagick`.  Here is an example:
```
convert_im = "/usr/local/Cellar/imagemagick/7.1.0-13/bin/convert"
for src in sources:
    print("  {}".format(src))
    new_png = "{}_combo.png".format(src)
    survey = "dss2blue"
    os.system("{} {}_mom0_{}.png {}_mom0.png {}_snr.png {}_mom1.png +append temp.png".format(convert_im, src, survey, src, src, src, src))
    os.system("{} {}_spec.png -resize 125% temp2.png".format(convert_im, src))
    os.system("{} {}_specfull.png -resize 125% temp3.png".format(convert_im, src))
    os.system("{} temp2.png temp3.png {}_pv.png +append temp4.png".format(convert_im, src, src, src))
    os.system("{} temp.png temp4.png -append {}".format(convert_im, new_png))
    os.system('rm -rf temp.png temp2.png temp3.png temp4.png')
```
On your system, you will likely be able to replace `/usr/local/Cellar/imagemagick/7.1.0-13/bin/convert` with simply `convert`.  Note `+` or `-` in front of `append` controls if the images are combined horizontally or vertically.
This example places all images (except the PanSTARRS overlay) on two rows, with the spatial plots on the top row and the spectral plots on the bottom row. 

Requirements
------------
Tested with: \
`Python 3.9.7` \
`Astropy 5.0.1 (and 4.2.1)` \
`Astroquery 0.4.5 (and 0.4.1)` \
`Numpy 1.20.2` \
`Reproject 0.7.1` \
`pvextractor 0.2`

Only guaranteed to run with the SoFiA 2.4.1 (bba8c43) or later.  SoFiA-2 needs to have `parameter.wcs = True` for units in images to make sense. To be on the safe side, also set `parameter.physical = True`. Also, if you run SoFiA-2 on a subregion, must have `parameter.offset = True`.


Installation
------------
Copy all the files to your favorite directory or do a git clone (nothing fancy supported yet.)

Problems
--------
Probably.

Documentation
-------------
Haha. You're looking at it.

Version history
---------------
* SIP 0.1.0
    * Release TBD
    * Still under development 
    
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
