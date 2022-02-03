SoFiA Image Pipeline (SIP)
=====

Introduction
------------
Take a SoFiA generated source catalog and produce images for publication or quick inspection.  Images include moments, pv-diagrams, optical overlay with HI contours and spectra, with and without noise.

Usage
-----
Get help on the parameters: 
* `python3 image_pipeline.py -h`

Examples:
* Use the xml file with output images in pdf format: \
`python3 image_pipeline.py -c <path/to/catalog.xml> -s pdf`

* Use ascii file with output images in png format and specify original data set: \
`python3 image_pipeline.py -c <path/to/catalog.txt> -s png -o <path/to/original_cube.fits>`

Requirements
------------
Tested with: \
`Python 3.9.7` \
`Astropy 5.0.1 (and 4.2.1)` \
`Astroquery 0.4.5 (and 0.4.1)` \
`Numpy 1.20.2` \
`Reproject 0.7.1`

Only guaranteed to run with the SoFiA 2.4.1 (bba8c43) for now, but may work with earlier verions of SoFiA-2.  SoFiA-2 needs to have `parameter.wcs` = True for units in images to make sense. 

Does not work with current version of SoFiA-1 due to some discrepancies in output.  


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