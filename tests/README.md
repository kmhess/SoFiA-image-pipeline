# Tests

We provide a test case so that users are able to verify the installation of SIP has worked as expected.

### 1. Setup

The first step to run the test case is to generate a catalogue file. This will require a local installation of [`sofia`](https://github.com/SoFiA-Admin/SoFiA-2) (installation instructions available on in the repository. Once installed, it can be run on the provided image cube (`UGC7012.fits`) from within the `tests/data` subdirectory. The command to run this, from `tests/data` is:

```
sofia sofia2_ugc7012.par
```

### 2. Install CLI

We provide a CLI tool for running the `sofia_image_pipeline`. Once you have cloned this repository you can install it in a development environment by running

```
python3 setup.py develop
```

from the base directory of the repository. You can verify the installation was successful by running

```
sofia_image_pipeline
```

which should then print usage and configuration instructions to your screen

```
usage: sofia_image_pipeline [-h] -c CATALOG [-x SUFFIX] [-o ORIGINAL] [-b BEAM] [-i IMAGE_SIZE] [-snr SNR_RANGE SNR_RANGE] [-s [SURVEYS [SURVEYS ...]]] [-m [IMAGEMAGICK]] [-ui USER_IMAGE] [-ur USER_RANGE USER_RANGE]
sofia_image_pipeline: error: the following arguments are required: -c/--catalog
```

### 3. Run on test case

The final step is to run the pipeline on the test data. This can be done from the base directory with

```
sofia_image_pipeline -c tests/data/UGC7012_cat.xml
```

The output will indicate whether or not it ran successfully. You should be able to see the output products of the pipeline in the new folder at `tests/data/UGC7012_figures`.
