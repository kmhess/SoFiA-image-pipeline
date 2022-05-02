# Tests

We provide a test case so that users are able to verify the installation of SIP has worked as expected. First step to run the tests is to set up the catalogue data. Next, we run the end-to-end test.

## 1. Setup

Run [`sofia`](https://github.com/SoFiA-Admin/SoFiA-2) on the image cube provided (`UGC7012.fits`) from within the `tests/data` subdirectory. This will generate the `UGC7012_cat.xml` file required for running SIP. Once `sofia` is installed run:

```
sofia sofia2_ugc7012.par
```

## 2. Run test

1. Run `image_pipeline.py` as Python script on the catalogue. To do this:

```
python3 src/image_pipeline.py -c tests/data/UGC7012_cat.xml
```

2. Use unittest framework to run the test. Once the catalogue has been generated, you can run:

```
python3 -m unittest tests/test_image_pipeline.py 
```