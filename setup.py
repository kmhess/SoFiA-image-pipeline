from setuptools import setup

setup(
    name='sofia_image_pipeline',
    version='1.3.10',
    description='SIP takes a SoFiA generated source catalog and produce images for publication or quick inspection. Images include HI contours overlaid on multiwavelength images, HI moment maps, pixel-by-pixel SNR maps, pv-diagrams with SoFiA mask, and spectra with and without noise.',
    url='https://github.com/kmhess/SoFiA-image-pipeline',
    author='Kelley M. Hess',
    author_email='kmhess.astro@gmail.com',
    packages=[
        "src",
        "src.modules",
    ],
    install_requires=[
        "astropy >= 5.0.2",
        "astroquery >= 0.4.7",
        "matplotlib >= 3.5.1",
        "numpy >= 1.22.0",
        "Pillow >= 10.0.1",
        "pvextractor >= 0.4",
        "requests >= 2.31.0",
        "setuptools >= 65.5.1",
        "sphinx_rtd_theme >= 1.0.0",
        "xmltodict >= 0.12.0"
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
    entry_points={
        'console_scripts': [
            'sofia_image_pipeline=src.image_pipeline:main',
            'download_usr_fig=src.download_usr_fig:main'
        ],
    },
)
