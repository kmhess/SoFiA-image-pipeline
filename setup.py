from setuptools import setup

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='sofia_image_pipeline',
    version='0.0.1',
    description='SIP takes a SoFiA generated source catalog and produce images for publication or quick inspection. Images include HI contours overlaid on multiwavelength images, HI moment maps, pixel-by-pixel SNR maps, pv-diagrams with SoFiA mask, and spectra with and without noise.',
    url='https://github.com/kmhess/SoFiA-image-pipeline',
    author='Kelley Hess',
    author_email='kmhess.astro@gmail.com',
    packages=[
        "src",
        "src.modules",
    ],
    install_requires=required,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python :: 3',
    ],
    entry_points={
        'console_scripts': [
            'sofia_image_pipeline=src.image_pipeline:main',
        ],
    },
)
