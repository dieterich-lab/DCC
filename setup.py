"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from codecs import open
from os import path
from setuptools import setup

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst')) as f:
    long_description = f.read()

setup(
    name='DCC',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.5.0',

    description='Detect circRNAs from chimeras',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/dieterich-lab/DCC',

    # Author details

    maintainer='Tobias Jakobi',
    maintainer_email='Tobias.Jakobi@med.Uni-Heidelberg.DE',

    author='Jun Cheng',
    author_email='s6juncheng@gmail.com',


    # Choose your license
    license='License :: OSI Approved :: GNU General Public License (GPL)',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License (GPL)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # 'Programming Language :: Python :: 2',
        # 'Programming Language :: Python :: 2.6',
        # 'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.2',
        # 'Programming Language :: Python :: 3.3',
        # 'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',

    ],

    # What does your project relate to?
    keywords='circRNA detection and quantification',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=['DCC'],

    # setup_requires=['Cython','pysam','matplotlib'],

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
     install_requires=[
         'HTSeq',
    #     'pysam >= 0.13',
    #     'numpy',
    #     'pandas',
    #     'Cython'
    ],

    #install_requires=read('requirements.txt').splitlines(),

    # python_requires='<3',

    # List additional groups of dependencies here (e.g. development
    # dependencies). You can install these using the following syntax,
    # for example:
    # $ pip install -e .[dev,test]
    # extras_require={
    #     'dev': ['check-manifest'],
    #     'test': ['coverage'],
    # },

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    package_data={
        'DCC': ['data/DCC.Repeats'],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('my_data', ['data/data_file'])],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'DCC=DCC:main'
        ],
    },
    scripts=[
        'scripts/DCC',
    ],

    project_urls={  # Optional
        'Bug Reports': 'https://github.com/dieterich-lab/DCC/issues',
        'Dieterich Lab': 'https://dieterichlab.org',
        'Source': 'https://github.com/dieterich-lab/DCC',
        'Documentation': 'http://docs.circ.tools'
},
)
