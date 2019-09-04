import setuptools
from setuptools import Extension

with open("README.md", "r") as fh:
    long_description = fh.read()

ext = Extension(
    'romspy/interpolation/vertical/linear',
    sources=['romspy/interpolation/vertical/linear.c'],
    extra_compile_args=['-fopenmp'],
    extra_link_args=['-lgomp'])

setuptools.setup(
    name="romspy",
    version="1.0.dev1",
    author="Nicolas Munnich",
    description="Preprocessing files for use in ROMS",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    ext_modules=[ext],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Programming Language :: C",
        "Topic :: Scientific/Engineering",
        "Typing :: Typed"
    ],
    python_requires='>=3.7',
)
