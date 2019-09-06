import setuptools
import os
import subprocess
# from setuptools import Extension
from setuptools.command.develop import develop
from setuptools.command.install import install


class PostDevelopCommand(develop):
    """Post-installation for development mode."""

    def run(self):
        # PUT YOUR POST-INSTALL SCRIPT HERE or CALL A FUNCTION
        develop.run(self)


class PostInstallCommand(install):
    """Post-installation for installation mode."""

    def run(self):
        filepath = os.path.realpath(__file__)
        vertical_path = os.path.join(filepath[0], "romspy/interpolation/vertical/")
        subprocess.run(["gcc", "-Wall", "-Wextra", "-fPIC", "-fopenmp", "-shared", vertical_path + "linear.c", "-o",
                        vertical_path + "linear.so"])
        install.run(self)


with open("README.md", "r") as fh:
    long_description = fh.read()

# ext = Extension(
#     'romspy/interpolation/vertical/linear',
#     sources=['romspy/interpolation/vertical/linear.c'],
#     extra_compile_args=['-fopenmp'],
#     extra_link_args=['-lgomp'])

setuptools.setup(
    name="romspy",
    version="1.0.dev3",
    author="Nicolas Munnich",
    author_email="nmdm20@bath.ac.uk",
    description="Preprocessing files for use in ROMS",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://https://github.com/Saixos/romspy",
    packages=setuptools.find_packages(include=["romspy", "romspy.*"]),
    install_requires=[
        'numpy>=1.17.1',
        'xarray>=0.12.3',
        'cdo>=1.4.3',
        'netCDF4>=1.5.1.2',
    ],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Programming Language :: C",
        "Topic :: Scientific/Engineering",
        "Typing :: Typed"
    ],
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
    },
    keywords='preprocessing ROMS roms Regional Ocean Modelling System cdo Climate Data Operators netcdf4 ',
    python_requires='>=3.6',
)
