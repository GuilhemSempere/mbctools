# Use this file to generate and publish a PyPi package for each new release
#
# Procedure:
#   - make sure __version__ variable is correctly set in mbctools.py, and that GitHub tag and release exist for the same version number
#   - enter command "python3 setup.py sdist bdist_wheel" to build the package
#   - enter command "twine upload dist/*", which will prompt for your PyPi token


import re
from setuptools import setup, find_packages

version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('mbctools.py').read(),
    re.M
).group(1)

with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")

setup(
    name='mbctools',
    version=version,    
    description='mbctools: A User-Friendly Metabarcoding and Cross-Platform Pipeline for Analyzing Multiple Amplicon Sequencing Data across a Large Diversity of Organisms',
    long_description = long_descr,
    long_description_content_type="text/markdown",
    url='https://github.com/GuilhemSempere/mbctools/blob/dev/mbctools.py',
    author='Christian Barnabe, Guilhem Sempere',
    author_email='guilhem.sempere@cirad.fr',
    license='MIT',
    packages=find_packages(),
    install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    py_modules=['mbctools'],
    entry_points={
        "console_scripts": [
            'mbctools = mbctools:main'
        ]
    },
)