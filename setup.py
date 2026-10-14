# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 18:09:24 2024

@author: Nate Whitsett
"""

from setuptools import setup, find_packages

setup(
   name='ardor',
   version='1.0',
   url='https://github.com/AstroMusers/ardor.git',
   description='A flare pipeline designed to search for induced flares in exoplanet host stars.',
   author='Nathan Whitsett',
   author_email='whitsett.n@wustl.edu',
   license='MIT',
   packages=find_packages(where='src'),
   package_dir={'': 'src'},
   install_requires=[
       'numpy',
       'scipy',
       'matplotlib',
       'pandas',
       'astropy',
   ],
   extras_require={
      'pcat': ['pcat @ git+https://github.com/tansudaylan/pcat.git'],
      'allesfitter': ['allesfitter'],
      'roman': ['tdpy @ git+https://github.com/tansudaylan/tdpy.git', 'emcee', 'dynesty'],
   },
)
