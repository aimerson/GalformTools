#! /usr/bin/env python

from setuptools import setup, find_packages

datafiles = ['data/Simulations/*.xml','data/stellarAstrophysics/Vega/*.xml',\
                 'data/dust/*/*.xml','data/powerSpectra/*.dat']

setup(name='galform',
      version='0.1',
      description='User tools for the Galform semi-annalytical model',
      url='http://bitbucket.org/aimerson/galformtools',
      author='Alex Merson',
      author_email='alex.i.merson@gmail.com',
      license='MIT',
      packages=find_packages(),
      package_dir={'galform':'galform'},
      package_data={'galform':datafiles},
      zip_safe=False)

