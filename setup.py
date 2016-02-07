#!/usr/bin/python

from setuptools import setup
import sys, site, os

setup(name='fasticlip',
      version='0.9.3',
      description='Fully Automated and Standardized iCLIP (FAST-iCLIP) is a fully automated tool to process iCLIP data.',
      
	  #install_requires = ['numpy>=1.7', 'pandas>=0.14', 'matplotlib_venn>=0.11', 'matplotlib>=1.3.1'],
	  install_requires = ['numpy>=1.7', 'pandas>=0.14', 'matplotlib>=1.3.1'],

	  packages=['fasticlip', 'fasticlip'],
	  package_dir={'fasticlip': 'fasticlip'},
	  
	  entry_points = {'console_scripts': [
	  'fasticlip = fasticlip.fasticlip:main']},
	  
	  author='Brian Do',
      author_email='bdo@stanford.edu',
      url='https://www.github.com/ChangLab/FAST-iCLIP',
	  license = "GPL2",
     )
