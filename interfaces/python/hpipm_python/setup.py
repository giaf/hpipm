from setuptools import setup, find_packages

import hpipm_python

setup(name='hpipm-python',
   version='0.1',
   description='Python interface to HPIPM',
   url='http://github.com/giaf/hpipm',
   author='Andrea Zanelli - Gianluca Frison',
   license='LGPL',
   packages = find_packages(),
   zip_safe=False)
