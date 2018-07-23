from setuptools import setup, find_packages
from pkgutil import walk_packages

import hpipm_python

# def find_packages(path=__path__, prefix=""):
#     yield prefix
#     prefix = prefix + "."
#     for _, name, ispkg in walk_packages(path, prefix):
#         if ispkg:
#             yield name

setup(name='hpipm-python',
   version='0.1',
   description='Python interface to HPIPM',
   url='http://github.com/giaf/hpipm',
   author='Andrea Zanelli - Gianluca Frison',
   license='LGPL',
   packages = find_packages(),
   zip_safe=False)
