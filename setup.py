#*****************************************************************************
#  Copyright (C) 2017 Lee Worden <worden dot lee at gmail dot com>
#
#  Distributed under the terms of the GNU General Public License (GPL) v.2
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# http://python-packaging.readthedocs.org/en/latest/minimal.html
from setuptools import setup

setup(name='boxmodel',
      version='0.1',
      description='Sage classes for compartmental ("box") models',
      url='http://github.com/tcporco/SageBoxModels',
      author='Lee Worden',
      author_email='worden.lee@gmail.com',
      license='GPL 2.0',
      packages=['boxmodel'],
      zip_safe=False)
