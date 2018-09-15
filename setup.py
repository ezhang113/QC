#!/usr/qc/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='qc',
    version='0.0.2',
    url='github.com/byee4/QC',
    license='',
    author='gpratt, byee4',
    author_email='bay001@ucsd.edu',
    description='Gathers metrics files from eCLIP pipelines and produces a tabbed QC file.',
    packages=['qc'],
    include_package_data=True,
    package_dir={
        'qc': 'qc',
    },
    entry_points = {
        'console_scripts': [
            'eclip_qc = qc.qcsummary_eclip:main',
            'rnaseq_qc = qc.qcsummary_rnaseq:main',
        ]
    }
)
