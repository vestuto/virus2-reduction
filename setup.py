#!/usr/bin/env python
"""
Purpose: enable standard pip install into standard sys.path() search path list
"""
from setuptools import setup, find_packages

setup(
    name='v2reduce',
    description='VIRUS2 Software for Data Reduction',
    version='2025.6.16',
    author='Jason Vestuto',
    author_email='jason.vestuto@austin.utexas.edu',
    packages=find_packages(include=['v2reduce', 'v2reduce.*']),
    package_data={'v2reduce': ['data/*.txt'],},
    include_package_data=True,
    scripts=['scripts/v2reduce_process.py',],
    install_requires=[
        'astropy',
        'numpy',
        'matplotlib',
        'pandas',
        'scikit-learn',
        'scipy',
        'seaborn'
    ],
    python_requires='>=3.11',
)

