import codecs
import os.path
from setuptools import setup

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name='Utmos',
    version=get_version('utmos/__init__.py'),
    author="ACEnglish",
    author_email="acenglish@gmail.com",
    url="https://github.com/acenglish/utmos",
    packages=['utmos'],
    license='MIT',
    description="Maximum-coverage algorithm for genomic variants",
    long_description=open('README.md', encoding='UTF-8').read(),
    long_description_content_type='text/markdown',
    entry_points={
      'console_scripts': [
         'utmos = utmos.__main__:main'
      ]
    },
    install_requires=[
        "truvari>=3.5.0",
        "scikit-allel==1.3.5",
        "h5py==3.7.0",
    ],
)
