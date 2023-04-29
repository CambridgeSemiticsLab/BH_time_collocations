import setuptools
from distutils.core import setup
from pathlib import Path

REQUIREMENTS_FILE = 'requirements.txt'
requirements = Path(REQUIREMENTS_FILE).read_text().splitlines()

setup(
    name='kingham_thesis',
    version='1.0',
    description='Analysis package for Biblical Hebrew Time Collostructions',
    author='Cody Kingham',
    license='MIT',
    url='https://github.com/CambridgeSemiticsLab/BH_time_collocations',
    package_dir={'': 'src'},
    packages=setuptools.find_packages(include=['kingham_thesis']),
    install_requires=requirements,
)

