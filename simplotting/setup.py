# setup.py

from setuptools import setup, find_packages

setup(
    name='simplotting',
    version='0.1',
    description='Package to read and plot outputs from the McGill-SIM using the mu-phi rheology',
    author='Félix St-Denis',
    author_email='felix.st-denis@mail.mcgill.ca',
    packages=find_packages(include=['simplotting', 'simplotting.*']),
    install_requires=[
        # Add any dependencies here
    ],
    entry_points={
        'console_scripts': [
            # Define any command line scripts here
        ],
    },
)