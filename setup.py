from distutils.core import setup

setup(
    name='husky',
    description='huskython: Husky in Python ',
    version='0.1-dev',
    url='https://github.com/NicoRenaud/huskython',
    package_dir = {
    'husky' : 'husky',
    'husky.hamiltonian' : 'husky/hamiltonian',
    'husky.hamiltonian.huckel' : 'husky/hamiltonian/huckel',
    'husky.hamiltonian.model' : 'husky/hamiltonian/model',
    'husky.transport' : 'husky/transport'
    },
    packages=['husky',
              'husky.hamiltonian','husky.hamiltonian.huckel','husky.hamiltonian.model',
    		  'husky.transport']
)
