from setuptools import setup
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info
import subprocess as sp
import platform

def make_hkl():
    '''
    Compile the C huckel code and turns it in a dynamic library
    '''
    print('--> Compile huckel library')
    src_path = './husky/hamiltonian/huckel/src/'
    osname = platform.system()
    if osname == 'Linux':
        cmd = 'gcc -shared -o ../hkl.so -fPIC huckel.c'
    elif osname == 'Darwin':
        cmd = 'gcc -dynamiclib huckel.c -o ../hkl.so'
    else:
        raise ValueError('Environement %s not supported.' %osname)

    sp.check_call(cmd,cwd=src_path,shell=True)

class hklinstall(install):
    '''
    custom hadler for the install command
    '''
    def run(self):
        make_hkl()
        super().run()

class hkldevelop(develop):
    '''
    custom hadler for the install command
    '''
    def run(self):
        make_hkl()
        super().run()

class hklegg(egg_info):
    '''
    custom hadler for the install command
    '''
    def run(self):
        make_hkl()
        super().run()

setup(
    name='husky',
    description='huskython: Husky in Python ',
    version='0.1-dev',
    url='https://github.com/NicoRenaud/huskython',
    packages=['husky'],
    install_requires=[
    'numpy>=1.13.3',
    'scipy>=1.0.0'
    ],
    extras_require={
    'test':['nose','coverage']
    },
    cmdclass={'install' : hklinstall, 
              'develop' : hkldevelop,
              'egg_info': hklegg }
)
