from setuptools import setup, find_packages

setup(
    name='Correct_TDS',
    version='0.1.13',    
    description='Correct TDS software for statistical gap analysis',
    url='https://github.com/THzbiophotonics/Correct-TDS-in-development/tree/master',
    author='PERETTI Romain',
    author_email='"romain.peretti@univ-lille.fr',
    license='IEMN',
    packages=find_packages(),
    python_requires='>=3.9, <3.13',
    install_requires=[
        'mpi4py>=2.0',
        'numpy',  
        'PyQt5>=5.15.0',
        'matplotlib>=3.5.0',
        'h5py>=3.6.0',
        'scipy>=1.8.0',
        'scikit-learn>=1.0.0',
        'numba>=0.57.0, <0.61.0',
        'pyswarm'
    ],
    entry_points={
        'console_scripts': [
            'Correct-TDS-interface = Correct_TDS.interface:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
)