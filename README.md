# Correct-TDS




## With the pip install package :

- Install pyOpt from https://github.com/madebr/pyOpt (optional)

- Install anaconda (https://docs.anaconda.com/anaconda/install/index.html)

- Create an environment

- Open the environment terminal

- To install the package, write in the command line :

	```
	pip install Correct-TDS
	```
- To run Correct@TDS, write in the command line :

	```
 	Correct-TDS-interface
	```
- The gui can freeze while computing, it's normal, it's still in development

## Manual instalation

- Install anaconda (https://docs.anaconda.com/anaconda/install/index.html)

- Install mpi4py with the command

	```
	conda install -c conda-forge mpi4py
	```
- Install pyswarm with the command

	```
	pip install pyswarm
	```
	
- Install numba with the command

	```
	pip install numba
	```	
- Install pyOpt from https://github.com/madebr/pyOpt (optional)

- To apply Tuckey window with alpha= 0.05, modify "apply_window = 0" to "apply_window = 1" in interface.py
	
- To run Correct@TDS, write in the command line:

	```
	python interface.py
	```
- The gui can freeze while computing, it's normal, it's still in development

## Usage

A THz software that performs statistical analysis and error correction on .he files obtained with TeraHertz Time Domain Spectroscopy (THz-TDS).
