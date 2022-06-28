# Correct-TDS-in-development

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
- To apply Tuckey window with alpha= 0.05, modify "apply_window = 0" to "apply_window = 1" in interface.py
	
- To run Correct@TDS, write in the command line:

	```
	python interface.py
	```
- The gui can freeze while computing, it's normal, it's still in development
