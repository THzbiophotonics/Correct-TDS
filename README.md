# Correct-TDS

| Information       | Links                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|:------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Project**       | [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) ![GitHub](https://img.shields.io/github/license/dodge-research-group/thztools) ![PyPI - Status](https://img.shields.io/pypi/status/thztools) ![PyPI - Version](https://img.shields.io/pypi/v/thztools) [![Anaconda-Server Badge](https://anaconda.org/conda-forge/thztools/badges/version.svg)](https://anaconda.org/conda-forge/thztools) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/thztools) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) [![Common Changelog](https://common-changelog.org/badge.svg)](https://common-changelog.org) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](code_of_conduct.md) [![DOI](https://zenodo.org/badge/569133241.svg)](https://zenodo.org/doi/10.5281/zenodo.10100093) [![pyOpenSci Peer-Reviewed](https://pyopensci.org/badges/peer-reviewed.svg)](https://github.com/pyOpenSci/software-review/issues/209) [![DOI](https://joss.theoj.org/papers/10.21105/joss.07542/status.svg)](https://doi.org/10.21105/joss.07542) |
| **Build**         | ![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/dodge-research-group/thztools/sphinx.yml?label=build%3Adocs) ![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/dodge-research-group/thztools/pytest-with-coverage.yml?label=build%3Atests%20(conda)) ![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/dodge-research-group/thztools/test-pip.yml?label=build%3Atests%20(pip)) [![codecov](https://codecov.io/gh/dodge-research-group/thztools/branch/dev/graph/badge.svg?token=U8PLKTQ7AH)](https://codecov.io/gh/dodge-research-group/thztools)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |

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
