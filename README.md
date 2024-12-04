# Correct-TDS

| Information       | Links                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|:------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Project**       | [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) ![GitHub](https://img.shields.io/github/license/dodge-research-group/Correct-TDS) ![PyPI - Status](https://img.shields.io/pypi/status/Correct-TDS) ![PyPI - Version](https://img.shields.io/pypi/v/Correct-TDS) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/Correct-TDS) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) [![Common Changelog](https://common-changelog.org/badge.svg)](https://common-changelog.org) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](code_of_conduct.md) [![DOI](https://zenodo.org/badge/569133241.svg)](https://zenodo.org/doi/10.5281/zenodo.10100093) [![pyOpenSci Peer-Reviewed](https://pyopensci.org/badges/peer-reviewed.svg)](https://github.com/pyOpenSci/software-review/issues/209) [![DOI](https://joss.theoj.org/papers/10.21105/joss.07542/status.svg)](https://doi.org/10.21105/joss.07542) |


## With the pip install package :

<ol>
  <li>Install pyOpt from https://github.com/madebr/pyOpt (optional)</li>
	
  <li>Install anaconda (https://docs.anaconda.com/anaconda/install/index.html)</li>
  
  <li>Create an environment with Anaconda</li>
	
  <li>Open the environment terminal</li>
	
  <li>Check the python version with the command line :</li>
	

	python --version


It should be between 3.9 and 3.12.

  <li>To install the package, write in the command line :</li>


	pip install Correct-TDS


  <li>To run Correct@TDS, write in the command line :</li>


 	Correct-TDS-interface

The gui can freeze while computing, it's normal, it's still in development

</ol>


## Manual instalation

<ol>
  <li>Install Anaconda: <a href="https://docs.anaconda.com/anaconda/install/index.html" target="_blank">https://docs.anaconda.com/anaconda/install/index.html</a></li>
  
  <li>Install mpi4py with the following command:</li>
  <pre><code>pip install mpi4py</code></pre>
  
  <li>Install pyswarm with the following command:</li>
  <pre><code>pip install pyswarm</code></pre>
  
  <li>Install numba with the following command:</li>
  <pre><code>pip install numba</code></pre>
  
  <li>Install pyOpt from <a href="https://github.com/madebr/pyOpt" target="_blank">https://github.com/madebr/pyOpt</a> (optional)</li>
  
  <li>To apply the Tuckey window with alpha = 0.05, modify <code>apply_window = 0</code> to <code>apply_window = 1</code> in <code>interface.py</code></li>
  
  <li>To run Correct@TDS, use the following command:</li>
  <pre><code>python interface.py</code></pre>
  
  The GUI can freeze while computing. This is normal; the application is still in development.
</ol>

## Usage

A THz software that performs statistical analysis and error correction on .he files obtained with TeraHertz Time Domain Spectroscopy (THz-TDS).

## Related software
Below is a list of other software projects that address related tasks in
THz-TDS analysis, with summaries taken from the project documentation.
- https://github.com/THzbiophotonics/Fit-TDS
  - "Python code aiming at the retrieving of material parameter from a
    TeraHertz time domain spectroscopy (TDS) measurements from a fit in the time
    domain."
- https://github.com/puls-lab/phoeniks
  - "A free and open-source (FOSS) Python class to extract the refractive index
    and absorption coefficient from time traces of a THz Time Domain
    Spectrometer (THz-TDS)."
- https://github.com/dotTHzTAG/CaTSper
  - "The CaTSper tool extracts the frequency-dependent optical constants from
    terahertz time-domain waveforms." (uses MATLAB)
- https://github.com/YaleTHz/nelly
  - "Nelly is a package for extracing refractice indices and conductivities from
    time-domain terahertz spectroscopy data." (uses MATLAB)
- https://github.com/dodge-research-group/thztools
  - "THzTools is an open-source Python package for data analysis in
    terahertz time-domain spectroscopy (THz-TDS)."
