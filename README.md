# dp_spectra
This repository contains codes, notebooks, and links to data for "Mesoscale to submesoscale wavenumber spectra in Drake Passage" by C. B. Rocha, T. K. Chereskin, S. T. Gille, and D. Menemenlis (submitted to JPO). An index notebook to guide you through is available [here](http://nbviewer.ipython.org/github/crocha700/dp_spectra/blob/master/index.ipynb). 

Codes developed and used in this project are available in the **src** directory. The most important codes were collected into a legit Python package, [**pyspec**](https://github.com/crocha700/pyspec). The repository for this package contains documentation, tests, and IPython notebooks showcasing the use of the most common functions.

# Abstract
   To determine the dominant governing dynamics at meso to submeso scales, we calculate upper ocean (0-200 m) horizontal wavenumber spectra in the Antarctic Circumpolar Current in Drake Passage from 13 years of shipboard ADCP measurements, altimeter data, and a new internal-tide-resolving and submesoscale-admitting Massachusetts Institute of Technology general circulation model simulation. At scales between 10 and 200 km, the  ADCP kinetic energy spectra approximately follow a $k^{-3}$ power law. While the observed flows are more energetic at the surface, the shape of the kinetic energy spectra is independent of depth. These characteristics are reminiscent of predictions of isotropic interior quasigeostrophic turbulence. The ratio of across-track to along-track kinetic energy spectra, however, significantly departs from the expectation of isotropic interior quasigeostrophic turbulence. The inconsistency is dramatic at scales smaller than 40 km. A Helmholtz decomposition of the one-dimensional ADCP spectra shows that ageostrophic flows account for the discrepancy between observed spectra and isotropic interior quasigeostrophic turbulence predictions. This conclusion is supported by analyses of synthetic and numerical model data. In Drake Passage, unbalanced motions account for about half of the near-surface kinetic energy at scales between 10 and 40 km. Model results indicate that ageostrophic flows imprint on the sea surface. At scales beteween 10 and 40 km, unbalanced motions account for about half of the model sea-surface height variance.

# Authors contributions

* Cesar B Rocha <<crocha@ucsd.edu>>
* Teresa K. Chereskin
* Sarah T. Gille
* Dimitris Menemenlis

Chereskin and Gille conceived the project, and performed preliminary analysis. Menemenlis performed the global numerical simulations. Rocha analyzed the ADCP and model data supervised by Chereskin and Gille. Gille analyzed the altimeter data. Rocha led the paper write up, with contributions from Chereskin and Gille. Menemenlis contributed to the model setup description.

# Funding
This study was funded by the NASA Ocean Surface Topography Science Team (NNX13AE44G) and the NSF Polar Program (PLR-1341431). C. B. Rocha was partially supported by  NSF (OCE 1357047).
