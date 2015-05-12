# dp_spectra
This repository contains codes, notebooks, and links to data for "[Mesoscale to submesoscale wavenumber spectra in Drake Passage](http://crocha700.github.io/pdfs/dp_spectra_submitted_twocols.pdf)" by C. B. Rocha, T. K. Chereskin, S. T. Gille, and D. Menemenlis (submitted to JPO). An index notebook to guide you through is available [here](http://nbviewer.ipython.org/github/crocha700/dp_spectra/blob/master/index.ipynb). 

Codes developed and used in this project are available in the **src** directory. The most important codes were collected into a legit Python package, [**pyspec**](https://github.com/crocha700/pyspec). The repository for this package contains documentation, tests, and IPython notebooks showcasing the use of the most common functions.

# Abstract
To determine tne the dominant governing dynamics at mesoscales to submesoscales, we calculate upper ocean
(0-200 m) horizontal wavenumber spectra in the Antarctic Circumpolar Current in Drake Passage from 13
years of shipboard ADCP measurements, altimeter data, and a new internal-tide-resolving and submesoscaleadmitting
Massachusetts Institute of Technology general circulation model simulation. At scales between 10
and 200 km, the ADCP kinetic energy spectra approximately follow a $k^{-3}$ power law. While the observed
flows are more energetic at the surface, the shape of the kinetic energy spectra is independent of depth. These
characteristics resemble predictions of isotropic interior quasigeostrophic turbulence. The ratio of across-track
to along-track kinetic energy spectra, however, significantly departs from the expectation of isotropic interior
quasigeostrophic turbulence. The inconsistency is dramatic at scales smaller than 40 km. A Helmholtz decomposition
of the one-dimensional ADCP spectra shows that ageostrophic flows account for the discrepancy
between the observed spectra and predictions of isotropic interior quasigeostrophic turbulence. This conclusion
is supported by analyses of synthetic and numerical model data. In Drake Passage, unbalanced motions
account for about half of the near-surface kinetic energy at scales between 10 and 40 km. Model results indicate
that ageostrophic flows imprint on the sea surface. At scales between 10 and 40 km, unbalanced motions
account for about half of the model sea surface height variance.

# Authors contributions

* [Cesar B Rocha](crocha700.github.io) <<crocha@ucsd.edu>>
* [Teresa K. Chereskin](http://tryfan.ucsd.edu)
* [Sarah T. Gille](http://www-pord.ucsd.edu/~sgille/)
* [Dimitris Menemenlis](https://science.jpl.nasa.gov/people/Menemenlis/)

Chereskin and Gille conceived the Drake Passage wavenumber spectra project, and performed preliminary analysis. Menemenlis performed the global numerical simulations together with other MITgcm developers. Rocha analyzed the ADCP, synthetic, and model data supervised by Chereskin and Gille. Gille analyzed the altimeter data. Rocha led the paper write up, with significant contributions from Chereskin and Gille. Menemenlis reviewed the text of Section 6 and the discussion relative to the model results.

# Funding
This study was funded by the NASA Ocean Surface Topography Science Team (NNX13AE44G) and the SWOT Science Definition Team (NNX13AE85G), and the NSF Polar Program (PLR-1341431). C. B. Rocha was partially supported by NSF (OCE 1357047).

