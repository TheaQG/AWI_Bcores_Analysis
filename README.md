# Code for the thesis "Laki to Tambora: Signal Restoration and Pattern Recognition in Ultra High Resolution Volcanic and Isotopic Signals" (2021) by T. Quistgaard
**DEAD LINKS**: This repository is no longer maintained, and the links to the data may be dead. The code is provided for reference and educational purposes only.

Link to the thesis:
- https://nbi.ku.dk/english/theses/masters-theses/thea-quistgaard/Thea_Quistgaard_MasterThesis.pdf

## Main modules implemented (in `MainModules/`):
- `HL_AnalyticThea_class.py`: Herron-Langway densification model for ice cores
- `DiffusionProfiles_calculations.py` and `Diffusivity.py`: Diffusion length profile model
- `Interpolation_Class.py`: Interpolation methods for isotopic signals
- `SignalAttenuation.py`: Signal attenuation and layer thickness estimation
- `Decon.py`: Spectral transforms and analysis, aling with general deconvolution/back diffusion methods
- `BackDiffuse_LT.py`: The final optimization module, along with the constrained peak detection method
- `sigmaSolver.py` and `TemperatureEstimates.py`: Final temperature estimates and sigma solver for the back diffusion process

All modules are described in depth in the thesis, and the code is available for open science purposes.


All other modules are auxiliary and used for testing, plotting, or data processing. The main modules are the ones listed above.