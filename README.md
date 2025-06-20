# Code for the thesis "Laki to Tambora: Signal Restoration and Pattern Recognition in Ultra High Resolution Volcanic and Isotopic Signals" (2021) by T. Quistgaard
**DEAD LINKS**: This repository is no longer maintained, and the links to the data may be dead. The code is provided for reference and educational purposes only.

Link to the thesis:
- https://nbi.ku.dk/english/theses/masters-theses/thea-quistgaard/Thea_Quistgaard_MasterThesis.pdf

Project is inspired by the great volcanic eruptions of the passed, detectable in ice core signals - but also in literature, art and economics:
<p align="center" width="100%">
  <img width="30%" alt="Turner and Byron" src="https://github.com/user-attachments/assets/90f07289-0c8e-4541-8374-4e312cd66777" />
</p>

Example of how the isotopic signal between volcanic eruptions is located by matching up electric conductivity measurements and isotopes.

Focus is on six ice cores, drilled in southern Greenland (left). Measurements are focused on electric conductivity measurements and stable water isotopes (right). The timing of the two volcanic eruptions of Laki (1783-1784) and Tambora (1815) are well determined in literature and can be used as an optimizing goal for this project. The volcanoes are detected in f.ex. electric conductivity measurements, where spikes can be seen clearly when a volcanic eruption allows for more specific chemicals in the atmosphere. This is clearly seen on (right, below). 
<p align="center" width="100%">
  <img width="20%" alt="Greenland drill sites" src="https://github.com/user-attachments/assets/b70ba5f9-4804-4176-a7ac-1206ebba7f9e"/>
  <img width="50%" alt="Electric and isotopic example" src="https://github.com/user-attachments/assets/70001bae-5246-40ec-a3b0-e4c3de84a7f8" />
</p>



## Main modules implemented (in `MainModules/`):
General idea of the project: Restoring isotopic diffused signals by tuning the diffusion length until the expected number of years in that section of the timeseries is detectable in the depth interval.
<p align="center" width="100%">
  <img width="60%" alt="General algorithm" src="https://github.com/user-attachments/assets/2b7f00a8-fcc3-429f-94d8-04ec50324cfc" >
</p>


- `HL_AnalyticThea_class.py`: Herron-Langway densification model for ice cores
- `DiffusionProfiles_calculations.py` and `Diffusivity.py`: Diffusion length profile model
- `Interpolation_Class.py`: Interpolation methods for isotopic signals
- `SignalAttenuation.py`: Signal attenuation and layer thickness estimation
- `Decon.py`: Spectral transforms and analysis, aling with general deconvolution/back diffusion methods
- `BackDiffuse_LT.py`: The final optimization module, along with the constrained peak detection method
- `sigmaSolver.py` and `TemperatureEstimates.py`: Final temperature estimates and sigma solver for the back diffusion process

All modules are described in depth in the thesis, and the code is available for open science purposes.

All other modules are auxiliary and used for testing, plotting, or data processing. The main modules are the ones listed above.

## Timeseries and data
Example of all cores analysed, full time series. Blue section corresponds to signal located between the two eruptions of Laki and Tambora.

<p align="center" width="100%">
  <img width="50%" alt="Timeseries, all cores" src="https://github.com/user-attachments/assets/663b158f-69f7-4305-954b-516b584ead50" /> 
</p>


Example of all cores analysed, timeseries within the Laki and Tambora eruptions.

<p align="center" width="100%">
  <img width="50%" alt="Timeseries, all cores, zoom" src="https://github.com/user-attachments/assets/237b3f15-31f9-42ba-8f5c-81167b2cc742" />
</p>






## Algorithms
<p align="center" width="100%">
    <img width="33%" alt="Algorithm 1" src="https://github.com/user-attachments/assets/1711fb13-8706-4ed8-8f6f-ca0e7fc1cd63">
    <img width="33%" alt="Algorithm 2a" src="https://github.com/user-attachments/assets/b3341200-d31b-43c4-a3da-3b8b0d6d3402">
    <img width="33%" alt="Algorithm 2b" src="https://github.com/user-attachments/assets/69aa6004-0c60-41ac-bee9-0f6c2c89382e">
</p>

Example of fitting to power spectrum, with signal, noise and combined signal and noise fits.

<p align="center" width="100%">
  <img width="65%" alt="Fitting to PSD" src="https://github.com/user-attachments/assets/3004d927-0acb-487b-b5ca-df2134c93a10" />
</p>


Example of final back-diffused series, with corresponding 33 years between eruptions.
<p align="center" width="100%">
  <img width="50%" alt="Final back-diffused series" src="https://github.com/user-attachments/assets/178a6b38-d147-4853-a912-e09c2bdef8a7" />
</p>


Tables with final diffusion length estimates, based on different algorithmic choices.
<p align="center" width="100%">
    <img width="49%" alt="Table with results 1" src="https://github.com/user-attachments/assets/11f30b04-9736-4178-8098-ba51a92efcf4">
    <img width="49%" alt="Table with results 2" src="https://github.com/user-attachments/assets/a2464495-569c-4023-b152-9ce801e087e5">
</p>


Number of estimated peaks as a function of diffusion length, for all cores.
<p align="center" width="100%">
  <img width="50%" alt="Final back-diffused series" src="https://github.com/user-attachments/assets/de150af3-821e-4b0e-ac36-54b3610b637f" />
</p>



