# Simulation Software of the detection process of the directional wave spectrum of SWIM
This is a software for simulatiing the detection process SWIM. The process steps follows the instructions in Hauser et al. 2001. For more details about data flow and processing methods please refer to this paper.
* Hauser D, Soussi E, Thouvenot E, et al. SWIMSAT: A real-aperture radar to measure directional spectra of ocean waves from space--main characteristics and performance simulation. Journal of Atmospheric & Oceanic Technology, 2001, 18(3):421-437.


## Configuration File `ParaConfig.json`:
Revising the radar parameters in `ParaConfig.json` such as flying height of radar, 3dB width of radar, and the radar gates to change the size of simulated ocean area. Revising the Experiment parameters to change the exp rounds or the scanning interval of SWIM.

## Running
After preparation of the configuration file, you can run the software simply by entering `main` in the command line of MATLAB. 3 figures would be shown at the end of running. As follows:
