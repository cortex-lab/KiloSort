#Fast spike sorting for hundreds of channels#

Implements an integrated template matching framework for detecting and clustering spikes from multi-channel electrophysiological recordings. Very fast when a GPU is available, but can also run on the CPU side. Described in this paper, available as preprint

Pachitariu M, Steinmetz NA, Kadir S, Carandini M and Harris KD (2016). Kilosort: realtime spike-sorting for extracellular electrophysiology with hundreds of channels. 
bioRxiv dx.doi.org/10.1101/061481, [link](http://biorxiv.org/content/early/2016/06/30/061481). 

### Installation ###
If you are running on the GPU, you must run mexGPUall in the CUDA folder after setting up mexcuda in Matlab ([instructions](http://uk.mathworks.com/help/distcomp/mexcuda.html)). More detailed instructions for installing and running the software are provided in the Docs folder.

You can verify that the code has been installed correctly by running master_eMouse inside the eMouse folder. See first readme_eMouse.txt. You can also use these scripts to understand how to pass the right settings into Kilosort (will depend on your probe, channel map configuration etc), and what you should be seeing in Phy during manual cleanup of Kilosort results. There are many parameters of the simulation which you can tweak to make it harder or easier, and perhaps more similar to your own data. 

### Integration with Phy GUI ###
Kilosort provides a results file called "rez", where the first column of rez.st3 are the spike times and the second column are the cluster identities. However, the best way to use this software is together with [Phy](https://github.com/kwikteam/phy), which provides a manual clustering interface for refining the results of the algorithm. 

NOTE that you need to use a special branch of Phy with Kilosort. Instructions in Docs/phy_installation_with_templates.txt 

You also need to install [npy-matlab](https://github.com/kwikteam/npy-matlab), to provide read/write functions from Matlab to Python, because Phy is written in Python.

Detailed instructions for interpreting results are provided [here](https://github.com/kwikteam/phy-contrib/blob/master/docs/template-gui.md).

### Questions ###

Please create an issue for bugs / installation problems. 

### Licence ###

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

