The paper which describes the algorithms implemented here is available at: http://biorxiv.org/content/early/2016/06/30/061481

This code was written by Marius Pachitariu. It is provided here with no warranty. Please direct all questions and requests to marius10patgmaildotcom. 

Instructions for installing and running the software are provided in the Docs folder. 

Kilosort provides a results file called "rez", where the first column of rez.st3 are the spike times and the second column are the cluster identities. 

However, the best way to use this software is together with Phy (https://github.com/kwikteam/phy), which provides a manual clustering interface for refining the results of the algorithm. 
*** NOTE that you need to use a special branch of Phy with Kilosort. Instructions in Docs/phy_installation_with_templates.txt ***

You also need to install npy-matlab (https://github.com/kwikteam/npy-matlab), to provide read/write functions from Matlab to Python, because Phy is written in Python.

You can verify that the code has been installed correctly by running master_eMouse inside the eMouse folder. See first readme_eMouse.txt. You can also use these scripts to understand how to pass the right settings into Kilosort (will depend on your probe, channel map configuration etc), and what you should be seeing in Phy during manual cleanup of Kilosort results. There are many parameters of the simulation which you can tweak to make it harder or easier, and perhaps more similar to your own data. 