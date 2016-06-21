The paper which describes the algorithms implemented here is available at:

This code was written by Marius Pachitariu. It is provided here with no warranty. Please direct all questions and requests to marius10patgmaildotcom. 

Instructions for installing and running the software are provided in the Docs folder. 

Kilosort provides a results file called "rez", where the first column of rez.st3 are the spike times and the second column are the cluster identities. 

However, the best way to use this software is together with Phy (https://github.com/kwikteam/phy), which provides a manual clustering interface for refining the results of the algorithm.

You also need to install npy-matlab (https://github.com/kwikteam/npy-matlab), to provide read/write functions from Matlab to Python, because Phy is written in Python.