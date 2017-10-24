Gadget P(k) reader

This code reads the raw unbinned P(k) from Gadget and performs the binnning and the overlapping of the top mesh and the folded mesh.
It can also find the optimal binning and overlapping parameters, such that the ratio of power spectra has a smooth derivative.
There are also some python scripts for plotting.

To test if code is working, run in this order:
     python PKbinneroverlapSingle.py     
     python PKbinneroverlapRatioOptimizer.py
     python PKbinneroverlapSingle.py    (this time, parameters will be read from file produced above)
     python Plots/plot.PowerSpectra.py
     python Plots/plot.ManyPowerSpectra.py

    pkbinner.py:
    contains the needed functions for this. (Inspired in the original IDL script)

    PKbinneroverlapSingle.py:
    performs the binning and overlapping for the given set of parameters (at the top of the file) or it can also read parameters from the bestOverlapParams.txt

    PKbinneroverlapRatioOptimizer.py:
    performs the binning on two sets of simulation models (model L and model E) for each redshift (snapshot number), changing the parameters of the overlapping, such that the ratio of power spectra E/L has the smoothest derivative possible. The set of parameters which perform the best are saved on a file bestOverlapParams.txt

Plot routines:
    
    tools.py: some tools for ratios and derivatives
    
    plot.PowerSpectra.py:  performs the log-log plot on a single power spectrum specified through the simulation and snapshot index.
    plot.ManyPowerSpectra.py:  useful to plot several power spectra in a single figure, for example two sets of simulations and many redshifts. All spectra have to be in same folder.

Default Parameters of binning and parameters of overlapping can be changed in the PKbinner*.py files themselves. The utilities in pkbinner.py should be general enough.
Folder names and simulation snapshots can be also changed in the files themselves.


Acknowledgments:

This script was done for my Master thesis at the University of Heidelberg, under the supervision of Prof. Luca Amendola and Dr. Marco Baldi. It was used for the publication: Fitting and forecasting coupled dark energy in the non-linear regime, JCAP 1601 (2016) no.01, 045, e-Print: arXiv:1508.07208 [astro-ph.CO].

The data in this repository comes from the CoDECS database of cosmological simulations: http://www.marcobaldi.it/web/CoDECS_summary.html

The script is inspired on the original IDL scripts used in Gadget2. Thanks also to Prof. Volker Springel for useful comments.
 

