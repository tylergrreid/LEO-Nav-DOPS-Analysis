# LEO Navigation Geometry Analysis #

This code analyzes the satellite navigation geometry of the GNSS core constellations and Low Earth Orbit (LEO) Walker constellations. This was developed during my PhD thesis undertaken in the GPS Research Lab in the Department of Aeronautics and Astronautics at Stanford University: 

[1]	T. G. R. Reid, "Orbital Diversity for Global Navigation Satellite Systems," Doctor of Philosophy, Aeronautics and Astronautics, Stanford University, Stanford, CA, 2017.

This thesis is available at the following link: https://purl.stanford.edu/dc409wn9227

Other papers based on this code can be found here: 

[2]	T. G. R. Reid, A. M. Neish, T. Walter, and P. Enge, “Broadband LEO Constellations for Navigation,” Navigation, Journal of The Institute of Navigation, 2018. [link: https://onlinelibrary.wiley.com/doi/pdf/10.1002/navi.234]

[3]	T. Reid, A. Neish, T. Walter, and P. Enge, “Leveraging Commercial Broadband LEO Constellations for Navigation,” in Proceedings of the 29th International Technical Meeting of the Satellite Division of The Institute of Navigation (ION GNSS+ 2016), Portland, OR, 2016. [link: http://web.stanford.edu/group/scpnt/gpslab/pubs/papers/Reid_IONGNSS_2016_LEO.pdf]


## How to Use ##

'MAIN_GDOP_Comparison.m' is the main simulation code. This creates the grid of user locations on Earth and simulates satellite visibility and geometry (Dilution of Precision). This analysis is done for the GNSS core constellations of GPS, GLONASS, Galileo, and BeiDou as well as a specified group of Low Earth Orbit (LEO) Walker constellations. The output is a series of *.mat files in .../results/simulation_data/. 

'MAIN_GDOP_Comparison_Plot.m’ analyzes and plots the simulation data produced by 'MAIN_GDOP_Comparison.m’. The plots are saved under .../results/plots. 

## More Info ## 

For more info, please refer to the Stanford University GPS Research Lab in the Department of Aeronautics and Astronautics: https://gps.stanford.edu

This also uses some functionality of the Stanford GPS Lab’s Matlab Algorithm Availability Simulation Tool (MAAST) available for download here: https://gps.stanford.edu/resources/tools/maast