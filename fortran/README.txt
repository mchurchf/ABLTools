This set of tools was written by Sang Lee and given to me (Paula Doubrawa) in 2016. They convert *.xy files to a TurbSim *.bts file.

The *.xy files need to be in individual time directories. They are composed of N rows (where N is the number of points) and 6 columns which give the X, Y, Z, U, V, W values of each point.

I modified the original code slighlty to be able to handle fractional time directories, and it is hard-coded for directories with 3 decimal points in precision.

It is also hard-coded to read in a file named "tsConv.inp" which should be in the parent directory hosting the individual time sub-directories. The structure of the "tsConv.inp" file is as follows:

    0                   : theta
    107                 : numgrid_y (see TurbSim manual to understand this)
    107                 : numgrid_z (see TurbSim manual to understand this)
    1.06                : gridres_z (see TurbSim manual to understand this) 
    1.5                 : gridres_y (see TurbSim manual to understand this)
    1                   : timestep
    100                 : how many times to process
    9.0                 : reference wind speed (see TurbSim manual to understand this)
    80.0                : height of reference wind speed (see TurbSim manual to understand this)
    owez                : outfile prefix
    sampleArray1_U.xy   : name of xy files in each time directory
    15201               : first time (==name of first directory where to find *.xy file) 

To compile, simply type
	$ ./makeof2fast

To clean it all and undo the compilation, simply type
	$ ./cleanobj
