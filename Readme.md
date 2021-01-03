### Read me

## Instruction for the simulations  
Make sure the packages are installed on the machine using `C00.00.install.R`.  
Then run `C02.01.simulation-Basic-ijk-batch.R` with the vector `args` properly defined.  
For example, this specification is done in the SLURM job submission code, `simrun.sh`.   


# C00 and C01 - preliminaries
C00.00.install.R - The R packages needed to be installed for these simulations.  
C01.01.*.R  - These are preliminary code for parameter selection based on the ZOE and the IBD data, and thus cannot be implemented without the data.  
C01.02.simulation.setup.R  -  Setting up the basic parameters for the simulations.  

# C02 - simulation implementation  
C02.01.simulation-Basic-ijk-batch.R - a modular code for simulations. This runs a specific setting (the baseline distribution, the disease effects, and the batch effects.)  
C02.02.*.R - Once all the simulations for each setting is completed, the outputs are summarized as a figure.

# C03 - Application  
C03.01.*.R - The tests based on the ZOE data. This code cannot be implemented without the data.  

# D0* and E0* - Not relevant.

# F0* - Functions 
These are needed to run the above code. They are sourced and called within each of the above code.   
