Our method can be run on any machine capable of running MATLAB2021a with the Parallel Computing Toolbox. Instalation is trivial. Simply place the 'gene_exp_nonpara' folder in any directory accesible to your MATLAB installation.



This repository contains a folder, './Datasets', which contains all the data files that were used in the work:


Each data file is a .mat file containing the information required by the data. See any of the examples listed above for details of information required by the method. We note the fact that the variable named 'ground' is not required and is only used in the case of synthetic data, where ground truth parameters are known. In general, the method requires hand-input data with the following information:

data.M    : The maximum number of RNA obeserved across ALL time points.

data.t_f  : The duration of the experiment.

data.obs_t: 1D array of collection times of RNA counts (note that these two quantities will determine the units of all rates output by the method. i.e. if times are provided in seconds, rates will be in seconds^-1).

data.J    : The number of cells per time point. If this varies across time points, provide as a 1D array of the same shape as data.obs_t.

data.K    : The number of observation times.

units     : A structure containing relevant units. This is not used for computation, only for visualization purposes.


This repository contains 4 main MATLAB scripts, contined in './code':

run.m              : will run the method for 10000 iterations on an example synthetic dataset. (Expected demo wall time: ~5days)

run_continue.m     : will continue the method for an additional 10000 iterations on the same example synthetic dataset.

run_experimental.m : will run the method for 10000 iterations on an example experimental dataset.

visualize_results.m: will visualize results from the './results' folder.

Each of these scripts can be run on user-provided data, as long as it is formatted to match the examples provided.

During the execution of each scripts, a set of files is saved in the local directory './Results' with names determined by runID periodically. Due to the Parallel Tempering (PT) scheme involved in the sampler, results are saved as a folder containing a 'PTChainI.mat' file for PT chain I. For typical use, only the first chain 'PTChain1.mat' is required.

In addition to the Parallel Computing Toolbox, this code relies on two user-uploaded MATLAB codes retrieved from the MATLAB File Exchange:


Victor M. Garcia Molla (2022). Sensitivity Analysis for ODEs and DAEs (https://www.mathworks.com/matlabcentral/fileexchange/1480-sensitivity-analysis-for-odes-and-daes), MATLAB Central File Exchange. Last Retrieved June 28, 2022.

      This code is used during the Hamiltonian Monte Carlo proposal scheme of the Gibbs sampler.
      
Yi Cao (2022). Munkres Assignment Algorithm (https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm), MATLAB Central File Exchange. Retrieved June 28, 2022.
    
    
      This code is used to visualize the results of our PT scheme.
