# _ardor_
_ardor_ is a multi-tiered flare detection pipeline that aims to extract and characterize stellar flares from TESS 2-minute and 20-second cadence data.

All relevant functions used in the below functions are found in _Flare.py_, and an example test case can be found in `example.py`.

## Tier 0 - Detrending
The tier0() function only requires the .fits file of interest: time, PDCSAP_Flux, Detrended_Flux, PDCSAP_Error = tier0(TESS_fits_file). The output is four numpy arrays: time (BJD), PDCSAP flux (electrons/s), detrended PDCSAP flux (electron/s, median centered), and the PDCSAP Flux Error, respectively.

## Tier 1 - Flare Candidate ID
The tier1() function requires the numpy array of the detrended flux (assumed to be the output from Tier 0, but feel free to use your own detrending methods), as well as the sensitivity of the flares you are interested in, given in terms of standard deviation away from the baseline flux. The default is 3. 
The syntax is: tier1(detrended_flux, sigma=3) = flare_index, flare_length. The output gives two numpy arrays. The first is a list of flare candidate **_indices_**, stating at what index in the pdcsap flux array the flare candidate occurs. To get the BJD time of the candidate flare, simply use the output of Tier 0, and pass time[flares[0]], etc. The second output is a 1:1 list that approximates the flare duration by measuring how many data points after the initial peak the data takes to return to baseline.

## Tier 2 - Coarse Model Fitting
The tier2() function requires a few more inputs, and has the following syntax: tier2(time, flux, pdcsap_error, flares, lengths, output_dir, host_name = 'My_Host', T = 4000, host_radius = 1), and returns .csv files of flare candidates. The time, flux (NOT detrended_flux), and pdcsap_error are all from the tier0 output. The flares and lengths lists are from the tier1 output. The output_dir is the folder in which you want to put the output .csv files. The host_name is a string that represents what the target is called, and writes the csv file titles accordingly. The T value is the effective temperature of the host, in Kelvin. The host_radius is the radius of the target, in solar radii. The output csv will have three unlabeled columns. The first column is "Time since Flare Peak (min)", the second is "Median Relative Flux", and the third is "Median relative PDCSAP Error". These csv files are the input for tier3.

## Tier 3 - _allesfitter_ MCMC
The tier3() function requires a few directory inputs, but otherwise requires the same information from tier2. The syntax is 
tier3(tier_2_output_dir, tier_3_working_dir, tier_3_output_dir, settings_template_dir, params_template_dir, host_name = 'My_Host', T = 4000, host_radius = 1, MCMC_CPUS = 1). The first directory should be the same as that passed to the tier2 function. The tier_3_working directory should be a new folder, as files will be copied and deleted in and out as the MCMC runs from file to file. If an error arises, be sure to clear this directory out before rerunning. The tier_3_output directory is the folder where you want the final MCMC products to be put. The `settings_templare_dir` is the directory where the settings.csv file template is found. Example templates are found in the 'templates' folder in the repo. A similar story for the params.csv file. Some settings are easily modified from these templates (i.e., MCMC walkers, MCMC steps), which you can feel free to mess around with. If you want to create your own custom settings/params files and want to know the different settings available from _allesfitter_, check out the _allesfitter_ documentation at https://www.allesfitter.com/ (Günther & Daylan, 2019 and 2021; https://github.com/MNGuenther/allesfitter). The host name, T value, and host radius should be the same values passed in tier2. Finally, the MCMC_CPUS argument is the integer number of CPUs you want to run the MCMC on. ((NOTE: Multicore processing is NOT supported on Windows; this can be bypassed by using WSL))



# Acknowledgements

The development of DarkInferno has been supported by the McDonnell Center for Space Sciences at Washington University in St. Louis.
