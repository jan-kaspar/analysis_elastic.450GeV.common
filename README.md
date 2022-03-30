# Basic usage

The following instructions have been tested in the LXPLUS7 environment.

Preparation: edit the `environment` file
  * source an appropriate CMSSW environment (tested with CMSSW_10_6_17),
  * set the `Elegent_DIR` variable to point to an installation of the [Elegent](https://elegent.hepforge.org/) package (used for the CNI calculations).

Initialisation:
```
bash --rcfile environment  # intialise environment
./configure  # prepare makefile and clang compilation dictionary
mk  # run compilation
```

Run fits to data:
```
./run_fits_data
```

Run MC study:
```
./run_simu_chain
```


# Structure of the repository

The repo is composed of four parts:
  * common code in the `classes` directory,
  * tool to import the differential cross-section results (script `run_input_import`),
  * tool to fit the experimental differential cross-section (script `run_fits_data`),
  * tool to study the fit performance with MC (script `run_simu_chain`).



# Common code

The shared code is in the `classes` directory. It includes the following modules:
  * `HadronicFitModel`: an implementation of hadronic model used in CNI fits, it supports
    * exponential modulus with up to 5 slope parameters
    * various phase t dependences (including constant)
  * `Result`: a structure holding the results of a CNI fit
  * `Stat`: a basic statistics toolbox
  * `command_line_tools`: routines to parse the command line arguments



# Differential-cross section importing

The input to the fits includes the differential cross-section from the beta* = 100m dataset and it may include also the one from the 11m dataset. For each dataset, the imported data consist of the central values, statistical uncertainties and systematics.

Note that the latest validated cross-section data are already included in this repo (committed to git), thus the rest of this section is only relevant if you wish to import new data.

As preparation, edit the `run_input_import` file to set up the input directories. Then you may execute
```
./run_input_import -version <tag>
```
This will import the data to the directory `data/input/<tag>`.



# Fitting

Fits of the experimental data are executed by
```
./run_fits_data -input <input version> -output <output version> -types <fit types>
```
It will take the differential cross-section input (cf. previous section) from the directory `data/input/<input version>` and save the output to the directory `data/fits/<output version>`. The last argument specifies which fit types (see later) shall be applied. For each fit type a separate sub-directory is made.

The fit types are defined in the `run_standard_fit` script. There are two major fit approaches:
  * _single_: all fit parameters are optimised in a single fit
  * _sequential_: a series of fits is made, each fit optimising a sub-set of parameters over a dedicated t region

The other configurable parameters include t fit range, degree of the B polynomial etc.

The fitting code can be found in these programs:
  * `make_fit`: systematics are treated via a covariance matrix
  * `make_fit_nuis`: systematics are treated via nuissance parameters

The structure of both programs is similar. Let us comment on the `make_fit` program.
  * `Model`: encapsulation of CNI modelling. Beyond others, it has a cache for the CNI corrections (`psi` in the code). This correction weakly depends on the fit function, thus it makes sense to update it only in between fit iterations. This greatly improves the fit speed.
  * `Dataset`: all data related to a dataset. There can be up to two datasets (beta* = 100 and 11m). Some of the data members may need updating during the fit iterations (see the `Update` method):
    * bin representative points
    * conversion from relative to absolute systematic uncertainties
  * `Data`: encapsulation of all data which are compared to the `Model`
  * `Metric`: class to evaluate the match between `Model` and `Data`, see especially the `CalculateChi2` method.
  * `Minimization`: a wrapper around ROOT routines to minimise the `Metric`.
  * `main`: runs multiple fit iterations and updates relevant parameters in between iterations

For each fit there are multiple outputs:
  * `result` file: summary of the most important parameters
  * `root` file: graphical, monitoring and debugging results
  * `log` file: logging information



# MC simulation

The MC study comprises the following steps.
  * For each dataset (beta* = 100 and 11m), random simulation of experimental differential cross-section obseravations = central values + fluctuations reflecting statistical, systematic and normalisation uncertainties (with due correlations where applicable). This is done by the `simu_make_histograms` program.
  * For each simulated cross-section, apply the selected fits. This uses the same script as the data fits, `run_standard_fit`.
  * Repeat the simulation multiple times and collect the results.
     * Analyze the fluctuations in the simulated histograms, essentially to validate if they follow the input distributions. This is done by the `simu_process_fits`.
     * Analyze the fit results, see the `simu_process_fits`. This allows to assess on the performance of various fit types, i.e. bias, RMS, ...

By running the simulation chain with different perturbations (i.e. error models), one can assess the reaction of the fits to various sources of errors: statistical uncertainties, systematics, ...

The full simulation chain is executed with
```
./run_simu_chain
```
Add `--help` flag to get the usage hint: you can specify the number of repetitions, the error model, the input dsigma/dt model, etc.

# Plots

For many standard plots, the generating scripts are placed in the `plots` directory. They are split by topic into a hierarchy of subdirectories. The scripts written in the Asymptote language. The similarity to C++ should make them at least readable.



# Testing

The repo includes a test suite in the `test` directory. The tests are triggered by executing
```
./run_test
```
By default the test will check
  * that there are no left FIXME and TODOs in the code
  * that compilation is successfull
  * run static analysis
  * run a differential cross-section import
  * run a series of data fits
  * run the simulation chain (a short version)

Altogether, the test takes few minutes.

It can also compare plots (main analysis and systematics) wrt. another version - for that pass `-ref <directory>` to the `run_test` script. The reference plots are carefully selected to cover most steps in the analysis.

The list of test steps can be configured with `-steps <list of steps>` flag.

Tip: I've set up `acron` to run the tests every night and send me an email with the results. I also triggered the tests manually at every significant change of the code to avoid any regression.