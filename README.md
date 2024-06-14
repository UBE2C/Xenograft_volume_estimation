# Xenograft volume estimation

## About the repository

This repository contains the scripts for **SynthXenoGen** and **XenoVol**, two CLI applications/tools written in R.

The general purpose of these applications is to help researchers to estimate the volumes of various xenografts using the obtained caliper measurements
combined with reference measurements from µCT or MRI imaging.


## The Applications

### SynthXenoGen

**Description**

SynthXenoGen (Synthetic Xenograft Data Generator) is designed to generate synthetic sample data for the XenoVol script.
It produces Length and Width measurements (mimicking caliper measurements), LxW measurement products, and full volume measurements (mimicking µCT measurements)
and the associated growth curve plots. Note that the total volume calculation can be set to assumes an ellipsoid or hemi-ellipsoid shape.
Although it's main purpose is to create sample data, I hope  that it will also be useful for other forms of simple xenograft growth modeling under various control and treatment conditions.

**Requirements**

The application requires the user to have R version 4.3.2 or later installed on their computer.
It should automatically dwnload and install any required packages when first launched.

**How to use**

Clone this repository or download the SynthXenoGen R script, and lunch the application from the command line by using 'Rscript' followed by the script name and any required/desired option flags.

**Call examples**

'''
To see the available option flags and descriptions type:
Rscript SynthXenoGen_v1.0.0.r --help

To generate a dataset with 6 measurements (with generated dates) from 10 samples , with a growth variation of 30-60% using an exponential growth model type:
Rscript SynthXenoGen_v1.0.0.r -m 6 -s 10 --growth_model exponential  --request_dates TRUE --growth_variation_min 0.3 --growth_variation_max 0.6
'''

### XenoVol

**Description**

XenoVol (Xenograft Volume Estimator) is designed to streamline the estimation of tumor xenograft volumes.
It implements the f-constant formula established by Feldman JP et al. (2010) and expanded upon by Sápi et al. (2015).
The primary purpose of XenoVol is to automate volume calculations, particularly for tumors where obtaining accurate height measurements is challenging or impractical,
such as for xenografts growing on the flank region. Instead, it relies on reference µCT measurements obtained throughout the experiments.

**Requirements**

The application requires the user to have R version 4.3.2 or later installed on their computer.
It should automatically dwnload and install any required packages when first launched.

Additionally the input data must have the .csv extension and must be formatted in the following way:

Caliper measurements
- the name of the .csv file must contain the word caliper (not case sensitive) such as 'Control_caliper_measurements.csv'
- the measurement columns must contain the LxW designation to indicate that they contain the product of the Length and Width measurements

µCT measurements
- the name of the .csv file must contain the word uCT (not case sensitive) such as 'Control_uCT_measurements.csv'
- the measurement columns must contain the uCT designation to indicate that they contain the reference volumes.
- there must be at least two µCT measurements available for the application to run.

For both .csv files
- the first 3 columns of the .csv files are fixed ID columns with the following designations:
  Column 1: 'Treatment_group_ID'
    This column should contain the degination of the experimental group e.g.: Ctrl for the control group.
  Column 2: 'Treatment'
    This column contain the treatment the gorup recieved, e.g.: control/none for the control group
  Column 3: 'Mouse_ID'
    This This column contain the designation of the mice, found in each group, e.g.: S1 (for sample 1)






