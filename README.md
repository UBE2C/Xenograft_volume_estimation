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

The application requires the user to have R version 4.3.2 or later installed on their computer.<br>
*NOTE: For the ease of use the application should automatically dwnload and install any required packages when first launched, so don't be alarmed :).*


**How to use**

Clone this repository or download the SynthXenoGen R script, and lunch the application from the command line by using 'Rscript' followed by the script name and any required/desired option flags.

*NOTE: for Windows useres, you will need to add the R executable file path to your global 'Path' environment in order call 'Rscript' from the command line. You will find an amazing guide on 
how to do this on the documentation of the this.path R package ([here](https://www.rdocumentation.org/packages/this.path/versions/0.4.4/topics/Running.R.from.the.command-line)).
After you set up your 'Path' environment I would recommend to use PowerShell for app calls.*


**Call examples**


To see the available option flags and descriptions type:
```
Rscript SynthXenoGen_v1.0.0.r --help
```

To generate a dataset with 6 measurements (with generated dates) from 10 samples , with a growth variation of 30-60% using an exponential growth model type:
```
Rscript SynthXenoGen_v1.0.0.r -m 6 -s 10 --growth_model exponential  --request_dates TRUE --growth_variation_min 0.3 --growth_variation_max 0.6
```


### XenoVol


**Description**

XenoVol (Xenograft Volume Estimator) is designed to streamline the estimation of tumor xenograft volumes.
It implements the f-constant formula established by Feldman JP et al. (2010) and expanded upon by Sápi et al. (2015).
The primary purpose of XenoVol is to automate volume calculations, particularly for tumors where obtaining accurate height measurements is challenging or impractical,
such as for xenografts growing on the flank region. Instead, it relies on reference µCT measurements obtained throughout the experiments.


**Requirements**

The application requires the user to have R version 4.3.2 or later installed on their computer.<br>
*NOTE: Just like with SynthXenoGen the application should automatically dwnload and install any required packages when first launched.*

Additionally the input data must have the .csv extension and must be formatted in the following way:

Caliper measurements
- the name of the .csv file must contain the word caliper (not case sensitive) such as 'Control_caliper_measurements.csv'
- the measurement columns must contain the 'LxW' designation to indicate that they contain the product of the Length and Width measurements

µCT measurements
- the name of the .csv file must contain the word uCT (not case sensitive) such as 'Control_uCT_measurements.csv'
- the measurement columns must contain the 'uCT' designation to indicate that they contain the reference volumes.

For both .csv files
- The first 3 columns of the .csv files are fixed ID columns with the following designations:
  - Column 1: 'Treatment_group_ID'
    This column should contain the degination of the experimental group e.g.: Ctrl for the control group.
  - Column 2: 'Treatment'
    This column contain the treatment the gorup recieved, e.g.: control/none for the control group
  - Column 3: 'Mouse_ID'
    This This column contain the designation of the mice, found in each group, e.g.: S1 (for sample 1)
- There should be no spaces in the column names, preferably use underscores (_) instead.
- In the case of dates in the column names avoid using dot (2024.01.01) or forward slash (2024/01/01) notations and use hyphens (2024-01-01) or underscores (2024_01_01) instead.

An example for a µCT input file snippet:
|Treatment_group_ID|Treatment|Mouse_ID|uCT_2024-01-01|...|
|---|---|---|---|---|
|Ctrl|Control|Mouse_1|2.1234|...|

An example LxW input file snippet:
|Treatment_group_ID|Treatment|Mouse_ID|LxW_2024-01-01|...|
|---|---|---|---|---|
|Ctrl|Control|Mouse_1|1.4321|...|

Example .csv files can be found in the [Example_data](Example_data/) folder. Please use them as a reference for input file formatting.


**How to use**

Like with SynthXenoGen, clone this repository or download the XenoVol R script, and lunch the application from the command line by using 'Rscript' followed by the script name and any required/desired option flags.

*NOTE: again, for Windows useres, you will need to add the R executable file path to your global 'Path' environment in order call 'Rscript' from the command line.*


**Call examples**


To see the available option flags and descriptions type:
```
Rscript XenoVol_v1.0.1.r --help
```

To process the provided example data:
```
Rscript XenoVol_v1.0.1.r -i Input_files -o Output_files -v TRUE --quiet FALSE -r TRUE --outlier_handling detect -t mZscore_test -p rmse -c TRUE -m linear_interpolation --plot_theme dark
```


