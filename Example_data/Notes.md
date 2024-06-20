# A note for the example data provided in this folder
<br>
<br>

## Experimental group designations
<br>

**Control** - Stands for a set of control measurements, representing an untreated group of xenografts. Therefore, unimpeded growth is expected here.
<br>

**Group 1 (G1)** - Stands for a set of treatment measurements, representing a group of xenografts receiving continuous treatment throughout the experiment.
Therefore, a continuous reduction in tumor mass is expected here.
<br>

**Group 2 (G2)** - Stands for a set of treatment measurements, representing a group of xenografts receiving a single dose of treatment at an early timepoint, with no further doses for the rest of the experiment.
Therefore, a slight reduction in tumor mass is expected here, coupled with recovery, albeit reaching a lower end-volume.
<br>
<br>

## Sample generation
<br>

To generate the aforementioned samples, SynthXenoGen was run using the following flags:
<br>

Control samples
```
Rscript SynthXenoGen_v1.0.0.r -m 6 -s 10 --growth_model exponential  --request_dates TRUE --growth_variation_min 0.3 --growth_variation_max 0.6
```
<br>

G1 samples
```
Rscript SynthXenoGen_v1.0.0.r -m 6 -s 10 --growth_model exponential --request_dates TRUE --treatment_requested TRUE --treatment_variation_min -0.3 --treatment_variation_max -0.1 --continuous_treatment TRUE --growth_variation_min 0.3 --growth_variation_max 0.6 --treatment_group_id G1 --treatment_type continuous
```
<br>

G2 samples
```
Rscript SynthXenoGen_v1.0.0.r -m 6 -s 10 --growth_model exponential --request_dates TRUE --treatment_requested TRUE --treatment_variation_min -0.6 --treatment_variation_max -0.1 --continuous_treatment FALSE --growth_variation_min 0.3 --growth_variation_max 0.6 --treatment_group_id G2 --treatment_type Single
```
