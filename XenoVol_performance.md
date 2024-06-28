# Introduction

This document will focus on the predictive performance of XenoVol by comparing the results of the original f-constant <br/>
formula and the self developed correction methods, and by comparing all these to the results of other, commnly used volume <br/>
prediction formulas. <br>

This is going to be a relatively long read, so if you have no time or desire to go through the whole document, the TL;DR is: <br/>
Based on the aforementioned tests XenoVol offers a more precise, easy to use way to predict xenograft volumes compared to the <br/>
other tested formulas, and the correction methods will increase the prediction accuracy even further compared to the original <br/>
f-constant formula :). <br>

<br>

# Correct or not to correct... <br>
<br>

## Real mean tumor volumes <br>

First in order to see how well the original model and the applied 'homebrew' correction methods predict xenograft volumes, <br/>
I calculated the mean tumor volumes of all collection time points, for all three generated datasets (Control dataset - <br/>
unimpeded growth, Continous-treatment dataset - continous, effective treatment and Single-treatment dataset - a single <br/>
treatment event). Figure 1 A-C depicts the visualized real mean tumor volumes. <br>

### Figure 1A <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_ctrl_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_ctrl_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_ctrl_p_dark.png">
</picture>
<br/>

*The figure shows the mean real xenograft volumes and their change over time in the Control dataset. Error bars are showing <br/>
the standard error of the means* <br>

### Figure 1B <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G1_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G1_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G1_p_dark.png">
</picture>
<br/>

*The figure shows the mean real xenograft volumes and their change over time in the Continous-treatment dataset. Error bars <br/>
are showing the standard error of the means* <br>

### Figure 1C <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G2_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G2_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G2_p_dark.png">
</picture>
<br/>

*The figure shows the mean real xenograft volumes and their change over time in the Single-treatment dataset. Error bars <br/>
are showing the standard error of the means* <br>

<br>

## Real vs predicted mean tumor volumes <br>

Then I overlayed the predicted tumor volumes I obtained by running XenoVol with two reference dataset, with no correction <br/>
(volumes predicted by the original formula) and either with the 'mean_correction' or the 'linear_interpolation' methods. <br/>
By doing this I was  hoping to see how closely the predicted volumes (either with or woithout correction follow the real real <br/>
xenograft volumes over time, and to see if the correction methods improve the prediction accuracy. Figure 2 A-C depicts the <br/>
visualized real and predicted mean tumor volumes. <br>

### Figure 2A <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_ctrl_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_ctrl_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_ctrl_p_dark.png">
</picture>
<br/>

*The figure shows the mean real and predicted xenograft volumes overlayed, and their change over time in the Control <br/>
dataset. Error bars are showing the standard error of the means* <br>

### Figure 2B <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G1_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G1_p_light.png">
<img alt="Shows the light or dark mode version of the same plot."
  src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G1_p_dark.png">
</picture>
<br/>

*The figure shows the mean real and predicted xenograft volumes and their change over time in the Continous-treatment <br/> 
dataset. Error bars are showing the standard error of the means* <br>

### Figure 2C <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G2_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G2_p_light.png">
<img alt="Shows the light or dark mode version of the same plot."
  src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G2_p_dark.png">
</picture>
<br/>

*The figure shows the mean real and predicted xenograft volumes and their change over time in the Single-treatment <br/> 
dataset. Error bars are showing the standard error of the means* <br>

The resulting graphs showed that the predicted volumes track the real volumes quite closely and based on the error bars <br/>
there seems to be no significant difference between the predicted and real tumor volumes (although a proper significance test <br/>
was not done to confrm this). Nonetheless it is important to note that it seems that applying one of the correction methods <br/>
seems to increase the prediction precision. <br>

<br>

## The effect of correction on prediction precision <br>

In order to see if the previous observations hold true, and if the correction methods have any tangible effect on the <br/>
precision of the predictions, I calculatated the Mean Absolute Percentage Error (MAPE) for all three XenoVol runs, and <br/>
compared them to each other. Figure 3 A-C depicts the visualized MAPE values of the estimations with and without correction. <br>

### Figure 3A <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_ctrl_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_ctrl_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_ctrl_p_dark.png">
</picture>
<br/>

*The figure shows the MAPE values of XenoVol ran with and without correction on the Control dataset.* <br>

### Figure 3B <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_G1_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_G1_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_G1_p_dark.png">
</picture>
<br/>

*The figure shows the MAPE values of XenoVol ran with and without correction on the Continous-treatment dataset.* <br>

### Figure 3C <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_G2_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_G2_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_G2_p_dark.png">
</picture>
<br/>

*The figure shows the MAPE values of XenoVol ran with and without correction on the Single-treatment dataset.* <br>



