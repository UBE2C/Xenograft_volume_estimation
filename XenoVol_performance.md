# XenoVol <br>
<br>

## Mean tumor volumes <br>

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










