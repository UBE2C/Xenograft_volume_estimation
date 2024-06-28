# Introduction <br>

This document will focus on the predictive performance of XenoVol by comparing the results of the original f-constant <br/>
formula and the self-developed correction methods, and by comparing all these to the results of other, commonly used volume <br/>
prediction formulas. <br>

This is going to be a relatively long read. If you don't have time to go through the whole document, here is the TL;DR: <br/>
Based on the aforementioned tests XenoVol offers a more precise, easy to use way to predict xenograft volumes compared to the <br/>
other tested formulas. Additionally, the correction methods further increase prediction accuracy compared to the original <br/>
f-constant formula. :). <br>

**NOTE: All the tests were performed using synthetically generated data; therefore, all the results and conclusions in this <br/>
document are also based on synthetically generated data. It is possible that real-life experimental data will generate very <br/>
different results with much worse estimations and precision, so please keep that in mind before using XenoVol.** <br>

<br>

# XenoVol - Correct or not to correct... <br>
<br>

## Real mean xenograft volumes <br>

First, in order to see how well the original model and the homebrew correction methods predict xenograft volumes, <br/>
I calculated the mean tumor volumes at all collection time points, for all three generated datasets (Control dataset - <br/>
unimpeded growth, Continuous-treatment dataset - continuous, effective treatment, and Single-treatment dataset - a single <br/>
treatment event). Figure 1A-C depicts the visualized real mean tumor volumes. <br>

### Figure 1A <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_ctrl_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_ctrl_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_ctrl_p_dark.png">
</picture>
<br/>

*The figure shows the mean real xenograft volumes and their change over time in the Control dataset. Error bars represent <br/>
the standard error of the means* <br>

### Figure 1B <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G1_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G1_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G1_p_dark.png">
</picture>
<br/>

*The figure shows the mean real xenograft volumes and their change over time in the Continuous-treatment dataset. Error bars <br/>
are represent the standard error of the means* <br>

### Figure 1C <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G2_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G2_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/real_vol_G2_p_dark.png">
</picture>
<br/>

*The figure shows the mean real xenograft volumes and their change over time in the Single-treatment dataset. Error bars <br/>
are represent the standard error of the means* <br>

<br>

## Real vs predicted mean xenograft volumes <br>

Then, I overlaid the predicted tumor volumes I obtained by running XenoVol with two reference datasets, with no correction <br/>
(volumes predicted by the original formula) and either with the 'mean_correction' or the 'linear_interpolation' methods. <br/>
By doing this I was  hoping to see how closely the predicted volumes (either with or without correction) follow the real <br/>
xenograft volumes over time, and to see if the correction methods improve the prediction accuracy. Figure 2A-C shows the <br/>
visualized real and predicted mean tumor volumes. <br>

### Figure 2A <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_ctrl_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_ctrl_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_ctrl_p_dark.png">
</picture>
<br/>

*The figure shows the mean real and predicted xenograft volumes overlaid, and their change over time in the Control <br/>
dataset. Error bars represent the standard error of the means* <br>

### Figure 2B <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G1_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G1_p_light.png">
<img alt="Shows the light or dark mode version of the same plot."
  src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G1_p_dark.png">
</picture>
<br/>

*The figure shows the mean real and predicted xenograft volumes and their change over time in the Continuous-treatment <br/> 
dataset. Error bars represent the standard error of the means* <br>

### Figure 2C <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G2_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G2_p_light.png">
<img alt="Shows the light or dark mode version of the same plot."
  src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_methods_G2_p_dark.png">
</picture>
<br/>

*The figure shows the mean real and predicted xenograft volumes and their change over time in the Single-treatment <br/> 
dataset. Error bars represent the standard error of the means* <br>

The resulting graphs showed that the predicted volumes track the real volumes quite closely and based on the error bars <br/>
there seems to be no significant difference between the predicted and real tumor volumes (although a proper significance test <br/>
was not done to confirm this). Nonetheless, it is important to note that it seems that applying one of the correction methods <br/>
seems to increase the prediction precision. <br>

<br>

## The effect of correction on prediction precision <br>

In order to see if the previous observations hold true, and if the correction methods have any tangible effect on the <br/>
precision of the predictions, I calculated the Mean Absolute Percentage Error (MAPE) for all three XenoVol runs, and <br/>
compared them. Figure 3A-C depicts the visualized MAPE values of the estimations with and without correction. <br>

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

*The figure shows the MAPE values of XenoVol ran with and without correction on the Continuous-treatment dataset.* <br>

### Figure 3C <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_G2_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_G2_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/XV_method_MAPE_G2_p_dark.png">
</picture>
<br/>

*The figure shows the MAPE values of XenoVol ran with and without correction on the Single-treatment dataset.* <br>

The resulting graphs suggest that both correction methods considerably reduce prediction error compared to the base method. <br/>
However, there appears to be no significant difference between their performances (note that no formal test of significance <br/>
was conducted). Nonetheless, of the two correction methods, 'linear interpolation' seems to have a slightly better effect <br/>
on reducing prediction error, despite some limitations. First, it is limited to samples containing no sporadic NAs (or to <br/>
sample sets where these samples were removed) and can only be used if at least two reference measurements are available. <br/>
Meanwhile, 'mean correction' appears to reduce prediction errors to a lesser extent, but it is a more robust method that <br/>
can be used even if only a single reference measurement is available and in the presence of samples with sporadic NAs. <br>

Overall, both correction methods can noticably reduce prediction error and increase volume prediction accuracy. <br/>
As there is no significant difference between the two methods, their use should be left to the user's discretion. <br/>
Nonetheless, in general, the use of 'mean correction' is advised unless conditions specifically allow for the use <br/>
of 'linear interpolation'. <br>

<br>

# XenoVol against the world... <br>
<br>

## Predicted mean xenograft volumes across various methods

After exploring the effects of the afroementioned correction methods on volume prediction and prediction accuracy I <br/>
wanted to explore, how XenoVol stacks up againts other, commonly used volume prediction methods mentioned by Sápi et al. (2015).<br/>
In order to do that I calculated the mean tumor volumes of all collection time points, for all three generated datasets <br/>
using both XenoVol and multiple other formulas. Then, I overlaid the predicted tumor volumes to see how closely the <br/>
predicted volumes follow the real xenograft volumes just like I did in the previous section.  <br/>
Figure 4A-C depicts the visualized rpredicted mean tumor volumes across various methods. <br>

### Figure 4A <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_ctrl_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_ctrl_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_ctrl_p_dark.png">
</picture>
<br/>

*The figure shows the mean predicted xenograft volumes across various prediction methods overlaid, and their change <br/>
over time in the Control dataset. Error bars represent the standard error of the means.* <br>

### Figure 4B <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_G1_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_G1_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_G1_p_dark.png">
</picture>
<br/>

*The figure shows the mean predicted xenograft volumes across various prediction methods overlaid, and their change <br/>
over time in the Continuous-treatment dataset. Error bars represent the standard error of the means.* <br>

### Figure 4B <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_G2_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_G2_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_G2_p_dark.png">
</picture>
<br/>

*The figure shows the mean predicted xenograft volumes across various prediction methods overlaid, and their change <br/>
over time in the Single-treatment dataset. Error bars represent the standard error of the means.* <br>

The resulting graphs showed that as before, the volumes predicted by XenoVol track the real volumes quite closely, however, <br/>
the volumes predicted by the other tested methods seem to differ from the real volumes, and based on the error bars, <br/>
this difference seems to be significant (although a proper significance test was not done to confirm this).  <br/>
Nonetheless, it is important to note that all methods preserve the trend and trajectory of the volume changes. <br>

<br>

## Prediction precision across various methods <br>

Although the previous paragraph suggests high prediction error for the comparison models, I wanted to further explore <br/>
how the precision of the different methods compares. Therefore, as before, I calculated the MAPE for all <br/>
tested methods, and compared them. Figure 5A-C depicts the visualized MAPE values for the estimations with the various methods. <br>

### Figure 5A <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_MAPE_ctrl_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_MAPE_ctrl_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_MAPE_ctrl_p_dark.png">
</picture>
<br/>

*The figure shows the MAPE values of all tested prediction methods on the Control dataset.* <br>

### Figure 5B <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_MAPE_G1_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_MAPE_G1_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_MAPE_G1_p_dark.png">
</picture>
<br/>

*The figure shows the MAPE values of all tested prediction methods on the Continuous-treatment dataset.* <br>

### Figure 5C <br/>

<picture>
<source media="(prefers-color-scheme: dark)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_MAPE_G2_p_dark.png">
<source media="(prefers-color-scheme: light)" srcset="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_MAPE_G2_p_light.png">
<img alt="Shows the light or dark mode version of the same plot." src="https://github.com/UBE2C/Xenograft_volume_estimation/blob/main/Performance_plots/all_methods_MAPE_G2_p_dark.png">
</picture>
<br/>

*The figure shows the MAPE values of all tested prediction methods on the Single-treatment dataset.* <br>

As expected, the resulting figures indicate that all tested prediction methods have a higher prediction error <br/>
compared to XenoVol. An interesting note is that the Xenograft-model formula shows a more erratic behavior in the second <br/>
dataset, while the other two methods (Ellipsoid formula and Hemi-ellipsoid formula) show a very consistent error rate <br/>
in all three samples. <br>

Taken together, the results show that XenoVol outperforms the other tested methods in volume prediction accuracy. <br>

<br>

# References <br>

Sápi J, Kovács L, Drexler DA, Kocsis P, Gajári D, et al. (2015) Tumor Volume Estimation and Quasi-Continuous Administration for Most Effective Bevacizumab Therapy. <br>
PLOS ONE 10(11): e0142190. [https://doi.org/10.1371/journal.pone.0142190](https://doi.org/10.1371/journal.pone.0142190) <br>

Protocol Online (2005). [Xenograft tumor model protocol](https://www.protocol-online.org/prot/Protocols/Xenograft-Tumor-Model-Protocol-3810.html). <br>

Tomayko MM, Reynolds CP. Determination of subcutaneous tumor size in athymic (nude) mice. Cancer Chemother Pharmacol. 1989;24(3):148-54. doi: [10.1007/BF00300234](https://link.springer.com/article/10.1007/BF00300234). PMID: [2544306](https://pubmed.ncbi.nlm.nih.gov/2544306/). <br>

Sápi J, Drexler DA, Sápi Z, Kovács L (2014) Identification of C38 colon adenocarcinoma growth under bevacizumab therapy and without therapy. In: [CINTI 2014—15th IEEE International Symposium on Computational Intelligence and Informatics. pp. 443–448](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=88f75305a523a145a71f115626e5f85e68da7dd5). Budapest, Hungary. <br>

Heitjan DF, Manni A, Santen RJ. Statistical analysis of in vivo tumor growth experiments. [Cancer Res. 1993 Dec 15;53(24):6042-50](https://aacrjournals.org/cancerres/article/53/24/6042/499640/Statistical-Analysis-of-in-Vivo-Tumor-Growth). PMID: [8261420](https://pubmed.ncbi.nlm.nih.gov/8261420/).
