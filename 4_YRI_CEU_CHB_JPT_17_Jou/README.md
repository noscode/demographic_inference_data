## 4_YRI_CEU_CHB_JPT_17_Jou


| Number of populations | Number of parameters | Max log likelihood | Size of spectrum |
| --- | --- | --- | --- |
| 4 | 17 | -57730.966 | 40x40x40x40 |


### Model Description

Demographic model for four modern human populations: YRI, CEU, CHB and JPT. Model and data from Jouganous et al. 2019. Model with sudden growth of ancestral population size, followed by split into population YRI and common population of CEU, CHB and JPT, which experience bottleneck and split with exponential recovery of all populations - first formation of CEU population followed by split of ancestal population of CHB and JPT. Migrations between populations are symmetrical.

### Plots

Schematic model plot:

<img src="model_plot.png" height="500" />

Simulated allele frequency spectrum (projections):

<img src="fs_plot_projections.png" />


### Optimal parameter values

| Parameter | Value | Description |
| --- | --- | --- |
| `nuAf` | 2.101 | The ancestral population size after sudden growth and size of YRI population. |
| `nuB` | 0.251 | The bottleneck size of CEU+CHB common population. |
| `nuEu0` | 0.222 | The bottleneck size for CEU population. |
| `nuEu` | 2.809 | The final size of CEU population after exponential growth. |
| `nuAs0` | 0.090 | The bottleneck size for CHB population. |
| `nuAs` | 5.548 | The final size of CHB population after exponential growth. |
| `nuJp0` | 0.388 | The bottleneck size for JPT population. |
| `nuJp` | 20.731 | The final size of JPT population after exponential growth. |
| `mAfB` | 3.794 | The scaled symmetric migration rate between YRI and CEU+CHB populations. |
| `mAfEu` | 0.257 | The scaled symmetric migration rate between YRI and CEU populations. |
| `mAfAs` | 0.126 | The scaled symmetric migration rate between YRI and CHB populations. |
| `mEuAs` | 1.073 | The scaled symmetric migration rate between CEU and CHB populations. |
| `mChJp` | 0.745 | The scaled symmetric migration rate between CHB and JPT populations. |
| `TAf` | 0.363 | The scaled time between ancestral population growth and first split. |
| `TB` | 0.111 | The time between the first split and second. Time of CEU+CHB+JPT population existence. |
| `TEuAs` | 0.056 | The time between second split and third split. Time of CHB+JPT population existence. |
| `TEuAs` | 0.056 | The time between third split and present. |

