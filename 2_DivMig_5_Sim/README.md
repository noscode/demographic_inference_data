## 2_DivMig_5_Sim


| Number of populations | Number of parameters | Size of spectrum |
| --- | --- | --- |
| 2 | 5 | 20x20 |


### Model Description

Simple two populations model. Ancestral population of constant size splits into two subpopulations of constant size with asymetrical migrations.
### Plots

Schematic model plot:

<img src="model_plot.png" height="500" />

Simulated allele frequency spectrum:

<img src="fs_plot.png" height="200" />


### Optimal parameter values

| Parameter | Value | Description |
| --- | --- | --- |
| `nu1` | 1.0 | Size of subpopulation 1 after split. |
| `nu2` | 0.1 | Size of subpopulation 2 after split. |
| `m12` | 5 | Migration rate from subpopulation 2 to subpopulation 1. |
| `m21` | 2.5 | Migration rate from subpopulation 1 to subpopulation 2. |
| `T` | 0.05 | Time of split. |

