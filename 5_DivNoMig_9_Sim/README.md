## 5_DivNoMig_9_Sim


| Number of populations | Number of parameters | Max log likelihood | Size of spectrum |
| --- | --- | --- | --- |
| 5 | 9 | -55340.517 | 10x10x10x10x10 |


### Model Description

Simple demographic history of five populations without migrations. Model have 9 parameters. Each population number i has constant size of nui. Ancestral population splits (T1+T2+T3+T4) time ago to population 1 and 2. Then (T2+T3+T4) time ago population 2 splits into population 2 and 3. (T3 + T4) time ago population 3 split to populations 3 and 4. And finally T4 time ago population 4 split in populations 4 and 5.

### Plots

Schematic model plot:

<img src="model_plot.png" height="500" />

Simulated allele frequency spectrum (projections):

<img src="fs_plot_projections.png" />


### Optimal parameter values

| Parameter | Value | Description |
| --- | --- | --- |
| `nu1` | 1.000 | Size of population 1. |
| `nu2` | 2.000 | Size of population 2. |
| `nu3` | 1.500 | Size of population 3 after split from population 2. |
| `nu4` | 1.000 | Size of population 4 after split from population 3. |
| `nu5` | 0.500 | Size of population 5 after split from population 4. |
| `T1` | 0.050 | Time between ancestral population split and second split. |
| `T2` | 0.100 | Time between second and third splits. |
| `T3` | 0.150 | Time between third and fourth splits. |
| `T4` | 0.050 | Time of fourth split. |

