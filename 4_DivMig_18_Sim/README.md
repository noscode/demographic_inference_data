## 4_DivMig_18_Sim


| Number of populations | Number of parameters | Max log likelihood | Size of spectrum |
| --- | --- | --- | --- |
| 4 | 18 | -16452.849 | 10x10x10x10 |


### Model Description

Four population demographic history with 18 parameters. Ancestral population of constant size was split (T1 + T2 + T3) time ago into two new populations. First (1) population has constant size till nowdays. The second formed population turned out to be a common population to another three populations: (T2 + T3) time ago it splits and formed so-called population 2 and T3 time ago the second formed population divided into population 3 and population 4. All migrations between populations were symmetrical.

### Plots

Schematic model plot:

<img src="model_plot.png" height="500" />

Simulated allele frequency spectrum (projections):

<img src="fs_plot_projections.png" />


### Optimal parameter values

| Parameter | Value | Description |
| --- | --- | --- |
| `nu1` | 1.500 | Size of population 1 after split of ancestral population. |
| `nu234` | 0.800 | Size of common ancestor population of populations 2, 3 and 4 after split of ancestral population. |
| `nu2` | 1.000 | Size of population 2 after split of common ancestor population of populations 2, 3 and 4. |
| `nu34` | 0.500 | Size of common ancestor population of populations 3 and 4 after division of population 2 from their common ancestor population. |
| `nu3` | 0.200 | Size of population 3. |
| `nu4` | 0.300 | Size of population 4. |
| `m1_234` | 5.000 | Symmetric migration rate between population 1 and common ancestor population of 2, 3 and 4 populations (between first and second splits). |
| `m1_2` | 0.500 | Symmetric migration rate between population 1 and population 2 (between second and third splits). |
| `m1_34` | 1.000 | Symmetric migration rate between population 1 and common ancestor of population 3 and population 4 (between second and third splits). |
| `m2_34` | 3.000 | Symmetric migration rate between population 2 and common ancestor of population 3 and population 4 (between second and third splits). |
| `m1_3` | 0.400 | Symmetric migration rate between population 1 and population 3 (after the third split). |
| `m1_4` | 0.300 | Symmetric migration rate between population 1 and population 4 (after the third split). |
| `m2_3` | 1.200 | Symmetric migration rate between population 2 and population 3 (after the third split). |
| `m2_4` | 1.300 | Symmetric migration rate between population 2 and population 4 (after the third split). |
| `m3_4` | 2.000 | Symmetric migration rate between population 3 and population 4 (after the third split). |
| `T1` | 0.100 | Time between ancestral population split (population 1 formation) and next split. |
| `T2` | 0.150 | Time between ancestral population of populations 2, 3 and 4 split (population 2 formation) and next split. |
| `T3` | 0.050 | Time of ancestral population of populations 3 and 4 split which have led to formations of population 3 and population 4. |

