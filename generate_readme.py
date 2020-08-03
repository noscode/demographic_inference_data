import os
import sys
import importlib.util
import re


REPO = "https://bitbucket.org/noscode/demographic_inference_data/src/master/"

INITIAL_README = """# Data for demographic inference

This repo contains AFS data (simulated and real) for different kinds of
demographic inference, for example, it could be used for some algorithm
benchmark and so on.

Right now everything is for the python package
[moments](https://bitbucket.org/simongravel/moments/src/master).

## Structure of data folders

### Name of folder

All data is located in its folder each of which has special name:

`{N}_{Descr}_{M}_{Origin}`

where `N` is the **number of populations**, `Descr` stands for **simple
description** of the demographic model, `M` is the **number of parameters**
in model and `Origin` is three letters of data origin. For example:

`2_DivMig_5_Sim` means AFS data for 2 populations that was simulated (`Sim`)
under the demographic history with divergence and migration that has in total
5 parameters.

For real data `Origin` contains the first three letters of the first author of the
corresponding paper. For example:

`2_YRI_CEU_6_Gut` means AFS data and demographic model for YRI and CEU
populations from the **Gut**enkunst et al., 2009 paper.

### Folder structure

Each directory contains several files:

* demographic_model.py - file with demographic model code. Code is written
for [moments](https://bitbucket.org/simongravel/moments/src/master).

* `main_script.py` - script with all information about AFS data. If data is
simulated then during script run it will be saved to `fs_data.fs`.

* `fs_data.fs` - AFS data in fs dadi format.

* `model_plot.png` - schematic plot of the demographic model.

* `fs_plot.png` or `fs_plot_projections.png` - plots of AFS.

File `main.py` is a valid file with python code and contains following
information about data:

| Variable | Description |
| --- | --- |
| `n_pop` | Number of populations |
| `par_labels` | Labels of model parameters |
| `max_ll` | Maximum log composite likelihood |
| `popt` | Optimal values parameters |
| `Nanc` | Size of ancestral population |
| `lower_bound` | Lower bound for parameters |
| `upper_bound` | Upper bound for parameters |
| `ns` | Sample sizes of AFS |

### Units

Note that all parameters values are in genetic units (relative to size
`Nanc` of ancestral population). Size of populations are relative to `Nanc`,
time is in `2 * Nanc` generations abd migrations are in `1 / (2 * Nanc)`
units. For more information see [moments manual](https://bitbucket.org/
simongravel/moments/src/master/doc/manual/manual.pdf) section `5.2 Units`.

"""

def load_module(dirname, filename):
    save_dir = os.path.abspath(".")
    sys.path.append(dirname)
    module = importlib.import_module(
        os.path.join(dirname, filename).replace('/', '.').rstrip('.py'))
    if "demographic_model" in sys.modules:
        del sys.modules["demographic_model"]
    sys.path = sys.path[:-1]
    return module

def generate_model_info(dirname, working_dir=None):
    if working_dir is None:
        working_dir = dirname
    s = f"## {dirname}\n\n"

    dem_model_file = 'demographic_model.py'
    sim_file = 'main_script.py'

    model = load_module(dirname, dem_model_file)
    model_description = load_module(dirname,
                                    dem_model_file).model_func.__doc__
    sim_info = load_module(dirname, sim_file)

    # Fast info
    s += "\n| Number of populations | Number of parameters "\
         "| Size of spectrum |\n| --- | --- | --- |\n"
    sp_size = "x".join([str(x) for x in sim_info.ns])
    s += f"| {sim_info.n_pop} | {len(sim_info.par_labels)} | {sp_size} |\n\n"

    # Description
    s += "\n### Model Description\n\n"
    par_labels = sim_info.par_labels
    values = {par_name: value for par_name, value in zip(par_labels,
                                                         sim_info.popt)}
    split_ind = model_description.rfind("\n\n")
    descr = re.sub("(&!\n)\n(&!\n)", "", model_description[:split_ind])
    descr = " ".join([x.strip() for x in descr.split("\n ")]).strip()
    s += descr

    # Plots
    s += "\n### Plots\n\n"
    s += "Schematic model plot:\n\n"
    s += f'<img src="{working_dir}/model_plot.png" height="500" />\n\n'

    if sim_info.n_pop <= 3:
        s += "Simulated allele frequency spectrum:\n\n"
        s += f'<img src="{working_dir}/fs_plot.png" height="200" />\n\n'
    if sim_info.n_pop >= 3:
        s += "Simulated allele frequency spectrum (projections):\n\n"
        s += f'<img src="{working_dir}/fs_plot_projections.png" '\
             'height="200" />\n\n'

    # Parameters
    s += "\n### Optimal parameter values\n\n"
    s += "| Parameter | Value | Description |\n"
    s += "| --- | --- | --- |\n"
    for line in model_description[split_ind+1:].split(":param "):
        if line.strip() == "":
            continue
        spl_ind = line.find(":")
        par_name = line[:spl_ind]
        par_descr = line[spl_ind + 1:].strip()
        par_descr = " ".join([x.strip() for x in par_descr.split("\n ")])
        s += f"| `{par_name}` | {values[par_name]} | {par_descr} |\n"
    s += "\n"
    return s

sim_dirs = ['1_Bot_4_Sim', '2_DivMig_5_Sim', '3_DivMig_8_Sim',
            "4_DivMig_11_Sim"]

with open("README.md", "w") as f:
    f.write(INITIAL_README)
    f.write("# Simulated data\n\n")
    for data_dir in sim_dirs:
        info = generate_model_info(data_dir)
        f.write(info)
        with open(os.path.join(data_dir, "README.md"), "w") as loc_f:
            loc_f.write(info.replace(f"{data_dir}/", ""))
    
