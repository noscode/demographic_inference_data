import os
import sys
import numpy as np
import moments
import importlib


def load_module(dirname, filename, base_dir="."):
    """
    Loads python file as module. Enables access to all fields in file.

    :param dirname: Name of the directory with file in current base dir.
    :type dirname: str
    :param filename: Name of file to load.
    :type filename: str
    :param base_dir: Base directory of directories.
    :type base_dir: str

    :note: The problem is about names of loaded modules, so dirname should be\
           set independently in order to add variety in names.
    """
    save_dir = os.path.abspath(".")

    base_dir = os.path.abspath(base_dir)
    os.chdir(os.path.abspath(dirname))
    sys.path.append(base_dir)
    sys.path.append(os.path.join(base_dir, dirname))

    module_name = os.path.join(dirname, filename)
    module_name = module_name.replace('/', '.').rstrip('.py')
    module_name = module_name.strip(".")
    module = importlib.import_module(module_name)

    # this modules could import models so we should remove this imports
    if "demographic_model" in sys.modules:
        del sys.modules["demographic_model"]

    sys.path = sys.path[:-2]
    os.chdir(save_dir)
    return module


def load_module_from_path(path_to_file):
    """
    Loads python file as module from its path.
    Enables access to all fields in file.

    :param path_to_file: Path to directory with file.
    :type path_to_file: str
    """
    path, filename = os.path.split(path_to_file)
    return load_module(path, filename)


class DemInfData(object):
    """
    Keeps object with case of data for demographic inference.

    :param name: Name of data set, consists of the following parts:
                 <number of populations>_<short description>_<number of
                 parameters>_<origin (simulation or from paper)>
    :param afs_data: Allele frequency spectrum (AFS).
    :param sequence_length: Length of sequence that was used for AFS.
    :param model_func: Function of demographic model.
    :param mutation_rate: Mutation rate in demographic model.
    :param pop_labels: Labels of populations.
    :param par_labels: Labels of the demographic parameters.
    :param lower_bound: Lower bound on values of demographic parameters.
    :param upper_bound: Lower bound on values of demographic parameters.
    :param popt: Known optimal values of the parameters.
    :param max_ll: Known values of maximal log-likelihood (on `popt`).
    """
    def __init__(self,
                 name,
                 afs_data,
                 sequence_length,
                 model_func,
                 mutation_rate,
                 pop_labels,
                 par_labels,
                 lower_bound,
                 upper_bound,
                 popt,
                 max_ll):
        self.name = name
        self.afs_data = afs_data
        self.sample_sizes = afs_data.sample_sizes
        self.model_func = model_func
        self.pop_labels = pop_labels
        self.number_of_populations = len(self.pop_labels)
        self.par_labels = par_labels
        self.number_of_parameters = len(self.par_labels)
        self.lower_bound = np.array(lower_bound)
        self.upper_bound = np.array(upper_bound)
        self.popt = popt
        self.max_ll = max_ll
        self.mutation_rate = mutation_rate
        self.sequence_length = sequence_length

        # checks
        assert self.number_of_populations == int(self.name.split("_")[0])
        assert self.number_of_parameters == len(self.lower_bound)
        assert self.number_of_parameters == len(self.upper_bound)
        assert self.number_of_parameters == len(self.popt)
        assert np.all(self.lower_bound <= self.upper_bound)

    @staticmethod
    def from_name(name):
        """
        Loads all required data from its name.
        """
        base_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.join(base_dir, "..")
        data_dir = os.path.join(base_dir, name)
        data_dir = os.path.abspath(data_dir)

        afs_data = moments.Spectrum.from_file(os.path.join(data_dir,
                                                        "fs_data.fs"))
        model_func = load_module(name,
                                 "demographic_model.py",
                                 base_dir).model_func
        model_info = load_module(name, "main_script.py", base_dir)
        
        obj = DemInfData(name=name,
                         afs_data=afs_data,
                         sequence_length=model_info.L,
                         model_func=model_func,
                         mutation_rate=model_info.mu,
                         pop_labels=model_info.pop_labels,
                         par_labels=model_info.par_labels,
                         lower_bound=model_info.lower_bound,
                         upper_bound=model_info.upper_bound,
                         popt=model_info.popt,
                         max_ll=model_info.max_ll)
        return obj

    def get_known_maximum_likelihood(self):
        """
        Returns maximum value of log-likelihood. It could be exact value of
        optimum or known estimation.
        """
        return self.max_ll

    def is_maximum_likelihood_exact(self):
        """
        Returns if value of log-likelihood is exact value or estimation.
        It is exactly value if origin of data set is simulation and estimation
        otherwise.
        """
        return self.name.lower().endswith("sim")
