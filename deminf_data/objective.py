import os
import numpy as np
import moments
from .load_data import DemInfData
from .transform import get_transform_by_name, identical


def wrap_dem_inf_function(model_func, afs_data):
    """
    Returns function that will return log-likelihood of model and data.

    :param model_func: Function that returns demographic model.
    :param afs_data: Data of AFS.
    """
    def wrapper(params):
        ns = afs_data.sample_sizes
        afs_model = model_func(params, ns)
        ll = moments.Inference.ll_multinom(afs_model, afs_data)
        return ll
    return wrapper


class Objective(object):
    """
    Class for keeping objective function.

    :param obj_function: Objective function.
    :type obj_function: func
    :param n_params: Number of parameters of function.
    :type n_params: int
    :param lower_bound: Lower bound of parameters.
    :type lower_bound: list
    :param upper_bound: Upper bound of parameters.
    :type upper_bound: list
    :param transform: Transformation of parameter space.
    :type transform: func
    :param inv_transform: Inverse transformation of parameter space.
    :type inv_transform: func
    :param negate: If True then __call__ will return value of -obj_function(p).
    :type negate: bool
    :param name: Name for objective
    :type name: str

    :note: Lower and upper bounds are given in usual parameter space.\
           They are transformed for new parameter space during initialization.
    """
    def __init__(self, obj_function, n_params,
                 par_labels, lower_bound, upper_bound,
                 transform=identical, inv_transform=identical, negate=False,
                 name=None):
        self.obj_function = obj_function
        self.n_params = n_params
        self.transform = transform
        self.inv_transform = inv_transform
        self.par_labels = list(par_labels)
        self.lower_bound = np.array(lower_bound)
        self.lower_bound = self.transform(self.lower_bound)
        # check if there was 0 in lower bound and log-transform was applied
        self.lower_bound[self.lower_bound == -np.inf] = np.log(1e-15)
        self.upper_bound = np.array(upper_bound)
        self.upper_bound = self.transform(self.upper_bound)
        self.negate = negate
        self.name = name

        # checks
        assert len(self.lower_bound) == self.n_params
        assert len(self.upper_bound) == self.n_params
        assert len(self.par_labels) == self.n_params

    @staticmethod
    def from_name(name, negate=False, type_of_transform=None):
        """
        Creates objective function for data set with given name. Objective
        function will take values of parameters and return value of
        log-likelihood between demographic history and data.

        Transform could be one of three types:

        - `None` means identical or no transform
        - `logarithm` for logarithm transformation of parameter space
        - `custom_logarithm` will apply log for several parameters\
          (all except migrations)

        :param name: Name of data set.
        :type name: str
        :param negate: If True then -log-likelihood will be returned in
                       __call__ method.
        :type negate: bool
        :param type_of_transform: Name of transformation of parameter space.
        :type type_of_transform: str
        """
        obj_data = DemInfData.from_name(name)
        obj_function = wrap_dem_inf_function(obj_data.model_func,
                                             obj_data.afs_data)
        transform, inv_transform = get_transform_by_name(type_of_transform,
                                                         obj_data.par_labels)
        obj = Objective(obj_function=obj_function,
                        n_params=obj_data.number_of_parameters,
                        par_labels=obj_data.par_labels,
                        lower_bound=obj_data.lower_bound,
                        upper_bound=obj_data.upper_bound,
                        transform=transform,
                        inv_transform=inv_transform,
                        negate=negate,
                        name=name)
        return obj

    def __call__(self, params):
        """
        Call of objective function. Values of parameters are in transformed
        space so they are transformed back and objective function is called
        from them.
        """
        params = np.array(params)
        assert len(params) == self.n_params
        assert np.all(self.lower_bound <= params)
        assert np.all(params <= self.upper_bound)

        sign = -1 if self.negate else 1
        transf_params = self.inv_transform(params)
        return sign * self.obj_function(transf_params)

    def get_gadma_variables(self):
        """
        Returns list of objects :class:`gadma.Variable` that are variables
        in GADMA.
        """
        import gadma
        variables = []
        for label, lb, ub in zip(self.par_labels,
                                 self.lower_bound,
                                 self.upper_bound):
            if label.lower().startswith("n"):
               cls = gadma.PopulationSizeVariable
            elif label.lower().startswith("t"):
                cls = gadma.TimeVariable
            elif label.lower().startswith("m"):
                cls = gadma.MigrationVariable
            elif label.lower().startswith("s") or label.lower().startswith("f"):
                cls = gadma.FractionVariable
            else:
                raise ValueError(label)
            variables.append(cls(name=label, domain=[lb, ub]))
        return variables
