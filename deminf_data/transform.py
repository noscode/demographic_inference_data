import numpy as np
from functools import partial


def identical(params):
    """
    Identical transform.

    :param params: Parameters for transform.
    """
    return params


def logarithm(params):
    """
    Usual np.log with additional fix when params are np.array of dtype=object.
    """
    return np.log(params.astype(float))


def exponential(params):
    """
    Usual np.exp with additional fix when params are np.array of dtype=object.
    """
    return np.exp(params.astype(float))


def custom_logarithm(par_labels, params):
    """
    Custom transform: logarithm for all parameters except migrations.

    :param params: Parameters for transform.
    """
    transf_params = np.array(params)
    par_labels = np.array(par_labels)
    not_migs = (par_labels != 'm')
    transf_params[not_migs] = np.log(transf_params[not_migs].astype(float))
    return transf_params


def inv_custom_logarithm(par_labels, params):
    """
    Invariant transform of custom transform. Exponent is applied for all
    parameters except migrations.

    :param params: Parameters for transform.
    """

    transf_params = np.array(params)
    par_labels = np.array(par_labels)
    not_migs = (par_labels != 'm')
    transf_params[not_migs] = np.exp(transf_params[not_migs].astype(float))
    return transf_params


def get_transform_by_name(name, par_labels):
    """
    Returns transformation and inverse transformation as functions from name.

    :param name: name of the transformation.

    returns: tuple of transform and inverse transform
    """
    if name == None or name.lower() == "none":
        return identical, identical
    elif name == "logarithm":
        return logarithm, exponential
    elif name == "custom logarithm":
        return (partial(custom_logarithm, par_labels),
                partial(inv_custom_logarithm, par_labels))
    else:
        raise ValueError(f"Unknown name of the transformation: {name}")
