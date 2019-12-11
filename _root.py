# encoding: utf-8

import collections as col
import numpy as np
import _minpack


class OptimizationError(Exception):
    pass


def message(status, dimension):
    return {
        0: "Improper input parameters.",
        1: "Algorithm estimates that the relative error between x and the"
           " solution is at most tol.",
        2: "Number of calls to fcn has reached or exceeded {}."
           .format(200*(dimension + 1)),
        3: "Parameter tol is too small. No further improvement in the"
           " approximate solution x is possible.",
        4: "Iteration is not making good progress.",
    }[status]

OptimizationResult = col.namedtuple("OptimizationResult",
                                    "x success status message fun")


def root(fun, x0, args=(), tol=1e-10):
    """Calculate roots of a given function by using MINPACK's HYBRD algorithm.

    Args:
        fun (callback): the user-supplied subroutine which calculates the
            functions. It should be declared as follows.::

                def fun(x): return some_calculation(x)

            The value of flag should not be changed by fun unless the user wants
            to terminate execution of root. In this case flag must be
            set to a negative integer.
        x0 (numpy array): an array of the initial estimate of the solution
            vector.
        args (tuple, optional): Extra arguments passed to fun
        tol (float, optional): nonnegative input variable. Termination
            occurs when the algorithm estimates that the relative error between
            x and the solution is at most tol.

    Returns:
        OptimizationResult: containing members
            * x: result
            * success: whether the optimization was successful or not
            * status: integer flag
            * message: corresponding to status
            * fun: result of fun evaluated at the solution

    Raises:
        OptimizationError: This error is thrown if any error occurs in the
            HYBRD1 routine
    """

    def minpack_fun(x, flag):
        fvec = fun(x, *args)
        return flag, fvec

    _minpack.fun = minpack_fun

    try:
        x, status = _minpack.root(x0, tol)
    except Exception as e:
        raise OptimizationError(e)

    xfun = fun(x, *args)
    success = status == 1 or all(np.abs(val) < tol for val in xfun)

    return OptimizationResult(
        x, success, status, message(status, len(x0)), xfun
    )
