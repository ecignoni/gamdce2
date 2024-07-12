import cython
import numpy as np
cimport numpy as np


# declare external C functions
cdef extern from "c_interface.h":
    void c_reweight_ce2_1d(int *n_samples, double *cv, double *bias, int *n_bins, double *fes, double *fes_min, double *fes_max, double *beta)
    void c_reweight_ce2_2d(int *n_samples, double *cv1, double *cv2, double *bias, int *n_bins1, int *n_bins2,
                           double *fes, double *fes_min1, double *fes_max1, double *fes_min2, double *fes_max2, double *beta)

# python interface
def reweight_1d(
    int n_samples,
    double[::1] cv,
    double[::1] bias,
    int n_bins,
    double fes_min,
    double fes_max,
    double beta
):
    """reweighting with 2nd order cumulant expansion - 1D

    Performs the reweighting using 2nd order Cumulant Expansion.

        F(s) = F*(s) - E[ΔV] - (1/2) β Var[ΔV]

             = - (1/β) log P*(s) - E[ΔV] - (1/2) β Var[ΔV]

    Parameters
    ----------
    n_samples: int
        number of samples/frames
    cv: np.ndarray, shape (n_samples,)
        collective variable
    bias: np.ndarray, shape (n_samples,)
        bias potential
    n_bins: int
        number of bins to dump the free energy
    fes_min: int
        minimum value of the free energy grid
    fes_max: int
        maximum value of the free energy grid
    beta: float
        thermodynamic beta, 1. / kBT

    Returns
    -------
    fes: np.ndarray, shape (n_bins,)
        reweighted free energy
    """
    if cv.ndim != 1 or cv.shape[0] != n_samples:
        raise ValueError("cv must be given as 1D array of length n_samples.")
    if bias.ndim != 1 or bias.shape[0] != n_samples:
        raise ValueError("bias must be given as 1D array of length n_samples.")
    if n_bins <= 0:
        raise ValueError("negative/zero number of bins.")
    if fes_min >= fes_max:
        raise ValueError("fes_min >= fes_max.")
    if beta <= 0.0:
        raise ValueError("beta is not a positive number.")

    _fes = np.empty((n_bins,), dtype=np.float64, order="F")
    cdef double[::1] fes = _fes

    c_reweight_ce2_1d(&n_samples, &cv[0], &bias[0], &n_bins, &fes[0], &fes_min, &fes_max, &beta)

    return np.ascontiguousarray(fes.base)

def reweight_2d(
    int n_samples,
    double[::1] cv1,
    double[::1] cv2,
    double[::1] bias,
    int n_bins1,
    int n_bins2,
    double fes_min1,
    double fes_max1,
    double fes_min2,
    double fes_max2,
    double beta
):
    """reweighting with 2nd order cumulant expansion - 2D

    Performs the reweighting using 2nd order Cumulant Expansion.

        F(s) = F*(s) - E[ΔV] - (1/2) β Var[ΔV]

             = - (1/β) log P*(s) - E[ΔV] - (1/2) β Var[ΔV]

    Parameters
    ----------
    n_samples: int
        number of samples/frames
    cv1: np.ndarray, shape (n_samples,)
        first collective variable
    cv2: np.ndarray, shape (n_samples,)
        second collective variable
    bias: np.ndarray, shape (n_samples,)
        bias potential
    n_bins1: int
        number of bins for the first cv
    n_bins2: int
        number of bins for the second cv
    fes_min1: int
        minimum value of the free energy grid for the first cv
    fes_max1: int
        maximum value of the free energy grid for the first cv
    fes_min2: int
        minimum value of the free energy grid for the second cv
    fes_max2: int
        maximum value of the free energy grid for the second cv
    beta: float
        thermodynamic beta, 1. / kBT

    Returns
    -------
    fes: np.ndarray, shape (n_bins1, n_bins2)
        reweighted free energy
    """
    if cv1.ndim != 1 or cv1.shape[0] != n_samples:
        raise ValueError("cv1 must be given as 1D array of length n_samples.")
    if cv2.ndim != 1 or cv2.shape[0] != n_samples:
        raise ValueError("cv2 must be given as 1D array of length n_samples.")
    if bias.ndim != 1 or bias.shape[0] != n_samples:
        raise ValueError("bias must be given as 1D array of length n_samples.")
    if n_bins1 <= 0:
        raise ValueError("negative/zero number of bins for cv1.")
    if n_bins2 <= 0:
        raise ValueError("negative/zero number of bins for cv2.")
    if fes_min1 >= fes_max1:
        raise ValueError("fes_min1 >= fes_max1.")
    if fes_min2 >= fes_max2:
        raise ValueError("fes_min2 >= fes_max2.")
    if beta <= 0.0:
        raise ValueError("beta is not a positive number.")

    _fes = np.empty((n_bins1, n_bins2), dtype=np.float64, order="F")
    cdef double[::1, :] fes = _fes

    c_reweight_ce2_2d(
        &n_samples,
        &cv1[0],
        &cv2[0],
        &bias[0],
        &n_bins1,
        &n_bins2,
        &fes[0, 0],
        &fes_min1,
        &fes_max1,
        &fes_min2,
        &fes_max2,
        &beta
    )
    # return C-contiguous
    return np.ascontiguousarray(fes.base)
