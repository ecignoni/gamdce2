import cython
import numpy as np
cimport numpy as np

cdef extern from "c_interface.h":
    void c_histogram1d(double *arr, int *arr_len, int *n_bins, double *hist, double *hist_min, double *hist_max)
    void c_reweight_ce2_1d(int *n_samples, double *cv, double *bias, int *n_bins, double *fes, double *fes_min, double *fes_max, double *beta)

def histogram1d(double[::1] arr, int arr_len, int n_bins, double[::1] hist, double hist_min, double hist_max):
    c_histogram1d(&arr[0], &arr_len, &n_bins, &hist[0], &hist_min, &hist_max)
    return hist.base

def reweight_ce2_1d(int n_samples, double[::1] cv, double[::1] bias, int n_bins, double[::1] fes, double fes_min, double fes_max, double beta):
    c_reweight_ce2_1d(&n_samples, &cv[0], &bias[0], &n_bins, &fes[0], &fes_min, &fes_max, &beta)
    return fes.base
