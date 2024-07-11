module c_interface
    use iso_c_binding, only: c_int, c_double
    use mod_gamdce2, only: histogram1d, reweight_ce2_1d

    implicit none

contains

    subroutine c_histogram1d(arr, arr_len, n_bins, hist, hist_min, hist_max) bind(c)
        implicit none
        integer(c_int), intent(in) :: arr_len
        real(c_double), intent(in), dimension(arr_len) :: arr
        integer(c_int), intent(in) :: n_bins
        real(c_double), intent(inout), dimension(n_bins) :: hist
        real(c_double), intent(in) :: hist_min
        real(c_double), intent(in) :: hist_max
        call histogram1d(arr, arr_len, n_bins, hist, hist_min, hist_max)
    end subroutine c_histogram1d

    subroutine c_reweight_ce2_1d(n_samples, cv, bias, n_bins, fes, fes_min, fes_max, beta) bind(c)
        implicit none
        integer(c_int), intent(in) :: n_samples
        real(c_double), intent(in), dimension(n_samples) :: cv
        real(c_double), intent(in), dimension(n_samples) :: bias
        integer(c_int), intent(in) :: n_bins
        real(c_double), intent(inout), dimension(n_bins) :: fes
        real(c_double), intent(in) :: fes_min
        real(c_double), intent(in) :: fes_max
        real(c_double), intent(in) :: beta
        call reweight_ce2_1d(n_samples, cv, bias, n_bins, fes, fes_min, fes_max, beta)
    end subroutine c_reweight_ce2_1d

end module c_interface
