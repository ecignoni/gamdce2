module c_interface
    use iso_c_binding, only: c_int, c_double
    use mod_gamdce2, only: reweight_ce2_1d, reweight_ce2_2d

    implicit none

contains

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

    subroutine c_reweight_ce2_2d(n_samples, cv1, cv2, bias, n_bins1, n_bins2, &
            fes, fes_min1, fes_max1, fes_min2, fes_max2, beta) bind(c)
        implicit none
        integer(c_int), intent(in) :: n_samples
        real(c_double), intent(in), dimension(n_samples) :: cv1, cv2
        real(c_double), intent(in), dimension(n_samples) :: bias
        integer(c_int), intent(in) :: n_bins1, n_bins2
        real(c_double), intent(inout), dimension(n_bins1, n_bins2) :: fes
        real(c_double), intent(in) :: fes_min1, fes_max1, fes_min2, fes_max2
        real(c_double), intent(in) :: beta
        call reweight_ce2_2d(n_samples, cv1, cv2, bias, n_bins1, n_bins2, fes, fes_min1, fes_max1, fes_min2, fes_max2, beta)
    end subroutine c_reweight_ce2_2d

end module c_interface
