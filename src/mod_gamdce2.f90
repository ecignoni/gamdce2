module mod_gamdce2
    use mod_types, only: wp=>dp
    implicit none

contains

    subroutine histogram1d(arr, arr_len, n_bins, hist, hist_min, hist_max)
        implicit none
        integer, intent(in) :: arr_len
        real(wp), intent(in), dimension(arr_len) :: arr
        integer, intent(in) :: n_bins
        real(wp), intent(inout), dimension(n_bins) :: hist
        real(wp), intent(in) :: hist_min
        real(wp), intent(in) :: hist_max

        real(wp) :: width
        integer :: i
        integer :: bin_idx

        width = (hist_max - hist_min) / n_bins
 
        do i = 1, arr_len
            bin_idx = int((arr(i)-hist_min)/width) + 1
            hist(bin_idx) = hist(bin_idx) + 1
        end do
    end subroutine histogram1d

    subroutine reweight_ce2_1d(n_samples, cv, bias, n_bins, fes, fes_min, fes_max, beta)
        implicit none

        ! number of samples/frames
        integer, intent(in) :: n_samples
        ! collective variable
        real(wp), intent(in), dimension(n_samples) :: cv
        ! bias
        real(wp), intent(in), dimension(n_samples) :: bias
        ! number of bins
        integer, intent(in) :: n_bins
        ! fes
        real(wp), intent(inout), dimension(n_bins) :: fes
        ! minimum value of the fes
        real(wp), intent(in) :: fes_min
        ! maximum value of the fes
        real(wp), intent(in) :: fes_max
        ! temperature
        real(wp), intent(in) :: beta
        
        ! temporary arrays storing the mean and mean squared
        ! of the bias on the fes grid
        real(wp), dimension(n_bins) :: bias_mean
        real(wp), dimension(n_bins) :: bias_mean2
        integer, dimension(n_bins) :: bin_count

        ! loop variable
        integer :: i
        ! variable of bin index to update fes
        integer :: bin
        real(wp) :: bin_width
        real(wp) :: beta_inv

        ! initialize arrays
        fes = 0.0
        bias_mean = 0.0
        bias_mean2 = 0.0
        bin_count = 0

        bin_width = (fes_max - fes_min) / n_bins
        beta_inv = 1.d0 / beta

        ! loop over the samples, build the fes histogram (inside fes)
        ! and collect mean and mean squared of the bias
        do i = 1, n_samples
            bin = int((cv(i) - fes_min) / bin_width) + 1
            fes(bin) = fes(bin) + 1.0
            bias_mean(bin) = bias_mean(bin) + bias(i)
            bias_mean2(bin) = bias_mean2(bin) + bias(i) * bias(i)
            bin_count(bin) = bin_count(bin) + 1
        end do

        ! - logarithm of the histogram
        do i = 1, n_bins
            if (bin_count(i) > 5) then
                fes(i) = - beta_inv * log(fes(i))
                ! mean of the bias
                bias_mean(i) = bias_mean(i) / bin_count(i)
                ! variance of the bias
                bias_mean2(i) = (bias_mean2(i) / bin_count(i)) - &
                    (bin_count(i) / (bin_count(i) - 1.0d0)) * bias_mean(i) * bias_mean(i)
            else
                fes(i) = 0.0d0
                bias_mean(i) = 0.0d0
                bias_mean2(i) = 0.0d0
            end if
        end do

        fes = fes - bias_mean - 0.5d0 * beta * bias_mean2

    end subroutine reweight_ce2_1d

end module mod_gamdce2
