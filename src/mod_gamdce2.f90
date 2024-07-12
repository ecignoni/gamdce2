module mod_gamdce2
    use mod_types, only: wp=>dp
    implicit none

contains

    subroutine reweight_ce2_1d(n_samples, cv, bias, n_bins, fes, fes_min, fes_max, beta)
        implicit none

        ! number of samples/frames
        integer, intent(in) :: n_samples
        ! collective variable
        real(wp), intent(in), dimension(n_samples) :: cv
        ! bias potential
        real(wp), intent(in), dimension(n_samples) :: bias
        ! number of bins
        integer, intent(in) :: n_bins
        ! fes
        real(wp), intent(inout), dimension(n_bins) :: fes
        ! minimum and maximum value of the fes
        real(wp), intent(in) :: fes_min, fes_max
        ! thermodynamic beta
        real(wp), intent(in) :: beta
        
        ! temporary arrays storing the mean and mean squared
        ! of the bias on the fes grid
        real(wp), dimension(n_bins) :: bias_mean
        real(wp), dimension(n_bins) :: bias_mean2
        ! bin count to perform the average inside each bin
        integer, dimension(n_bins) :: bin_count

        ! loop variable
        integer :: i
        ! variable of bin index to update fes
        integer :: bin
        ! width of the bins
        real(wp) :: bin_width
        ! inverse of beta
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
            ! bounds check
            if ((cv(i) < fes_min) .or. (cv(i) > fes_max)) then
                cycle
            end if
            bin = int((cv(i) - fes_min) / bin_width) + 1
            ! update
            fes(bin) = fes(bin) + 1.0d0
            bias_mean(bin) = bias_mean(bin) + bias(i)
            bias_mean2(bin) = bias_mean2(bin) + bias(i) * bias(i)
            bin_count(bin) = bin_count(bin) + 1
        end do

        ! - logarithm of the histogram
        do i = 1, n_bins
            ! only if we have some statistics
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

    subroutine reweight_ce2_2d(n_samples, cv1, cv2, bias, n_bins1, n_bins2, fes, fes_min1, fes_max1, fes_min2, fes_max2, beta)
        implicit none

        ! number of samples / frames
        integer, intent(in) :: n_samples
        ! collective variables
        real(wp), intent(in), dimension(n_samples) :: cv1, cv2
        ! bias potential
        real(wp), intent(in), dimension(n_samples) :: bias
        ! number of bins for first and second cv
        integer, intent(in) :: n_bins1, n_bins2
        ! fes
        real(wp), intent(inout), dimension(n_bins1, n_bins2) :: fes
        ! minimum and maximum values of the fes along the two cvs
        real(wp), intent(in) :: fes_min1, fes_max1, fes_min2, fes_max2
        ! thermodynamic beta
        real(wp), intent(in) :: beta

        ! temporary arrays storing the mean and mean squared
        ! of the bias on the fes grid
        real(wp), dimension(n_bins1, n_bins2) :: bias_mean
        real(wp), dimension(n_bins1, n_bins2) :: bias_mean2
        ! bin counts to perform averages inside each bin
        integer, dimension(n_bins1, n_bins2) :: bin_count

        ! loop variables
        integer :: i, j
        ! variables of bin index to update fes
        integer :: bin1, bin2
        ! width of the bins along the two dimensions
        real(wp) :: bin_width1, bin_width2
        ! inverse of beta
        real(wp) :: beta_inv

        ! initialize arrays
        fes = 0.0
        bias_mean = 0.0
        bias_mean2 = 0.0
        bin_count = 0

        bin_width1 = (fes_max1 - fes_min1) / n_bins1
        bin_width2 = (fes_max2 - fes_min2) / n_bins2
        beta_inv = 1.d0 / beta

        ! loop over the samples, build the fes histogram (inside fes)
        ! and collect mean and mean squared of the bias
        do i = 1, n_samples
            ! first cv out of bounds
            if ((cv1(i) < fes_min1) .or. (cv1(i) > fes_max1)) then
                cycle
            ! second cv out of bounds
            else if ((cv2(i) < fes_min2) .or. (cv2(i) > fes_max2)) then
                cycle
            end if
            bin1 = int((cv1(i) - fes_min1) / bin_width1) + 1
            bin2 = int((cv2(i) - fes_min2) / bin_width2) + 1
            ! update
            fes(bin1, bin2) = fes(bin1, bin2) + 1.0d0
            bias_mean(bin1, bin2) = bias_mean(bin1, bin2) + bias(i)
            bias_mean2(bin1, bin2) = bias_mean2(bin1, bin2) + bias(i) * bias(i)
            bin_count(bin1, bin2) = bin_count(bin1, bin2) + 1
        end do

        ! logarithm of the histogram
        do j = 1, n_bins2
            do i = 1, n_bins1

                ! only if we have some statistics
                if (bin_count(i, j) > 5) then
                    fes(i, j) = - beta_inv * log(fes(i, j))
                    ! mean of the bias
                    bias_mean(i, j) = bias_mean(i, j) / bin_count(i, j)
                    ! variance of the bias
                    bias_mean2(i, j) = (bias_mean2(i, j) / bin_count(i, j)) - &
                        (bin_count(i, j) / (bin_count(i, j) - 1.0d0)) * bias_mean(i, j) * bias_mean(i, j)
                else
                    fes(i, j) = 0.0d0
                    bias_mean(i, j) = 0.0d0
                    bias_mean2(i, j) = 0.0d0
                end if

            end do
        end do

        fes = fes - bias_mean - 0.5d0 * beta * bias_mean2

    end subroutine reweight_ce2_2d

end module mod_gamdce2
