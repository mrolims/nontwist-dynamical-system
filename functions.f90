module functions
    implicit none
    
    contains

    subroutine time_series(out, theta0, I0, params, N)
        implicit none
        real(8), intent(in) :: theta0, I0, params(:)
        integer, intent(in) :: N
        real(8), intent(out) :: out(:, :)
        ! --- Parameters --- !
        real(8), parameter :: pi = dacos(-1.0d0)
        ! --- Variables --- !
        integer :: j
        real(8) :: theta, I
        real(8) :: q_factor, E_r, v_parallel

        theta = theta0
        I = I0
        do j = 0, N
            ! Update I
            I = I + params(14) * dsin(2 * pi * theta)
            ! Functions
            q_factor = params(1) + params(2) * I ** 2 + params(3) * I ** 3
            E_r = params(4) * I + params(5) * dsqrt(dabs(I)) + params(6)
            v_parallel = params(7) + params(8) * dtanh(params(9) * I + params(10))
            ! Update theta
            theta = modulo(theta + params(13) * v_parallel * (params(11) / q_factor - params(12)) + params(15) * E_r / dsqrt(dabs(I)), 1.0d0)
            out(j, 1) = theta
            out(j, 2) = I
        end do
        
        return

    end subroutine time_series

    subroutine return_base_exp(b, e, num)
        implicit none
        integer, intent(out) :: b, e
        real(8), intent(in) :: num

        e = int(floor(dlog10(num)))
        b = int(num/10**(dble(e)))

        return

    end subroutine return_base_exp
    
end module functions