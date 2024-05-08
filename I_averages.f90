program I_averages
    ! I_averages.f90
    ! ==============
    ! 
    ! Calculates the average of the action and the square root of
    ! the squared averaged action as a function of time.
    !
    ! Usage
    ! -----
    ! This program takes on three parameters
    !
    !       ifx functions.f90 I_averages.f90 -o Iav.x
    !       ./Iav eps I_ini N
    !
    ! Author: Matheus Rolim Sales
    ! Email: matheusrolim95@gmail.com
    ! Last updated: 08/05/2024
    use functions
    implicit none
    ! --- Parameters --- !
    integer, parameter :: n_ic = int(1e6)
    real(8), parameter :: pi = dacos(-1.0d0)
    ! --- Variables --- !
    integer :: j
    integer :: N
    integer :: eN, bN, en_ic, bn_ic, eI_ini, bI_ini
    real(8) :: eps
    real(8), allocatable :: Imean(:), Irms(:)
    real(8) :: I_ini
    real(8) :: theta(n_ic), I(n_ic)
    real(8) :: params(15)
    real(8) :: q_factor(n_ic), E_r(n_ic), v_parallel(n_ic)
    real(8) :: S(n_ic)
    character :: path*10, datafile*250, arg*15

    path = "Data"

    call getarg(1, arg)
    read(arg, *)eps ! The perturbation
    call getarg(2, arg)
    read(arg, *)I_ini ! The initial action
    call getarg(3, arg)
    read(arg, *)N ! The maximum number of iterations
    ! Allocate the arrays
    allocate(Irms(1:N), Imean(1:N))
    ! Variables to format the datafile
    call return_base_exp(bN, eN, dble(N))
    call return_base_exp(bn_ic, en_ic, dble(n_ic))
    call return_base_exp(bI_ini, eI_ini, I_ini)
    ! Define the parameters
    params = (/ 5.0d0, -6.3d0, 6.3d0, 10.7d0, -15.8d0, 4.13d0, -9.867d0, 17.47d0, 10.1d0, -9.0d0, 15.0d0, 6.0d0, 1.83d-2, eps, -9.16d-1/)
    ! Write in the screen for the user to check if everything is correct
    write(*,*)
    write(*, "(a)")trim(path)
    write(*, "(a, f0.7)")"eps = ", eps
    write(*, "(a, i0, a, i0)")"I_ini = ", bI_ini, "e", eI_ini
    write(*, "(a, i0, a, i0)")"N = ", bN, "e", eN
    ! Define the file formats
    910 format(a, "Irms_eps=", f0.7, "_I0=", i0, "e", i0, "_nic=", i0, "e", i0, "_N=", i0, "e", i0, ".dat")
    write(unit=datafile, fmt=910)trim(path), eps, bI_ini, eI_ini, bn_ic, en_ic, bN, eN
    open(10, file=trim(datafile))
    911 format(a, "Imean_eps=", f0.7, "_I0=", i0, "e", i0, "_nic=", i0, "e", i0, "_N=", i0, "e", i0, ".dat")
    write(unit=datafile, fmt=911)trim(path), eps, bI_ini, eI_ini, bn_ic, en_ic, bN, eN
    open(11, file=trim(datafile))
    ! Generates a distribution of random numbers for theta
    call random_number(theta)
    ! Set the value for the action of all particles to I_ini
    I = I_ini
    ! Iterates the ensemble
    do j = 1, N
        ! Update I
        I = I + params(14) * dsin(2 * pi * theta)
        ! Functions
        q_factor = params(1) + params(2) * I ** 2 + params(3) * I ** 3
        E_r = params(4) * I + params(5) * dsqrt(dabs(I)) + params(6)
        v_parallel = params(7) + params(8) * dtanh(params(9) * I + params(10))
        ! Update theta
        theta = modulo(theta + params(13) * v_parallel * (params(11) / q_factor - params(12)) + params(15) * E_r / dsqrt(dabs(I)), 1.0d0)
        ! Computes the average root mean squared
        S = S + I ** 2
        Irms(j) = dsqrt((sum(S)/j)/n_ic)
        ! Computes the mean
        Imean(j) = sum(I)/n_ic
    end do
    ! Saves the data
    do j = 1, N
        write(10, "(*(f30.16))")dble(j), Irms(j)
    end do
    close(10)
    !
    write(11, "(*(f30.16))")dble(0), I_ini
    do j = 1, N
        write(11, "(*(f30.16))")dble(j), Imean(j)
    end do
    close(11)

end program I_averages