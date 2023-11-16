! read input, define parameters, allocate variables, and define standalone functions
! each of the other modules imports this module, as does the main program

module initialization

    implicit none
    save

        integer :: i, j, irow, icol, ir, ith, iter

        ! coordinate stuff
        double precision :: rsun, pi, rmin_frac, rmax_frac, thetamin_frac, thetamax_frac
        double precision :: rmin, rmax, thetamin, thetamax
        integer :: n
        double precision :: dr, dr_mid, dth, dth_mid

        ! time stuff
        double precision :: t, dt
        integer :: niter

        ! arrays
        double precision, allocatable :: r(:), r_mid(:), theta(:), theta_mid(:)
        double precision, allocatable :: afield(:,:), afield2(:,:)
        double precision, allocatable :: bfield(:,:), bfield2(:,:)
        double precision, allocatable :: alpha(:,:)

        ! section 2.1 stuff
        double precision, allocatable :: psi_mid(:,:), v_r_mid(:,:), v_theta_mid(:,:)
        double precision, allocatable :: omega(:,:), dror(:,:), drot(:,:)
        double precision, allocatable :: eta_a(:,:), eta_b(:,:)

        ! matrix coefficients
        double precision, allocatable :: aa(:,:), ba(:,:), ca(:,:), da(:,:), ea(:,:), fa(:,:)
        double precision, allocatable :: ab(:,:), bb(:,:), cb(:,:), db(:,:), eb(:,:), fb(:,:)

        ! tridi stuff
        double precision, allocatable :: subdiag(:), diag(:), superdiag(:), rhs(:), result(:)

        ! initialization
        character(len=32) :: init_mode

        ! output
        character(len=512) :: data_dir
        integer :: output_cadence

contains

    subroutine read_input()
        ! read in the input parameters from inputs.txt

        namelist /domain/ rsun, pi, rmin_frac, rmax_frac, thetamin_frac, thetamax_frac, n
        namelist /timestepping/ dt, iter, niter
        namelist /init/ init_mode
        namelist /output/ data_dir, output_cadence

        open(99, file='input.txt', status='old')

        read(99, nml=domain)

        rmin=rmin_frac*rsun
        rmax=rmax_frac*rsun
        dr=(rmax-rmin)/(n+1)
        dr_mid=(rmax-rmin)/(2*n+2)
        thetamin=thetamin_frac*pi
        thetamax=thetamax_frac*pi
        dth=(thetamax-thetamin)/(n+1)
        dth_mid=(thetamax-thetamin)/(2*n+2)

        read(99, nml=timestepping)
        read(99, nml=init)
        read(99, nml=output)

    end subroutine read_input


    subroutine allocate_variables()

        ! coordinates, dynamic variables
        allocate( r(0:n+1), theta(0:n+1) )
        allocate( r_mid(0:2*n+2), theta_mid(0:2*n+2) )
        allocate( afield(0:n+1,0:n+1), afield2(0:n+1,0:n+1) )
        allocate( bfield(0:n+1,0:n+1), bfield2(0:n+1,0:n+1) )
        allocate( alpha(0:n+1,0:n+1) )

        ! tridi stuff
        allocate( subdiag(n), diag(n), superdiag(n), rhs(n), result(n) )

        ! section 2.1
        allocate( psi_mid(0:2*n+2,0:2*n+2), v_r_mid(0:2*n+2,0:2*n+2), v_theta_mid(0:2*n+2,0:2*n+2) )
        allocate( omega(0:n+1,0:n+1), dror(0:n+1,0:n+1), drot(0:n+1,0:n+1) )
        allocate( eta_a(0:n+1,0:n+1), eta_b(0:n+1,0:n+1) )

        ! matrix coefficients
        allocate( aa(0:n+1,0:n+1), ba(0:n+1,0:n+1), ca(0:n+1,0:n+1), da(0:n+1,0:n+1), ea(0:n+1,0:n+1), fa(0:n+1,0:n+1) )
        allocate( ab(0:n+1,0:n+1), bb(0:n+1,0:n+1), cb(0:n+1,0:n+1), db(0:n+1,0:n+1), eb(0:n+1,0:n+1), fb(0:n+1,0:n+1) )

    end subroutine allocate_variables


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Standalone Functions
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine solve_tridiag(a,b,c,v,x,n)
        ! a - sub-diagonal
        ! b - the main diagonal
        ! c - super-diagonal
        ! v - right-hand side
        ! x - result
        ! n - number of equations

        implicit none

        integer,intent(in) :: n
        double precision,dimension(n),intent(in) :: a,b,c,v
        double precision,dimension(n),intent(out) :: x
        double precision,dimension(n) :: bp,vp
        double precision :: m
        integer i

        ! Make copies of the b and v variables so that they are unaltered by this sub

        bp(:) = b(:)
        vp(:) = v(:)

        ! The first pass (setting coefficients)

        do i = 2,n
            m = a(i)/bp(i-1)
            bp(i) = b(i) - m*c(i-1)
            vp(i) = v(i) - m*vp(i-1)
        end do

        x(n) = vp(n)/bp(n)

        ! The second pass (back-substition)

        do i = n-1, 1, -1
            x(i) = (vp(i) - c(i)*x(i+1))/bp(i)
        end do

    end subroutine solve_tridiag


end module initialization
