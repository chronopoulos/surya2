! definitions of arrays and fields that don't change in the main loop (they're static)

module staticfields

    use initialization

    implicit none
    save

contains

    subroutine define_coordinates()

        ! grid-point resolution
        do i=0,n+1
            r(i) = rmin+i*dr
            theta(i) = thetamin+i*dth
        enddo

        ! mid-point resolution
        do i=0,2*n+2
            r_mid(i) = rmin+i*dr_mid
            theta_mid(i) = thetamin+i*dth_mid
        enddo

    end subroutine define_coordinates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Section 2.1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! meridional circulation
    subroutine define_mc()

        implicit none
        double precision :: beta_1=1.5, beta_2=1.8, epsi=2.0000001, big_gamma=3.47d10   ! MR profile parameters
        double precision :: little_gamma=0.95, m=1.5                                    ! stratification parameters
        !double precision :: psi_0=1d0                                                   ! overall scaling
        double precision :: psi_0=5.472e6                                               ! overall scaling
        double precision :: rho_top_ish=2d-7   ! this is not what we want, but it will do
        double precision :: r_p, r_0
        double precision :: dpsi_dr, dpsi_dtheta, rho
        double precision :: p, q

        ! constants
        r_p=0.61*rsun
        r_0=(rsun-rmin)/4

        ! streamfunction (mid-point resolution, avoid the boundaries)
        open(82, file=TRIM(ADJUSTL(data_dir))//"/psi.dat", action='write', status='replace')
        do i=1,2*n+1
            do j=1,2*n+1

                if(theta_mid(j) .le. pi/2.) then
                    psi_mid(i,j) = psi_0 * (r_mid(i)-r_p) * dsin( pi*(r_mid(i)-r_p)/(rsun-r_p) ) &
                        * (1.-dexp(-beta_1*theta_mid(j)**epsi)) * (1.-dexp(beta_2*(theta_mid(j)-pi/2.))) &
                        * dexp(-((r_mid(i)-r_0)/big_gamma)**2) / (r_mid(i) * dsin(theta_mid(j)))
                else
                    psi_mid(i,j) = (-1.) * psi_0 * (r_mid(i)-r_p) * dsin( pi*(r_mid(i)-r_p)/(rsun-r_p) ) &
                        * (1.-dexp(-beta_1*(pi-theta_mid(j))**epsi)) * (1.-dexp(beta_2*(pi/2. - theta_mid(j)))) &
                        * dexp(-((r_mid(i)-r_0)/big_gamma)**2) / (r_mid(i) * dsin(theta_mid(j)))
                end if

                if (r_mid(i) .le. r_p) then
                    psi_mid(i,j) = 0.
                end if

                ! output
                write(82, '(3(e15.9,1x))') theta_mid(j), r_mid(i), psi_mid(i,j)

            enddo ! j
        enddo ! i
        close(82)
      
        ! velocity field as curl of the streamfunction
        do i=1,2*n+1
            do j=1,2*n+1

                dpsi_dr     = ( psi_mid(i+1,j)-psi_mid(i-1,j) )/(2.*dr_mid)
                dpsi_dtheta = ( psi_mid(i,j+1)-psi_mid(i,j-1) )/(2.*dth_mid)

                ! stratification
                rho = rho_top_ish * (rsun/r_mid(i) - little_gamma)**m

                v_r_mid(i,j)     =  ( dpsi_dtheta + psi_mid(i,j)/dtan(theta_mid(j)) ) / (rho*r_mid(i))
                v_theta_mid(i,j) = -( dpsi_dr + psi_mid(i,j)/r_mid(i) )/rho

            enddo ! j
        enddo ! i

    end subroutine define_mc


    ! differential rotation
    subroutine define_dr()

        implicit none
        double precision :: r_t_frac = 0.7, d_t_frac = 0.025
        double precision :: omega_rz = 2.719d-6, omega_eq = 2.895d-6
        double precision :: alpha_2 = -3.939e-7, alpha_4 = -4.218e-7
        double precision :: omega_scz, r_t, d_t, p, q

        r_t = r_t_frac * rsun
        d_t = d_t_frac * rsun
 
        open(25, file=TRIM(ADJUSTL(data_dir))//"/diffrot.dat", action='write', status='replace')
        do i=0,n+1
            do j=0,n+1
                omega_scz = omega_eq + alpha_2*dcos(theta(j))**2 + alpha_4*dcos(theta(j))**4
                omega(i,j) = omega_rz + 0.5*( 1. + erf( (r(i)-r_t)/d_t ) ) * (omega_scz - omega_rz)
                write(25, '(3(e15.9,1x))') theta(j), r(i), omega(i,j)
            enddo ! j
        enddo ! i

    end subroutine define_dr


    ! Diffusivities
    subroutine define_diffusivities()

        implicit none
        double precision :: eta_rz=2.2d8, eta_scz=2.6d12, eta_scz1=4d10 ! cgs
        double precision :: r_bcz, rprime_bcz
        double precision :: r_t_frac = 0.7, d_t_frac = 0.025
        double precision :: r_t, d_t

        r_t=r_t_frac*rsun
        d_t=d_t_frac*rsun

        r_bcz = 0.7*rsun
        rprime_bcz = 0.72*rsun

        do i=0,n+1
            do j=0,n+1
                eta_a(i,j) = eta_rz + 0.5*eta_scz1*( 1. + erf( (r(i)-rprime_bcz)/d_t ) )
                eta_b(i,j) = eta_rz + 0.5*eta_scz*( 1. + erf( (r(i)-r_bcz)/d_t ) )
            enddo ! j
        enddo ! i

    end subroutine define_diffusivities

    ! Diffusivities
    subroutine define_alpha()

        implicit none
        double precision :: alpha_0 = 2500.0

        do i=0,n+1
            do j=0,n+1
                alpha(i,j) = 0.25 * alpha_0 * dcos(theta(j)) * (1. + erf((r(i)-0.95*rsun)/(0.025*rsun))) * (1. - erf((r(i) - rsun) / (0.025*rsun)))
            enddo ! j
        enddo ! i

    end subroutine define_alpha

    ! Matrix Coefficients
    subroutine define_matrix_coefficients()

        implicit none
        integer :: imid, jmid
        double precision :: ur_ij, ur_iminus, ur_iplus, utheta_ij, utheta_jminus, utheta_jplus
        double precision :: sinplus, sinminus
        double precision :: urbar_ij, urbar_iminus, urbar_iplus

        aa=0d0
        ba=0d0
        ca=0d0
        da=0d0
        ea=0d0
        fa=0d0

        do i=1,n
            imid=2*i
            do j=1,n
                jmid=2*i
                ! repeated definitions

                ur_ij = v_r_mid(imid,jmid)/r_mid(imid)
                ur_iminus = v_r_mid(imid-1,jmid)/r_mid(imid-1)
                ur_iplus = v_r_mid(imid+1,jmid)/r_mid(imid+1)

                utheta_ij = v_theta_mid(imid,jmid)/(r(i)*dsin(theta_mid(jmid)))
                utheta_jminus = v_theta_mid(imid,jmid-1)/(r(i)*dsin(theta_mid(jmid-1)))
                utheta_jplus = v_theta_mid(imid,jmid+1)/(r(i)*dsin(theta_mid(jmid+1)))

                sinplus=dsin(theta(j)+dt/2.)
                sinminus=dsin(theta(j)-dt/2.)

                urbar_ij = ( ur_ij - (eta_b(i+1,j)-eta_b(i-1,j))/(2.*dr) )/r_mid(imid)
                ! these next two are kind of a hack to avoid defining eta with mid-point precision
                !  -> d(eta)/d(r) is defined as the first order forward difference
                urbar_iminus = ( ur_iminus - (eta_b(i,j)-eta_b(i-1,j))/dr )/r_mid(imid-1)
                urbar_iplus = ( ur_iplus - (eta_b(i+1,j)-eta_b(i,j))/dr )/r_mid(imid+1)

                ! matrix coefficients

                aa(i,j) = -( eta_a(i,j)/(r(i)*dth)**2 + utheta_ij/dth*sinminus &
                    * ( 0.5 + dt/(4.*dth)*v_theta_mid(imid,jmid-1)/r(i) ) ) * dt/2.

                ba(i,j) = ( eta_a(i,j)*( 2./(r(i)*dth)**2 + 1./(dtan(theta(j))*r(i)**2*dth) &
                    + 1./(2.*r(i)**2*dsin(theta(j))**2) ) + utheta_ij/dth*( &
                    sinplus *(0.5 + dt/(4.*dth)*utheta_jplus *dsin(theta(j)))  &
                    - sinminus*(0.5 - dt/(4.*dth)*utheta_jminus*dsin(theta(j))))) * dt/2.

                ca(i,j) = -( eta_a(i,j)*( 1./(r(i)*dth)**2  + 1./(dtan(theta(j))*r(i)**2*dth) ) &
                    - ur_ij*sinplus/dth * (0.5 - dt/(4.*dr)*v_theta_mid(imid,jmid+1)/(r(i))) ) * dt/2.

                da(i,j) = -( eta_a(i,j)/dr**2 + ur_ij*r_mid(imid-1)/dr &
                    * (0.5 + dt/(4.*dr)*v_r_mid(imid-1,j)) ) * dt/2.

                ea(i,j) = (eta_a(i,j) * (2./dr**2 + 2./(r(i)*dr) + 1./(2.*(r(i)*dsin(theta(j)))**2)) &
                    + ur_ij/dr*( dr/2. + dt*r(i)/(4.*dr)*(v_r_mid(imid+1,j) + v_r_mid(imid-1,j)) ) ) * dt/2.

                fa(i,j) = -( eta_a(i,j)*(1./dr**2 + 2./(r(i)*dr)) - ur_ij*r_mid(imid+1)/dr &
                    * (0.5 - dt/(4.*dr)*v_r_mid(imid+1,j)) ) * dt/2.

                ! this is reversed relative to choudhuri eqns 33-40
                ! replaced:
                !   db <-> ab
                !   eb <-> bb
                !   fb <-> cb

                db(i,j) = -( eta_b(i,j)/(r(i)*dth)**2 + v_theta_mid(imid,jmid-1)/(r(i)*dth) &
                    * (0.5 + dt/(4.*dth)*v_theta_mid(imid,jmid-1)/r(i)) ) * dt/2.

                eb(i,j) = ( eta_b(i,j)*( 2./(r(i)*dth)**2 &
                    + 1./(dtan(theta(j))*r(i)**2*dth) + 0.5/(r(i)*dsin(theta(j)))**2) &
                    + v_theta_mid(imid,jmid+1)/(r(i)*dth)*(0.5+dt/(4.*dth)*utheta_ij*dsin(theta(j))) &
                    - v_theta_mid(imid,jmid-1)/(r(i)*dth)*(0.5-dt/(4.*dth)*utheta_ij*dsin(theta(j)))) * dt/2.

                fb(i,j) = -( eta_b(i,j)*( 1./(r(i)*dth)**2 + 1./(dtan(theta(j))*r(i)**2*dth) ) &
                - v_theta_mid(imid,jmid+1)/(r(i)*dth) &
                * ( 0.5 - dt/(4.*dr)*utheta_ij*dsin(theta(j)))) * dt/2.

                ab(i,j) = -( eta_b(i,j)/dr**2 + urbar_ij/dr*r_mid(imid-1)*( 0.5 + dt/(4.*dr) &
                    * urbar_iminus*r_mid(imid-1))) * dt/2.

                bb(i,j) = ( eta_b(i,j)*( 2./dr**2 + 2./(r(i)*dr) + 0.5/(r(i)*dsin(theta(j)))**2 ) &
                    + urbar_ij + urbar_ij/dr*( dr/2. + r(i)*dt/(4.*dr)*(urbar_iplus*r_mid(imid+1) &
                    + urbar_iminus*r_mid(imid-1)))) * dt/2.

                cb(i,j) = eta_b(i,j)*( 1./dr**2 + 1./(r(i)*dr) ) - urbar_ij/dr*r_mid(imid+1) &
                    * (0.5 - dt/(4.*dr)*urbar_iplus*r_mid(imid+1))

            enddo ! j
        enddo ! i

   end subroutine define_matrix_coefficients

end module staticfields
