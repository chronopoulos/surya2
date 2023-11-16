! here are the routines contained in the main loop of the program

module loop

    use initialization

    implicit none
    save

contains

    subroutine initialize_fields()

        double precision :: harvest

        call random_seed()

        if (init_mode .eq. 'choudhuri') then
            do i = 0,n+1
                do j = 0,n+1
                    afield(i,j) = 0.0d0
                    bfield(i,j) = dsin(theta(j))*dsin(pi*(r(i)-rmin)/(rsun-rmin))
                end do
            end do
        else if (init_mode .eq. 'random') then
            do i=1,n
                do j=1,n
                    call random_number(harvest)
                    afield(i,j)=(harvest-0.5)*10.
                    call random_number(harvest)
                    bfield(i,j)=(harvest-0.5)*10.
                enddo ! j
            enddo ! i
        endif

    end subroutine initialize_fields


    subroutine advance_interior_a()

        ! first half-step
        do j=1,n

            ! compute RHS
            do i=1,n
                rhs(i) = -da(j,i)*afield(j-1,i) + (1.-ea(j,i))*afield(j,i) - fa(j,i)*afield(j+1,i) + alpha(j,i)*bfield(j,i)*dt/2.
            enddo ! i

            ! tridi solve
            subdiag = aa(j,1:n)     ! a-transpose
            diag = 1 + ba(j,1:n)    ! b-transpose
            superdiag = ca(j,1:n)   ! c-transpose
            call solve_tridiag(subdiag,diag,superdiag,rhs,result,n)

            ! store result
            afield2(j,1:n)=result

        enddo ! j

        ! second half-step
        do j=1,n

            ! compute RHS
            do i=1,n
                rhs(i) = -aa(i,j)*afield2(i,j-1) + (1.-ba(i,j))*afield2(i,j) - ca(i,j)*afield2(i,j+1) + alpha(i,j)*bfield(i,j)*dt/2.
            enddo ! i

            ! tridi solve
            subdiag = da(1:n,j)
            diag = 1 + ea(1:n,j)
            superdiag = fa(1:n,j)
            call solve_tridiag(subdiag,diag,superdiag,rhs,result,n)

            ! store result
            afield(1:n,j)=result

        enddo ! j

    end subroutine advance_interior_a

    subroutine set_boundary_conditions()
        ! bottom
        do j=0,n+1
            afield(0,j)=0.0d0
            bfield(0,j)=0.0d0
        end do
        ! top
            ! tmp
            !do nup=1,11
            !    do l=1,96
            !        ss1(l)=ss1_p(l)*ss(n,dq)
            !    end do
            !    do j=0,n+1
            !        call coeff(n,ss1,dq,j,cf)
            !        afield(n+1,j) = afield(n,j) + cf*dp
            !    end do
            !end do
        do j=0,n+1
            bfield(n+1,j)= 0.
            afield(n+1,j) = 0.  ! mp
        end do
    end subroutine set_boundary_conditions


    subroutine advance_interior_b()
        ! this is reversed relative to choudhuri eqns 33-40
        ! replaced:
        ! db <-> ab
        ! eb <-> bb
        ! fb <-> cb

        double precision :: br, bt

        ! first half-step
        do j=1,n

            ! compute RHS
            do i=1,n
                br = (afield(j,i+1)*dsin(theta(i)+dth) - afield(j,i-1)*dsin(theta(i)-dth)) * dt/(4*dth)
                bt = (afield(j-1,i)*(r(j)-dr) - afield(j+1,i)*(r(j)+dr)) * dsin(theta(i)) * dt/(4*r(j)*dr)
                rhs(i) = -db(j,i)*bfield(j-1,i) + (1.-eb(j,i))*bfield(j,i) - fb(j,i)*bfield(j+1,i) + br*dror(j,i) + bt*drot(j,i)
            enddo ! i

            ! tridi solve
            subdiag = ab(j,1:n)     ! a-transpose
            diag = 1 + bb(j,1:n)    ! b-transpose
            superdiag = cb(j,1:n)   ! c-transpose
            call solve_tridiag(subdiag,diag,superdiag,rhs,result,n)

            ! store result
            bfield2(j,1:n)=result

        enddo ! j

        ! second half-step
        do j=1,n

            ! compute RHS
            do i=1,n
                br = (afield(i,j+1)*dsin(theta(j)+dth) - afield(j,j-1)*dsin(theta(j)-dth)) * dt/(4*dth)
                bt = (afield(i-1,j)*(r(i)-dr) - afield(i+1,j)*(r(i)+dr)) * dsin(theta(j)) * dt/(4*r(i)*dr)
                rhs(i) = -ab(i,j)*bfield2(i,j-1) + (1.-bb(i,j))*bfield2(i,j) - cb(i,j)*bfield2(i,j+1) + br*dror(i,j) + bt*drot(i,j)
            enddo ! i

            ! tridi solve
            subdiag = db(1:n,j)
            diag = 1 + eb(1:n,j)
            superdiag = fb(1:n,j)
            call solve_tridiag(subdiag,diag,superdiag,rhs,result,n)

            ! store result
            bfield(1:n,j)=result

        enddo ! j

    end subroutine advance_interior_b

end module loop
