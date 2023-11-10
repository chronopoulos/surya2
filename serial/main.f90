program main

    use initialization
    use staticfields
    use loop

    implicit none

    ! Initialization
    call read_input()
    call allocate_variables()

    ! Static Fields
    call define_coordinates()
    call define_mc()
    call define_dr()
    call define_diffusivities()
    call define_matrix_coefficients() ! sets to zero

    ! Main Loop
    t=0d0
    call initialize_fields()
    do iter=1,niter

        print*, 'Timestep: ', iter, ' of ', niter
        call advance_interior_a()
        call advance_interior_b()
        t=t+dt

        if (iter/40*40 .eq. iter) then
            open(95, file='run.dat', status='unknown', access='append')
            write(95, '(6(e15.9,1x))') t, bfield(45,74), bfield(44,74), bfield(45,54), afield(120,122), afield(120,6)
            close(95)
        end if


    enddo ! iter

end program main
