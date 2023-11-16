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
    call define_alpha()
    call define_matrix_coefficients() ! sets to zero

    ! check courant condition
    print*, 'dr/dt = ', dr/dt
    print*, 'dth*rmin/dt = ', dth*rmin/dt

    ! Main Loop
    t=0d0
    call initialize_fields()
    open(95, file=TRIM(ADJUSTL(data_dir))//"/run.dat", status='replace', action='write')
    do iter=1,niter

        call advance_interior_a()
        call advance_interior_b()
        call set_boundary_conditions()
        t=t+dt

        ! output
        if (MOD(iter, output_cadence) .eq. 0) then
            ! run.dat
            if (iter .eq. 1) then
                open(95, file=TRIM(ADJUSTL(data_dir))//"/run.dat", action='write', status='replace')
            else
                open(95, file=TRIM(ADJUSTL(data_dir))//"/run.dat", action='write', position='append')
            end if
            write(95, '(6(e15.9,1x))') t, bfield(45,74), bfield(44,74), bfield(45,54), afield(120,122), afield(120,6)
            close(95)
            ! stdout
            print*, 'Timestep: ', iter, ' of ', niter
        end if


    enddo ! iter

end program main
