program main

   use initialization
   use staticfields
   use loop

   implicit none

!~~~~~~~~~~~~~~~~
!~~~~~~~~
!~~~~
!~~
!~
!


!%%%%%%%%%%%%%%%%%%%%%%
! Initialization
!%%%%%%%%%%%%%%%%%%%%%%

   call initialize_mpi()
   call read_input()
	call allocate_variables()


!%%%%%%%%%%%%%%%%%%%%%%
! Static Fields
!%%%%%%%%%%%%%%%%%%%%%%

	call define_coordinates()
   call define_mc()
   call define_dr()
   call define_diffusivities()
   call define_matrix_coefficients() ! sets to zero


!%%%%%%%%%%%%%%%%%%%%%%
! Main Loop
!%%%%%%%%%%%%%%%%%%%%%%

	t=0d0
   call initialize_fields()

   call MPI_Barrier(MPI_Comm_World, ierr)
   tstart=MPI_Wtime()

	do iter=1,niter
      if ((myrank .eq. 0) .and. (verbose)) then
         print*, 'Timestep: ', iter
      endif
		call advance_interior_a()
      !call advance_interior_b()
		t=t+dt
	enddo ! iter

   call MPI_Barrier(MPI_Comm_World, ierr)
   tstop=MPI_Wtime()

!%%%%%%%%%%%%%%%%%%%%%%
! Print Elapsed Wall Time
!%%%%%%%%%%%%%%%%%%%%%%

   if (myrank .eq. 0) then
      print*, 'Walltime: ', tstop-tstart
   endif

!%%%%%%%%%%%%%%%%%%%%%%
! Finalization
!%%%%%%%%%%%%%%%%%%%%%%

   call finalize_mpi()

end program main
