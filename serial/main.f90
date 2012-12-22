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
	do iter=1,niter
      print*, 'Timestep: ', iter
		call advance_interior_a()
      call advance_interior_b()
		t=t+dt
	enddo ! iter

end program main
