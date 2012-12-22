! here are the routines contained in the main loop of the program

module loop

	use initialization

	implicit none
	save

contains

	subroutine initialize_fields()

      double precision :: harvest

      if (init_mode .eq. 'choudhuri') then
         if (myrank .eq. 0) then
            do i = 0,n+1
            do j = 0,n/np
               afield(i,j) = 0.0d0
               bfield(i,j) = dsin(theta(j))*dsin(pi*(r(i)-rmin)/(rsun-rmin))
            end do
            end do
         else if (myrank .eq. np-1) then
            do i = 0,n+1
            do j = 1,n/np+1
               afield(i,j) = 0.0d0
               bfield(i,j) = dsin(theta((np-1)*n/np+j))*dsin(pi*(r(i)-rmin)/(rsun-rmin))
            end do
            end do
         else
            do i = 0,n+1
            do j = 1,n/np
               afield(i,j) = 0.0d0
               bfield(i,j) = dsin(theta(myrank*n/np+j))*dsin(pi*(r(i)-rmin)/(rsun-rmin))
            end do
            end do
         endif
      endif

	end subroutine initialize_fields


	subroutine advance_interior_a()

   integer :: ith_first, ith_last
   double precision, allocatable  :: cp(:)
   double precision, allocatable  :: a_local(:), b_local(:), d_local(:), cp_local(:), x_local(:)
   double precision, allocatable  :: subdiag(:), diag(:), superdiag(:), rhs(:), result(:)

   !%%%%%%%%%%%%

      ! determine local extent
      if (myrank .eq. 0) then
         ith_first=0
         ith_last=n/np
      else if (myrank .eq. np-1) then
         ith_first=1
         ith_last=n/np+1
      else
         ith_first=1
         ith_last=n/np
      endif

      ! allocate
      allocate ( cp(n) )
      allocate ( a_local(ith_first:ith_last), b_local(ith_first:ith_last), cp_local(ith_first:ith_last), d_local(ith_first:ith_last), x_local(ith_first:ith_last) )
      allocate ( subdiag(n), diag(n), superdiag(n), rhs(n), result(n) )

		! first half-step
      do j=1,n

         ! every proc computes all of c-prime in a recursive loop
         cp(1)=ca(j,1)/ba(j,1)
         do i=2,n
            cp(i) = ca(j,i) / ( ba(j,i) - cp(i-1)*aa(j,i) )
         enddo

         ! partition the data into local chunks
         a_local = aa( j, myrank*n/np+1 : (myrank+1)*n/np )
         b_local = ba( j, myrank*n/np+1 : (myrank+1)*n/np )
         cp_local = cp( myrank*n/np : (myrank+1)*n/np ) ! cp_local is one longer

			! compute RHS
			do i = ith_first, ith_last
				d_local(i) = -da(j,i)*afield(j-1,i) + (1.-ea(j,i))*afield(j,i) - fa(j,i)*afield(j+1,i)
			enddo ! i

			! tridi solve
         call tridi_pscheme( n,np,myrank, a_local, b_local, cp_local, d_local, x_local )

			! store result
         afield2(j,ith_first:ith_last)=x_local

		enddo ! j

		! second half-step
		do j=ith_first,ith_last

			! compute RHS
			do i=1,n
            ! should be this, but requires ghost column communication
				   !rhs(i) = -aa(i,j)*afield2(i,j-1) + (1.-ba(i,j))*afield2(i,j) - ca(i,j)*afield2(i,j+1)
				! do this for now, until you implement ghost columns
               rhs(i) = -aa(i,j)*afield2(i,j) + (1.-ba(i,j))*afield2(i,j) - ca(i,j)*afield2(i,j)
			enddo ! i

			! tridi solve
      	subdiag = da(1:n,j)
      	diag = 1 + ea(1:n,j)
      	superdiag = fa(1:n,j)
      	call tridi_serial(subdiag,diag,superdiag,rhs,result,n)

			! store result
			afield(1:n,j)=result

		enddo ! j

	end subroutine advance_interior_a


	subroutine advance_interior_b()
   ! this is reversed relative to choudhuri eqns 33-40
   ! replaced:
   ! db <-> ab
   ! eb <-> bb
   ! fb <-> cb

   double precision, allocatable  :: subdiag(:), diag(:), superdiag(:), rhs(:), result(:)

   allocate ( subdiag(n), diag(n), superdiag(n), rhs(n), result(n) )

		! first half-step
      do j=1,n

			! compute RHS
			do i=1,n
				rhs(i) = -db(j,i)*bfield(j-1,i) + (1.-eb(j,i))*bfield(j,i) - fb(j,i)*bfield(j+1,i)
			enddo ! i

			! tridi solve
      	subdiag = ab(j,1:n)		! a-transpose
      	diag = 1 + bb(j,1:n)		! b-transpose
      	superdiag = cb(j,1:n)	! c-transpose
      	call tridi_serial(subdiag,diag,superdiag,rhs,result,n)

			! store result
			bfield2(j,1:n)=result

		enddo ! j

		! second half-step
		do j=1,n

			! compute RHS
			do i=1,n
				rhs(i) = -ab(i,j)*bfield2(i,j-1) + (1.-bb(i,j))*bfield2(i,j) - cb(i,j)*bfield2(i,j+1)
			enddo ! i

			! tridi solve
      	subdiag = db(1:n,k)
      	diag = 1 + eb(1:n,j)
      	superdiag = fb(1:n,j)
      	call tridi_serial(subdiag,diag,superdiag,rhs,result,n)

			! store result
			bfield(1:n,j)=result

		enddo ! j

	end subroutine advance_interior_b

end module loop
