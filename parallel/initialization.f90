! read input, declare and allocate variables, define parameters and standalone functions
! each of the other modules imports this module, as does the main program

module initialization

      use mpi
   	implicit none
   	save  ! might need this to "preserve data values"

      integer :: i, j, k, irow, icol, ir, ith, iter

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

      ! section 2.1 stuff
      double precision, allocatable :: psi_mid(:,:), v_r_mid(:,:), v_theta_mid(:,:)
      double precision, allocatable :: omega(:,:)
      double precision, allocatable :: eta_a(:,:), eta_b(:,:)

      ! matrix coefficients
      double precision, allocatable :: aa(:,:), ba(:,:), ca(:,:), da(:,:), ea(:,:), fa(:,:)
      double precision, allocatable :: ab(:,:), bb(:,:), cb(:,:), db(:,:), eb(:,:), fb(:,:)

      ! Field initialization
      character(len=32) :: init_mode

      ! Misc input
      logical :: verbose

      ! MPI stuff
      integer :: myrank, np, ierr, stat(MPI_STATUS_SIZE)

      ! Timing stuff
      double precision :: tstart, tstop


contains

   subroutine initialize_mpi()

      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
         if (mod(n,np) .ne. 0) stop 'ERROR: np must evenly divide n'
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

   end subroutine initialize_mpi


   subroutine finalize_mpi()

      call MPI_Finalize(ierr)

   end subroutine finalize_mpi


   subroutine read_input()
   ! read in the input parameters from inputs.txt

      namelist /domain/ rsun, pi, rmin_frac, rmax_frac, thetamin_frac, thetamax_frac, n
      namelist /timestepping/ dt, iter, niter
      namelist /init/ init_mode
      namelist /misc/ verbose

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

      read(99, nml=misc)

   end subroutine read_input


	subroutine allocate_variables()

      integer :: ith_first, ith_last

      ! coordinates
      allocate( r(0:n+1), theta(0:n+1) )
      allocate( r_mid(0:2*n+2), theta_mid(0:2*n+2) )

      ! section 2.1
      allocate( psi_mid(0:2*n+2,0:2*n+2), v_r_mid(0:2*n+2,0:2*n+2), v_theta_mid(0:2*n+2,0:2*n+2) )
      allocate( omega(0:n+1,0:n+1) )
      allocate( eta_a(0:n+1,0:n+1), eta_b(0:n+1,0:n+1) )

      ! matrix coefficients
      allocate( aa(0:n+1,0:n+1), ba(0:n+1,0:n+1), ca(0:n+1,0:n+1), da(0:n+1,0:n+1), ea(0:n+1,0:n+1), fa(0:n+1,0:n+1) )
      allocate( ab(0:n+1,0:n+1), bb(0:n+1,0:n+1), cb(0:n+1,0:n+1), db(0:n+1,0:n+1), eb(0:n+1,0:n+1), fb(0:n+1,0:n+1) )

      ! dynamic variables: partitioned over theta
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
         allocate( afield(0:n+1,ith_first:ith_last), afield2(0:n+1,ith_first:ith_last) )
         allocate( bfield(0:n+1,ith_first:ith_last), bfield2(0:n+1,ith_first:ith_last) )


	end subroutine allocate_variables


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Standalone Functions
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   subroutine tridi_serial(a,b,c,v,x,n)
      ! a - sub-diagonal
      ! b - the main diagonal
      ! c - super-diagonal
      ! v - right-hand side
      ! x - result
      ! n - number of equations

      implicit none

      integer,intent(in) :: n
      real(8),dimension(n),intent(in) :: a,b,c,v
      real(8),dimension(n),intent(out) :: x
      real(8),dimension(n) :: bp,vp
      real(8) :: m
      integer i

      ! Make copies of the b and v variables so that they are unaltered by this sub

         bp(1) = b(1)
         vp(1) = v(1)

      !The first pass (setting coefficients)

         do i = 2,n
            m = a(i)/bp(i-1)
            bp(i) = b(i) - m*c(i-1)
            vp(i) = v(i) - m*vp(i-1)
         end do

         x(n) = vp(n)/bp(n)

      !The second pass (back-substition)

         do i = n-1, 1, -1
            x(i) = (vp(i) - c(i)*x(i+1))/bp(i)
         end do

   end subroutine tridi_serial


   subroutine tridi_pscheme( n,np,myrank, a_local, b_local, cp_local, d_local, x_local )
   ! n - problem size
   ! np - total number of procs
   ! myrank - my rank..
   ! a(n),b(n),c(n) - sub, super, and main diagonals (memory inefficient, but whatever)
   ! x(n/np) - local x (solution)
   ! d(n/np) - local

      implicit none

      ! inputs, outputs
      integer, intent(in) :: n, np, myrank
      double precision, intent(in)  :: a_local(n/np), b_local(n/np), d_local(n/np)
      double precision, intent(in)  :: cp_local(0:n/np) ! has length n/np+1, and cp_local(0) is ok
      double precision, intent(out) :: x_local(n/np)

      ! local variables
      integer :: ii, jj, kk
      double precision :: ss_local(n/np), tt_local(n/np), dp_local(n/np)
      double precision :: ghostpoint, prod1, prod2, summ

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! define ss_local and tt_local for the forward sweep
      do ii=1,n/np
         ss_local(ii) = a_local(ii) / ( cp_local(ii-1)*a_local(ii) - b_local(ii) )
         tt_local(ii) = d_local(ii) / (-cp_local(ii-1)*a_local(ii) + b_local(ii) )
      enddo

      ! forward sweep: d-prime
      if (myrank .eq. 0) then

         dp_local(1)=d_local(1)/b_local(1)

         prod1=1d0
         do kk=2,n/np
            prod1 = prod1*ss_local(kk)
         enddo

         summ=0d0
         do jj=2,n/np
            prod2=1d0
            do kk=jj+1,n/np
               prod2 = prod2*ss_local(kk)
            enddo
            summ = summ + tt_local(jj)*prod2
         enddo

         dp_local(n/np) = dp_local(1)*prod1 + summ
         ! this is eqn applepi in december 17, 2012 notes

         ! send it off
         if (np .gt. 1) then
            call MPI_Send( dp_local(n/np), 1, MPI_Double_Precision, myrank+1, 0, MPI_Comm_World, ierr )
         endif

         ! then go back and fill in the rest (recursively)
         do ii=2,n/np-1
            dp_local(ii) = ss_local(ii)*dp_local(ii-1) + tt_local(ii)
         enddo

      else ! myrank .ne. 0

         ! first receive the ghostpoint
         call MPI_Recv( ghostpoint, 1, MPI_Double_Precision, myrank-1, 0, MPI_Comm_World, stat, ierr )
         
         prod1=1d0
         do kk=1,n/np
            prod1 = prod1*ss_local(kk)
         enddo

         summ=0d0
         do jj=1,n/np
            prod2=1d0
            do kk=jj+1,n/np
               prod2 = prod2*ss_local(kk)
            enddo
            summ = summ + tt_local(jj)*prod2
         enddo

         dp_local(n/np) = ghostpoint*prod1 + summ
         ! this is eqn applepi in december 17, 2012 notes

         ! send it off
         if (myrank .ne. np-1) then
            call MPI_Send( dp_local(n/np), 1, MPI_Double_Precision, myrank+1, 0, MPI_Comm_World, ierr )
         endif

         ! then go back and fill in the rest (recursively)
         dp_local(1) = ss_local(1)*ghostpoint + tt_local(1)
         do ii=2,n/np-1
            dp_local(ii) = ss_local(ii)*dp_local(ii-1) + tt_local(ii)
         enddo

      endif ! end of forward sweep: dprime

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! define ss_local and tt_local for the back substitution
      do ii=1,n/np
         ss_local(ii) = -cp_local(ii)
         tt_local(ii) =  dp_local(ii)
      enddo

      ! back-substitution: x
      if (myrank .eq. np-1) then

         x_local(n/np)=dp_local(n/np)

         prod1=1d0
         do kk = (n/np-1), 1, -1
            prod1 = prod1*ss_local(kk)
         enddo

         summ=0d0
         do jj=(n/np-1), 1, -1
            prod2=1d0
            do kk = (jj-1), 1, -1
               prod2 = prod2*ss_local(kk)
            enddo
            summ = summ + tt_local(jj)*prod2
         enddo

         x_local(1) = x_local(n/np)*prod1 + summ
         ! this is eqn pumpkinpi in december 17, 2012 notes

         ! send it off
         if (np .gt. 1) then
            call MPI_Send( x_local(1), 1, MPI_Double_Precision, myrank-1, 0, MPI_Comm_World, ierr )
         endif

         ! then go back and fill in the rest (recursively)
         do ii = (n/np-1), 2, -1
            x_local(ii) = ss_local(ii)*x_local(ii+1) + tt_local(ii)
         enddo

      else ! myrank .ne. np-1

         ! first receive the ghostpoint
         call MPI_Recv( ghostpoint, 1, MPI_Double_Precision, myrank+1, 0, MPI_Comm_World, stat, ierr )

         prod1=1d0
         do kk = (n/np), 1, -1
            prod1 = prod1*ss_local(kk)                            ! ss should be ss_local! and the others
         enddo

         summ=0d0
         do jj=(n/np), 1, -1
            prod2=1d0
            do kk = (jj-1), 1, -1
               prod2 = prod2*ss_local(kk)
            enddo
            summ = summ + tt_local(jj)*prod2
         enddo

         x_local(1) = ghostpoint*prod1 + summ
         ! this is eqn pumpkinpi in december 17, 2012 notes

         ! send it off
         if (myrank .ne. 0) then
            call MPI_Send( x_local(1), 1, MPI_Double_Precision, myrank-1, 0, MPI_Comm_World, ierr )
         endif

         ! then go back and fill in the rest (recursively)
         x_local(n/np) = ss_local(n/np)*ghostpoint + tt_local(n/np)
         do ii = (n/np-1), 2, -1
            x_local(ii) = ss_local(ii)*x_local(ii+1) + tt_local(ii)
         enddo

      endif ! end of back substitution: x

   end subroutine tridi_pscheme


end module initialization
