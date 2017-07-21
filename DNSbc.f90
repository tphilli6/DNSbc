  module DNSbc

    implicit none

    integer,  parameter :: dp=selected_real_kind(15,307)
    real(dp), parameter :: pi=4._dp*atan(1._dp)

    ! User supplied input length scale and mesh spacing
    real(dp) :: Lx, Ly, Lz
    real(dp) :: dx, dy, dz

    integer :: My, Mz ! Mesh inflow plane dimensions
    integer :: Nx, Ny, Nz

    real(dp), allocatable, dimension(:,:,:) :: Rx, Ry, Rz, Ua
    real(dp), allocatable, dimension(:,:,:) :: bijk !Filter coefficients
  
    !if periodic, random field is also periodic
    logical :: pY=.false.
    logical :: pZ=.false.
   
    real(dp), allocatable, dimension(:,:) :: Rxbar, Rybar, Rzbar
    real(dp) :: timeOld, rold

contains

    !Initialize the DNS inflow filter
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine setupDNSFilter(LLx, LLy, LLz, ddx, ddy, ddz, MMy, MMz, ppY, ppZ)
    use mpi

      real(dp) :: LLx, LLy, LLz, ddx, ddy, ddz
      logical  :: ppY, ppZ

      integer  :: ii, jj, kk, MMy, MMz
      real(dp) :: nnx, nny, nnz
      real(dp) :: sqrt3=sqrt(3._dp)

      integer :: rank, ierr

      Lx=LLx
      Ly=LLy
      Lz=LLz
      dx=ddx
      dy=ddy
      dz=ddz
      My=MMy
      Mz=MMz
      pY=ppY
      pZ=ppZ

      ! Compute filter width
      nnx = Lx/dx
      nny = Ly/dy
      nnz = Lz/dz
      Nx = int(ceiling(2*nnx))
      Ny = int(ceiling(2*nny))
      Nz = int(ceiling(2*nnz))


      call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr)

      if (rank==0) then
        allocate( Rx(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), &
                  Ry(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), &
                  Rz(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), &
                  bijk(-Nx:Nx, -Ny:Ny, -Nz:Nz),         &
                  Ua(3, My, Mz),                        &
                  Rxbar(-Ny+1:Ny+My, -Nz+1:Nz+Mz),      &
                  Rybar(-Ny+1:Ny+My, -Nz+1:Nz+Mz),      &
                  Rzbar(-Ny+1:Ny+My, -Nz+1:Nz+Mz) )

        ! Initialize random number array
        call srand(20)

        do kk=-Nz+1,Nz+Mz
          do jj=-Ny+1,Ny+My
            do ii=-Nx,Nx
                  Rx(ii,jj,kk) = sqrt3*(2._dp*rand(0) - 1._dp)
                  Ry(ii,jj,kk) = sqrt3*(2._dp*rand(0) - 1._dp)
                  Rz(ii,jj,kk) = sqrt3*(2._dp*rand(0) - 1._dp)
            enddo
          enddo
        enddo
  
        call periodic(Rx)
        call periodic(Ry)
        call periodic(Rz)
 
        !Compute filter coefficients
        do kk=-Nz,Nz
          do jj=-Ny,Ny
            do ii=-Nx,Nx
              bijk(ii,jj,kk) = bk(ii,nnx,Nx)*bk(jj,nny,Ny)*bk(kk,nnz,Nz)
            enddo
          enddo
        enddo

        Ua=0._dp
        timeOld=0._dp
        rold=1._dp

        ! Intialize Ua
        do kk=1,Mz
          do jj=1,My
            call Ualpha(Ua(:,jj,kk), jj, kk)
          enddo
        enddo


      else

        allocate( Rx(1,1,1), Ry(1,1,1), Rz(1,1,1), bijk(1,1,1), Ua(3, My, Mz), &
                  Rxbar(1,1), Rybar(1,1), Rzbar(1,1) )
        Ua = 0._dp

      endif

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      ! computes the perturbation velocity
      call MPI_Bcast( Ua, 3*My*Mz, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)



    end subroutine setupDNSFilter

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine closeDNSFilter()

      deallocate( Rx, Ry, Rz, bijk, Ua)

    endsubroutine closeDNSFilter

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Compute the velocity perturbation
    subroutine Ualpha(Va, jj, kk)
      integer, intent(in) :: jj,kk
      real(dp), dimension(3), intent(out) :: Va

      Va(1) = sum(  bijk*Rx(:, -Ny+jj:Ny+jj, -Nz+kk:Nz+kk ) )
      Va(2) = sum(  bijk*Ry(:, -Ny+jj:Ny+jj, -Nz+kk:Nz+kk ) )
      Va(3) = sum(  bijk*Rz(:, -Ny+jj:Ny+jj, -Nz+kk:Nz+kk ) )

    end subroutine Ualpha


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Update the velocity perturbation array
    subroutine updateUalpha(time)
    use mpi

      real(dp), intent(in) :: time

      integer  :: rank, ierr
      integer  :: jj, kk
      integer  :: q
      real(dp) :: r, dt, tbar

      ! if the fractional update r<tol then consider it a full step
      real(dp) :: tol=1e-4_dp


      call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr)
      if (rank==0) then
        dt = dx ! just for clarity

        if (rold.ne.1) then
          call backstepRandomField(rold)
          timeOld = dt*( int(timeOld/dt)+1)
        end if

        q = int( (time - timeOld)/dt )
        if ( abs(q-(time-timeOld)/dt)<tol ) q=q-1
        do jj = 1,q
          call updateRandomField()
        enddo

        r = (time-timeOld)/dt - q
        if (abs(1._dp-r)<tol) then
          r=1._dp
        endif
        call updateRandomFieldPartial(r)

        timeOld=time
        rold=r

        do kk=1,Mz
          do jj=1,My
            call Ualpha(Ua(:,jj,kk), jj, kk)
          enddo
        enddo

      endif

      call MPI_Bcast( Ua, 3*My*Mz, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)


    end subroutine updateUalpha

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Update the random field (called after each time step)
    subroutine backstepRandomField(r)

      real(dp), intent(in) :: r

      integer :: ii

      ! Using the stored old random matrix solve for the previous step
      ! Rx_ii = Rx_ii*(1-r) + Rxbar*r
      Rx(Nx,:,:) = (Rx(Nx,:,:) - Rxbar*r)/(1._dp-r)
      Ry(Nx,:,:) = (Ry(Nx,:,:) - Rybar*r)/(1._dp-r)
      Rz(Nx,:,:) = (Rz(Nx,:,:) - Rzbar*r)/(1._dp-r)

      ! Work backwards resetting the random matrix
      ! Rx_ii = Rx_ii*(1-r) + Rx_ii+1*r
      do ii=Nx-1,-Nx,-1
        Rx(ii,:,:) = (Rx(ii,:,:) - Rx(ii+1,:,:)*r)/(1._dp-r)
        Ry(ii,:,:) = (Ry(ii,:,:) - Ry(ii+1,:,:)*r)/(1._dp-r)
        Rz(ii,:,:) = (Rz(ii,:,:) - Rz(ii+1,:,:)*r)/(1._dp-r)
      enddo


      ! Now that everything is back one time step
      ! complete a full shift getting to an integer time step
      do ii=-Nx,Nx-1
        Rx(ii,:,:) = Rx(ii+1,:,:)
        Ry(ii,:,:) = Ry(ii+1,:,:)
        Rz(ii,:,:) = Rz(ii+1,:,:)
      enddo

        Rx(Nx,:,:) = Rxbar
        Ry(Nx,:,:) = Rybar
        Rz(Nx,:,:) = Rzbar


    end subroutine backstepRandomField


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Update the random field (called after each time step)
    subroutine updateRandomFieldPartial(r)

      real(dp) :: r
      real(dp) :: sqrt3=sqrt(3._dp)
      integer  :: ii, jj, kk, cnt

      ! r is a fractional update. The idea is that the turbulent eddy of length Lx
      ! propagates through space at the rate of inflow velocity U. The digital filter
      ! uses atleast the grid spacing dx, dy, dz. Assuming the flow is in the
      ! x-direction, dy and dz are constant as an underlying equally spaced grid
      ! is assumed and the resulting velocity is interpolated to the actual grid.
      ! For the x-direction; however, the grid is related to time, not the physical
      ! mesh spacing. Therefore, dt<=Lx/U. In practice, dt is not constant, so to ensure
      ! that dt<=Lx/U, a smaller dt is chosen. If dt=Lx/U then the random field needs
      ! to be only updated once; however, if dt<Lx/U and dt evenly divides into Lx/U
      ! P number of times, then the random field should be updated P times. For
      ! mod(Lx/U, dt) /= 0 then a fractional update is done
      !
      !  Rnew(ii) = Rold(ii)*(1-r) + R(ii+1)*r
      !
      ! if r==1, then a full step update is done
      ! if r==0, then no update is done
      ! otherwise a linear interpolation between the two is completed
      !
      ! I only hypothsis that this will work, so now time to test it out

      !Shift random number field
      do ii=-Nx,Nx-1
        Rx(ii,:,:) = Rx(ii,:,:)*(1._dp-r) + Rx(ii+1,:,:)*r
        Ry(ii,:,:) = Ry(ii,:,:)*(1._dp-r) + Ry(ii+1,:,:)*r
        Rz(ii,:,:) = Rz(ii,:,:)*(1._dp-r) + Rz(ii+1,:,:)*r
      enddo

      !Update random numbers in last field
      cnt=0
      do kk=-Nz+1,Nz+Mz
        do jj=-Ny+1,Ny+My
          Rxbar(jj,kk) = sqrt3*(2._dp*rand(0) - 1._dp)
          Rybar(jj,kk) = sqrt3*(2._dp*rand(0) - 1._dp)
          Rzbar(jj,kk) = sqrt3*(2._dp*rand(0) - 1._dp)
          cnt=cnt+3
        enddo
      enddo

      !Update random numbers in last field
      ii=Nx
      Rx(ii,:,:) = Rx(ii,:,:)*(1._dp-r) + Rxbar*r
      Ry(ii,:,:) = Ry(ii,:,:)*(1._dp-r) + Rybar*r
      Rz(ii,:,:) = Rz(ii,:,:)*(1._dp-r) + Rzbar*r





    end subroutine updateRandomFieldPartial


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Update the random field (called after each time step)
    subroutine updateRandomField()

      real(dp) :: sqrt3=sqrt(3._dp)
      integer  :: ii, jj, kk, cnt


      !Shift random number field
      do ii=-Nx,Nx-1
        Rx(ii,:,:) = Rx(ii+1,:,:)
        Ry(ii,:,:) = Ry(ii+1,:,:)
        Rz(ii,:,:) = Rz(ii+1,:,:)
      enddo

      !Update random numbers in last field
      ii=Nx
      cnt=0
      do kk=-Nz+1,Nz+Mz
        do jj=-Ny+1,Ny+My
          Rx(ii,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
          Ry(ii,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
          Rz(ii,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
          cnt=cnt+3
        enddo
      enddo

      !update if periodic
      call periodic(Rx)
      call periodic(Ry)
      call periodic(Rz)


    end subroutine updateRandomField


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Update random field for periodic boundary conditions
    subroutine periodic(R)

    implicit none

    real(dp), dimension(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), intent(inout) :: R 

    !Periodic in y
    if (pY) then
      R(:, My+1:My+Ny,:) = R(:, 1:Ny, :)
      R(:, 1-Ny:0,:)    = R(:, My-Ny:My, :)
    endif

    if (pZ) then
      R(: ,:, Mz+1:Mz+Nz) = R(: , :, 1:Nz)
      R(: ,:, 1-Nz:0)    = R(: , :, Mz-Nz:Mz)
    endif
    
 
   


    end subroutine periodic
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Support functions for the DNS filter
    function btildek(k,n)
      integer,  intent(in) :: k
      real(dp), intent(in) :: n

      real(dp) :: btildek

      btildek = exp( -pi*0.5_dp*real(k*k,dp)/(n*n) )

    end function btildek

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    function bk(k,n,NN)
      integer, intent(in) :: k, NN
      real(dp), intent(in) :: n
      real(dp) :: bk

      real(dp) :: bksum2
      integer  :: i

      bksum2=0._dp
      do i=-NN,NN
        bksum2 = bksum2 + btildek(i, n)**2
      enddo

      bk = btildek(k,n)/sqrt( bksum2 )

    end function bk

end
