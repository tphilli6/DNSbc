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
                  Ua(3, My, Mz) )
   
         
        ! Initialize random number array
        call srand(20) !sets the random seed for reproduction of results
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

      else

        allocate( Rx(1,1,1), Ry(1,1,1), Rz(1,1,1), bijk(1,1,1), Ua(3, My, Mz) )
        Ua = 0._dp

      endif

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      ! computes the perturbation velocity
      call updateUalpha()
 


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
    subroutine updateUalpha()
    use mpi

      integer :: rank, ierr
      integer :: ii, jj, kk

      call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr)
      if (rank==0) then

        call updateRandomField()
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
    subroutine updateRandomField()

      real(dp) :: sqrt3=sqrt(3._dp)
      integer  :: ii, jj, kk
      
      !Shift random number field     
      do ii=-Nx,Nx-1
        Rx(ii,:,:) = Rx(ii+1,:,:)
        Ry(ii,:,:) = Ry(ii+1,:,:)
        Rz(ii,:,:) = Rz(ii+1,:,:)
      enddo

      !Update random numbers in last field
      ii=Nx
      do kk=-Nz+1,Nz+Mz
        do jj=-Ny+1,Ny+My
          Rx(ii,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
          Ry(ii,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
          Rz(ii,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
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
