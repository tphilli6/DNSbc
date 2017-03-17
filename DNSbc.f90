  module DNSbc
    integer,  parameter :: dp=selected_real_kind(15,307)
    real(dp), parameter :: pi=4._dp*atan(1._dp)

    ! User supplied input length scale and mesh spacing
    real(dp) :: Lx=0.125           
    real(dp) :: Ly=0.125
    real(dp) :: Lz=0.125

    real(dp) :: dx=0.06125
    real(dp) :: dy=0.06125
    real(dp) :: dz=0.06125

    integer :: My=10, Mz=10 ! Mesh inflow plane dimensions
    integer :: Nx, Ny, Nz

    real(dp), allocatable, dimension(:,:,:) :: Rx, Ry, Rz
    real(dp), allocatable, dimension(:,:,:) :: bijk !Filter coefficients
   
    real(dp), dimension(3) :: ubar = (/ 2.0_dp, 0._dp, 0._dp /)
    real(dp) :: uu=0.1_dp
    real(dp) :: vv=0.03_dp
    real(dp) :: ww=0.01_dp
    real(dp) :: uv=0.02_dp
    real(dp) :: uw=0.01_dp
    real(dp) :: vw=0.015_dp

contains

    !Initialize the DNS inflow filter
    subroutine setupDNSFilter(LLx, LLy, LLz, ddx, ddy, ddz, My, Mz)

      real(dp) :: LLx, LLy, LLz, ddx, ddy, ddz
      integer  :: ii, jj, kk
      real(dp) :: nnx, nny, nnz
      real(dp) :: sqrt3=sqrt(3._dp)

      write(*,'(A)') 'Setting up DNS Boundary Conditions...'
      Lx=LLx
      Ly=LLy
      Lz=LLz
      dx=ddx
      dy=ddy
      dz=ddz
      ! Compute filter width 
      nnx = Lx/dx
      nny = Ly/dy
      nnz = Lz/dz
      Nx = int(ceiling(2*nnx))
      Ny = int(ceiling(2*nny))
      Nz = int(ceiling(2*nnz))

      allocate( Rx(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz),&
                Ry(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz),&
                Rz(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz),&
                bijk(-Nx:Nx, -Ny:Ny, -Nz:Nz) )

      ! Initialize random number array
      do kk=-Nz+1,Nz+Mz
        do jj=-Ny+1,Ny+My
          do ii=-Nx,Nx
                Rx(ii,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
                Ry(ii,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
                Rz(ii,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
          enddo
        enddo
      enddo

      !Compute filter coefficients
      do kk=-Nz,Nz
        do jj=-Ny,Ny
          do ii=-Nx,Nx
            bijk(ii,jj,kk) = bk(ii,nnx,Nx)*bk(jj,nny,Ny)*bk(kk,nnz,Nz)
          enddo
        enddo
      enddo


    end subroutine setupDNSFilter

    subroutine closeDNSFilter()
            
      deallocate( Rx, Ry, Rz, bijk)

    endsubroutine closeDNSFilter
   
    !Primary routine to compute inflow velocity.
    subroutine ComputeDNSVelocity(vel, jj, kk)

      integer, intent(in) :: jj, kk
      real(dp), dimension(3), intent(out) :: vel

      real(dp) :: a11, a21, a22, a31, a32, a33
      real(dp), dimension(3) :: Ua

      a11 = sqrt(uu)
      a21 = uv/a11
      a22 = sqrt( vv - a21**2 )
      a31 = uw/a11
      a32 = ( vw - a21*a31 )/a22
      a33 = sqrt( ww - a31**2 - a32**2 )

      call Ualpha(Ua, jj, kk)

      vel(1) = ubar(1) + Ua(1)*a11
      vel(2) = ubar(2) + Ua(1)*a21 + Ua(2)*a22
      vel(3) = ubar(3) + Ua(1)*a31 + Ua(2)*a32 +Ua(3)*a33

    end subroutine ComputeDNSVelocity


    ! Compute the velocity perturbation
    subroutine Ualpha(Ua, jj, kk)
      integer, intent(in) :: jj,kk
      real(dp), dimension(3), intent(out) :: Ua

      Ua(1) = sum(  bijk*Rx(:, -Ny+jj:Ny+jj, -Nz+kk:Nz+kk ) )
      Ua(2) = sum(  bijk*Ry(:, -Ny+jj:Ny+jj, -Nz+kk:Nz+kk ) )
      Ua(3) = sum(  bijk*Rz(:, -Ny+jj:Ny+jj, -Nz+kk:Nz+kk ) )

    end subroutine Ualpha

    !Update the random field (called after each time step)
    subroutine updateRandomField()
      real(dp) :: sqrt3=sqrt(3._dp)
      
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



    end subroutine updateRandomField

    ! Support functions for the DNS filter
    function btildek(k,n)
      integer,  intent(in) :: k
      real(dp), intent(in) :: n

      real(dp) :: btildek
  
      btildek = exp( -pi*0.5_dp*real(k*k,dp)/(n*n) ) 
    
    end function btildek
  
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
