subroutine  DNSbcStart(Lx, Ly, Lz, dx, dy, dz, My, Mz, Ym, Zm, pY, pZ)

  use DNSbc, only : dp, setupDNSFilter
  implicit none

  real(dp) :: Lx, Ly, Lz, dx, dy, dz
  integer  :: My, Mz
  real(dp), dimension(My) :: Ym
  real(dp), dimension(Mz) :: Zm
  logical  :: pY, pZ

  call setupDNSFilter(Lx, Ly, Lz, dx, dy, dz, My, Mz, Ym, Zm, pY, pZ)

end subroutine

subroutine DNSbcStop()
  use DNSbc, only : closeDNSFilter
  implicit none

  call closeDNSFilter()

end subroutine

!Read turbulent flow properties
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!unit tested
subroutine readTurbProperties(turbFile)
  use DNSbc, only : yuu, uu, vv, ww, uv, uw, vw, yuu_read
  implicit none

  character(*), intent(in) :: turbFile
  integer :: i, n

  open(21, file=turbFile, status='old')
  read(21, *) n
  allocate( yuu(n), uu(n), vv(n), ww(n), uv(n), uw(n), vw(n) )
  do i=1,n
    read(21, '(7(E23.15))') yuu(i), uu(i), vv(i), ww(i), uv(i), uw(i), vw(i) 
  enddo
  close(21)
  yuu_read=.true.

end subroutine


!Read velocity profile
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!unit tested
subroutine readVelProfile(velFile)
  use DNSbc, only : yvel, velprof, yvel_read
  implicit none

  character(*), intent(in) :: velFile
  integer :: i, n

  open(21, file=velFile, status='old')
  read(21, *) n
  allocate( yvel(n), velprof(n,3) )
  do i=1,n
    read(21, '(4(E23.15))') yvel(i), velprof(i,:)
  enddo
  close(21)
  yvel_read=.true.

end subroutine

    
!Get velocity @ y
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!unit tested
subroutine getVelocity(y, vel)
  use DNSbc, only : dp, yvel, velprof, interpolateToY
  implicit none

  real(dp), intent(in) :: y
  real(dp), dimension(3), intent(out) :: vel
   
  vel(1) = interpolateToY(y, yvel, velprof(:,1)) 
  vel(2) = interpolateToY(y, yvel, velprof(:,2)) 
  vel(3) = interpolateToY(y, yvel, velprof(:,3)) 

end subroutine 

!Get turbulent properties @ y
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!unit tested
subroutine getTurbProperties(y, uuOut, vvOut, wwOut, uvOut, uwOut, vwOut)
  use DNSbc, only : dp, yuu, uu, vv, ww, uv, uw, vw, interpolateToY
      implicit none

  real(dp), intent(in) :: y
  real(dp), intent(out) :: uuOut, vvOut, wwOut, uvOut, uwOut, vwOut
   
  uuOut  = interpolateToY(y, yuu, uu) 
  vvOut  = interpolateToY(y, yuu, vv) 
  wwOut  = interpolateToY(y, yuu, ww) 
  uvOut  = interpolateToY(y, yuu, uv) 
  uwOut  = interpolateToY(y, yuu, uw) 
  vwOut  = interpolateToY(y, yuu, vw) 

end subroutine 

!Get DNS velocity
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
subroutine getDNSvelocity(y, z, vel)
  use DNSbc, only : dp 
  implicit none

  real(dp), intent(in) :: y, z
  real(dp), dimension(3), intent(out) :: vel 

  real(dp), dimension(3,3) :: Rij   
  real(dp) :: uu, vv, ww, uv, uw, vw

  real(dp), dimension(3) :: Vinf, velp

  call getVelocity(y, Vinf)
  call getTurbProperties(y, uu, vv, ww, uv, uw, vw)

  Rij(1,:) = [uu, uv, uw]
  Rij(2,:) = [uv, vv, vw]
  Rij(3,:) = [uw, vw, ww]
 
  call DNSVelocityPerturbationX(velp, y, z)
  call DNSVelocity(vel, velp, Vinf, Rij)


end subroutine 


!Read velocity profile
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
subroutine DNSVelocityPerturbation(vel, jj, kk)
  use DNSbc, only : dp, Ua

  implicit none

  integer, intent(in) :: jj, kk
  real(dp), dimension(3), intent(out) :: vel

  vel = Ua(jj,kk,:)

end subroutine

subroutine DNSVelocityPerturbationX(vel, y, z)
  use DNSbc, only : dp, Ua, My, Mz, Ymesh, Zmesh

  implicit none
 
  real(dp), dimension(3), intent(out) :: vel
  real(dp),               intent(in)  :: y, z
  real(dp), dimension(1,2) :: L
  real(dp), dimension(2,1) :: R
  real(dp), dimension(2,2) :: C
  real(dp), dimension(1,1) :: vtemp

  integer :: j0, k0, n, jj, kk
  real(dp) :: y0, z0, yb, zb
  real(dp) :: dy, dz

  y0=0._dp
  z0=0._dp
  yb=0._dp
  zb=0._dp
  dy=0._dp
  dz=0._dp
  j0=1
  k0=1

  do jj=1,My-1
    if (y>=Ymesh(jj)) then
      j0=jj 
      y0=Ymesh(jj)
      dy=Ymesh(jj+1)-Ymesh(jj)
    endif
  enddo

  do kk=1,Mz-1
    if (z>=Zmesh(kk)) then
      k0=kk 
      z0=Zmesh(kk)
      dz=Zmesh(kk+1)-Zmesh(kk)
    endif
  enddo


  !j0 = max(int(y/dy),1)
  !if (j0.gt.My-1) then
  !  j0=My-1
  !endif
  !y0 = real(j0-1,dp)*dy
  yb = (y-y0)/dy

  !k0 = max(int(z/dz),1)
  !if (k0.gt.Mz-1) then
  !  k0=Mz-1
  !endif
  !z0 = real(k0-1,dp)*dz
  zb = (z-z0)/dz

  L(1,1) = 1._dp - yb
  L(1,2) = yb
  R(1,1) = 1._dp - zb
  R(2,1) = zb

  do n=1,3
    C = Ua(j0:j0+1, k0:k0+1,n) 
    vtemp = matmul(matmul(L,C), R)
    vel(n) = vtemp(1,1)
  enddo 

end subroutine DNSVelocityPerturbationX


subroutine DNSVelocity(vel, vp, velAve, Rij)
  use DNSbc, only : dp
  implicit none

  real(dp), dimension(3,3), intent(in) :: Rij
  real(dp), dimension(3),   intent(in) :: vp
  real(dp), dimension(3),   intent(in) :: velAve
  real(dp), dimension(3),   intent(out):: vel

  real(dp) :: a11, a21, a22, a31, a32, a33

  a11 = sqrt(Rij(1,1))

  a21 = Rij(2,1)/a11
  a22 = sqrt( Rij(2,2) - a21**2 )

  a31 = Rij(3,1)/a11
  a32 = ( Rij(3,2) - a21*a31 )/a22
  a33 = sqrt( Rij(3,3) - a31**2 - a32**2 )

  vel(1) = velAve(1) + vp(1)*a11
  vel(2) = velAve(2) + vp(1)*a21 + vp(2)*a22
  vel(3) = velAve(3) + vp(1)*a31 + vp(2)*a32 + vp(3)*a33


end subroutine

subroutine DNSUpdate(currenttime)
  use DNSbc, only : dp, updateUalpha
  implicit none

  real(dp), intent(in) :: currenttime

  call updateUalpha(currenttime)


end subroutine
