subroutine  DNSbcStart(Lx, Ly, Lz, dx, dy, dz, My, Mz, pY, pZ)

  use DNSbc, only : dp, setupDNSFilter
  implicit none

  real(dp) :: Lx, Ly, Lz, dx, dy, dz
  integer  :: My, Mz
  logical  :: pY, pZ

  call setupDNSFilter(Lx, Ly, Lz, dx, dy, dz, My, Mz, pY, pZ)

end subroutine

subroutine DNSbcStop()
  use DNSbc, only : closeDNSFilter
  implicit none

  call closeDNSFilter()

end subroutine

subroutine DNSVelocityPerturbation(vel, jj, kk)
  use DNSbc, only : dp, Ua

  implicit none

  integer, intent(in) :: jj, kk
  real(dp), dimension(3), intent(out) :: vel

  vel = Ua(jj,kk,:)

end subroutine

subroutine DNSVelocityPerturbationX(vel, y, z)
  use DNSbc, only : dp, dy, dz, Ua, My, Mz

  implicit none
 
  real(dp), dimension(3), intent(out) :: vel
  real(dp),               intent(in)  :: y, z
  real(dp), dimension(1,2) :: L
  real(dp), dimension(2,1) :: R
  real(dp), dimension(2,2) :: C
  real(dp), dimension(1,1) :: vtemp

  integer :: j0, k0, n, jj, kk
  real(dp) :: y0, y00, z0, z00, yb, zb

  !y00=0._dp 
  !z00=0._dp

  !do jj=0,My-1
  !  if (y>=dy*jj) then
  !    yb=(y-dy*jj)/dy
  !    j0=jj
  !  endif
  !enddo

  !do kk=0,Mz-1
  !  if (z>=dz*kk) then
  !    zb=(z-dz*kk)/dz
  !    k0=kk
  !  endif
  !enddo


  j0 = max(int(y/dy),1)
  if (j0.gt.My-1) then
    j0=My-1
  endif
  y0 = real(j0-1,dp)*dy
  yb = (y-y0)/dy

  k0 = max(int(z/dz),1)
  if (k0.gt.Mz-1) then
    k0=Mz-1
  endif
  z0 = real(k0-1,dp)*dz
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
