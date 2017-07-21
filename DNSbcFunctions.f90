subroutine  DNSbcStart(Lx, Ly, Lz, dx, dy, dz, My, Mz)

  use DNSbc, only : dp, setupDNSFilter
  implicit none

  real(dp) :: Lx, Ly, Lz, dx, dy, dz
  integer  :: My, Mz

  call setupDNSFilter(Lx, Ly, Lz, dx, dy, dz, My, Mz)

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

  vel = Ua(:,jj,kk)

end subroutine

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

subroutine DNSUpdate(time)
  use DNSbc, only : dp, updateUalpha
  implicit none

  real(dp), intent(in) :: time

  call updateUalpha(time)

end subroutine
