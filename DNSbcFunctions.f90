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
  use DNSbc, only : yuu, uu, vv, ww, uv, uw, vw, yuu_read, yuu_constant
  implicit none

  character(*), intent(in) :: turbFile
  integer :: i, n, ierr

  if (.not.yuu_constant) then

    open(21, file=turbFile, status='old', iostat=ierr)
    read(21, *, iostat=ierr) n
    allocate( yuu(n), uu(n), vv(n), ww(n), uv(n), uw(n), vw(n) )
    do i=1,n
      read(21, '(7(E23.15))', iostat=ierr) yuu(i), uu(i), vv(i), ww(i), uv(i), uw(i), vw(i) 
    enddo
    close(21)
  
    if (ierr.ne.0) then
      write(*,'(A)') 'Error opening or reading ',turbFile
      write(*,'(A)') '--Writing example file'
  
      !Write sample file
      open(22, file='Sample-'//turbFile, status='unknown')
      write(22,'(A)') 'n [number of lines (integer)]'
      write(22,'(A)') 'y1  uu1  vv1  ww1  uv1  uw1  vw1 [(float)]'
      write(22,'(A)') 'y2  uu2  vv2  ww2  uv2  uw2  vw2 [(float)]'
      write(22,'(A)') '...'
      write(22,'(A)') 'yn  uun  vvn  wwn  uvn  uwn  vwn [(float)]'
      close(22)
  
      stop
    endif
  
    yuu_read=.true.

  else
    write(*,'(A)') 'Warning. Turbulent properties exists through setTurbProperties.'

  endif 

end subroutine


!Set turbulent flow properties
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
subroutine setTurbProperties(uuIn, vvIn, wwIn, uvIn, uwIn, vwIn)
  use DNSbc, only : dp, yuu, uu, vv, ww, uv, uw, vw, yuu_read, yuu_constant
  implicit none

  real(dp) :: uuIn, vvIn, wwIn, uvIn, uwIn, vwIn

  if (.not.yuu_read) then

    allocate( yuu(1), uu(1), vv(1), ww(1), uv(1), uw(1), vw(1) )
    uu = uuIn
    vv = vvIn
    ww = wwIn
    uv = uvIn
    uw = uwIn
    vw = vwIn

    yuu_constant = .true.
  else
    write(*,'(A)') 'Warning. Turbulent properties exists through readTurbProperties.'

  endif


end subroutine setTurbProperties

!Read velocity profile
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!unit tested
subroutine readVelProfile(velFile)
  use DNSbc, only : yvel, velprof, yvel_read, yvel_constant
  implicit none

  character(*), intent(in) :: velFile
  integer :: i, n, ierr

 
  if (.not.yvel_constant) then 
    open(21, file=velFile, status='old', iostat=ierr)
    read(21, *, iostat=ierr) n
    allocate( yvel(n), velprof(n,3) )
    do i=1,n
      read(21, '(4(E23.15))', iostat=ierr) yvel(i), velprof(i,:)
    enddo
    close(21)


    if (ierr.ne.0) then
      write(*,'(A)') 'Error opening or reading ',velFile
      write(*,'(A)') '--Writing example file'

      !Write sample file
      open(22, file='Sample-'//velFile, status='unknown')
      write(22,'(A)') 'n [number of lines (integer)]'
      write(22,'(A)') 'y1  u1  v1  w1 [(float)]'
      write(22,'(A)') 'y2  u2  v2  w2 [(float)]'
      write(22,'(A)') '...'
      write(22,'(A)') 'yn  un  vn  wn [(float)]'
      close(22)

      stop
    endif


    yvel_read=.true.

   else
    write(*,'(A)') 'Warning. Velocity profile already exists through setVelProfile.'

   endif

end subroutine

!Set velocity profile
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
subroutine setVelProfile(velProfIn)
  use DNSbc, only : dp, yvel, velProf, yvel_read, yvel_constant
  implicit none

  real(dp) :: velProfIn

  if (.not.yvel_read) then

    allocate( yvel(1), velprof(1,3) )
    velProf = velProfIn
    yvel_constant = .true.

  else
    write(*,'(A)') 'Warning. Velocity profile already exists through readVelProfile.'

  endif


end subroutine setVelProfile
    
!Get velocity @ y
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!unit tested
subroutine getVelocity(y, vel)
  use DNSbc, only : dp, yvel, velprof, interpolateToY, &
                    yvel_read, yvel_constant
  implicit none

  real(dp), intent(in) :: y
  real(dp), dimension(3), intent(out) :: vel
  
  if (yvel_read) then
    vel(1) = interpolateToY(y, yvel, velprof(:,1)) 
    vel(2) = interpolateToY(y, yvel, velprof(:,2)) 
    vel(3) = interpolateToY(y, yvel, velprof(:,3)) 

  elseif (yvel_constant) then 
    vel = velprof(1,:)

  else
    write(*,'(A)') 'Error! Velocity data not present.'
    stop
  endif

end subroutine 

!Get turbulent properties @ y
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!unit tested
subroutine getTurbProperties(y, uuOut, vvOut, wwOut, uvOut, uwOut, vwOut)
  use DNSbc, only : dp, yuu, uu, vv, ww, uv, uw, vw, interpolateToY, &
                    yuu_read, yuu_constant
      implicit none

  real(dp), intent(in) :: y
  real(dp), intent(out) :: uuOut, vvOut, wwOut, uvOut, uwOut, vwOut
  
  if (yuu_read) then
  uuOut  = interpolateToY(y, yuu, uu) 
  vvOut  = interpolateToY(y, yuu, vv) 
  wwOut  = interpolateToY(y, yuu, ww) 
  uvOut  = interpolateToY(y, yuu, uv) 
  uwOut  = interpolateToY(y, yuu, uw) 
  vwOut  = interpolateToY(y, yuu, vw) 

  elseif (yuu_constant) then
    uuOut = uu(1)
    vvOut = vv(1)
    wwOut = ww(1)
    uvOut = uv(1)
    uwOut = uw(1)
    vwOut = vw(1)

  else
    write(*,'(A)') 'Error! Turbulent data not present.'
    stop
  endif

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
