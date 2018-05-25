program main
  use dnsbc, only : setupDNSFilter, closeDNSFilter, My, Mz, interpolateToY
  use mpi

  implicit none

  integer, parameter  :: dp=selected_real_kind(15,307)
  real(dp), parameter :: pi = 4_dp*atan(1._dp)

  integer, parameter :: niter=100
  real(dp), parameter :: r=4.0_dp
  real(dp) :: dt=0.2_dp/r!0.003_dp ! should be at least half the smallest length scale so nx>=2
  real(dp) :: Ly=4.0_dp
  real(dp) :: Lz=4.0_dp!4.0_dp/6._dp*pi

  real(dp) :: LSx = .8_dp!0.3_dp !this should be the largest length scale
  real(dp) :: LSy = 1.6_dp!0.35_dp
  real(dp) :: LSz = 1.6_dp!0.35_dp


  integer :: n,ii,jj,kk
  real(dp), dimension(3) :: vel
  real(dp), dimension(3) :: vsum
  real(dp) :: uus, vvs, wws, uvs, uws, vws
  real(dp), dimension(3) :: vbar, ubar, vp
  real(dp) :: uut, vvt, wwt, uvt, uwt, vwt
  real(dp), dimension(3,3) :: Rij

  real(dp) :: uu, vv, ww, uv, uw, vw

  integer :: nave
  integer, parameter :: MMy=(11-1)*int(r)+1!65
  integer, parameter :: MMz=(11-1)*int(r)+1!65!115!215

  real(dp), dimension(MMy) :: LSx_test, LSz_test
  real(dp)                 :: LSy_test


  real(dp), dimension(MMy) :: YYmesh
  real(dp), dimension(MMz) :: ZZmesh

  integer :: rank, nproc, ierr

  real(dp) :: t0, t1, tstart, tfinal

  real(dp), dimension(MMy, MMz, niter) :: xvel, yvel, zvel
  real(dp), dimension(MMy, MMz, niter) :: xvelp, yvelp, zvelp


  interface
    subroutine setLengthScale(selectfile, L_uu, L_vv, L_ww, L_avg)
      use DNSbc, only : dp, yLx, yLy, yLz, Lxvec, Lyvec, Lzvec, yL_read, yL_constant

      real(dp), optional, intent(in) :: L_uu, L_vv, L_ww, L_avg
      integer,            intent(in) :: selectfile
    end subroutine
  end interface


   !MPI test components
   call mpi_init(ierr)
   call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr)
   call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr)

  if (rank==0) then
    call btildeTest
    call interpolateToYTest
    call readTurbPropertiesTest
    call readVelProfileTest
    call readLengthScaleTest
  endif


  ! Fast Test
  do ii=1,MMy
    !YYmesh(ii) = Ly*real(ii-1,dp)/real(MMy-1,dp)-1._dp
    YYmesh(ii) = cos( (MMy-ii)*pi/(MMy-1) )
  enddo
  do ii=1,MMz
    ZZmesh(ii) = Lz*real(ii-1,dp)/real(MMz-1,dp)
  enddo



  ! Expected
  ubar = (/ 2.0_dp, 0._dp, 0._dp /)
  uu = 0.1_dp
  vv = 0.03_dp
  ww = 0.01_dp
  uv = 0.02_dp
  uw = 0.01_dp
  vw = 0.015_dp
  Rij(1,:) = (/uu, uv, uw/)
  Rij(2,:) = (/uv, vv, vw/)
  Rij(3,:) = (/uw, vw, ww/)

  call setTurbProperties(uu, vv, ww, uv, uw, vw)
  call setVelProfile(ubar)
  call setLengthScale(1, L_avg=LSx)
  call setLengthScale(2, L_avg=LSy)
  call setLengthScale(3, L_avg=LSz)

  !call readTurbProperties('Channel_inflow_turbulence.txt')
  !call readVelProfile('Channel_inflow_velocity.txt')
  !call readLengthScale('Channel_LengthScaleX.txt',1)
  !call readLengthScale('Channel_LengthScaleY.txt',2)
  !!!call setLengthScale(2, L_avg=LSy)
  !call readLengthScale('Channel_LengthScaleZ.txt',3)




  call setupDNSFilter(LSx, LSy, LSz, dt, Ly/real(MMy-1,dp), Lz/real(MMz-1,dp), MMy, MMz, &
                      YYmesh, ZZmesh, .true., .true.)

  ! More realistic test but still fast-ish
  !call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .000001_dp, .005_dp/29._dp, 0.005_dp/29._dp, 30, 30, &
  !                    .true., .true.)

  ! Tests of scaling  ------------------------------------------------
  !call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .000001_dp, .005_dp/79._dp, 0.005_dp/79._dp, 80, 80, &
  !                    YYmesh, ZZmesh, .true., .true.)

  !call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .0000005_dp, .005_dp/159._dp, 0.005_dp/159._dp, 160, 160, &
  !                    YYmesh, ZZmesh, .false., .false.)

  !call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .00000025_dp, .005_dp/319._dp, 0.005_dp/319._dp, 320, 320, &
  !                    .false., .false.)

  !call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .000000125_dp, .005_dp/639._dp, 0.005_dp/639._dp, 640, 640, &
  !                    .false., .false.)




  nave=0

  vsum=0._dp
  uus =0._dp
  vvs =0._dp
  wws =0._dp
  uvs =0._dp
  uws =0._dp
  vws =0._dp

    call cpu_time(tstart)
    do n=1,niter
      call cpu_time(t0)
      do kk=1,Mz
        do jj=1,My

          !call ComputeDNSVelocity(vel, jj, kk)
          !call Ualpha(vel, jj, kk)
          call DNSVelocityPerturbation(vp, jj, kk)
          !call DNSVelocityPerturbationX(vp, YYmesh(jj), ZZmesh(kk))

          !call DNSVelocity(vel, vp, ubar, Rij)

          call getDNSvelocity(yymesh(jj), zzmesh(kk), vel)



          xvel(jj,kk,n) = vel(1)
          yvel(jj,kk,n) = vel(2)
          zvel(jj,kk,n) = vel(3)

          xvelp(jj,kk,n) = vp(1)
          yvelp(jj,kk,n) = vp(2)
          zvelp(jj,kk,n) = vp(3)

          vsum = vsum + vel
          nave = nave + 1

          uus = uus + moment(vel(1), vsum(1)/nave, vel(1), vsum(1)/nave)
          vvs = vvs + moment(vel(2), vsum(2)/nave, vel(2), vsum(2)/nave)
          wws = wws + moment(vel(3), vsum(3)/nave, vel(3), vsum(3)/nave)
          uvs = uvs + moment(vel(1),vsum(1)/nave, vel(2), vsum(2)/nave)
          uws = uws + moment(vel(1),vsum(1)/nave, vel(3), vsum(3)/nave)
          vws = vws + moment(vel(2),vsum(2)/nave, vel(3), vsum(3)/nave)
        enddo
      enddo


      vbar=vsum/nave
      uut = uus/nave
      vvt = vvs/nave
      wwt = wws/nave
      uvt = uvs/nave
      uwt = uws/nave
      vwt = vws/nave


      call DNSUpdate(dt*real(n,dp))

      call cpu_time(t1)

      if (mod(n,10)==0) &
      print*, "n=",n, "time=", t1-t0, "seconds ", "U=[", vbar,"]"


    enddo
    call cpu_time(tfinal)
    print*, 'Average time per iteration: ', (tfinal-tstart)/niter, 's'

    print*, 'Rand Test: ', rand(0)

    if (rank==0) then
      print*,
      print*, 'Expected output:'
      write(*,'(A, 3f8.4, 6f9.3)') 'Processor   :', ubar, uu, vv, ww, uv, uw, vw
      print*, '----------------------------------------------------------------'
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    write(*,'(A, I2, A1, 3f8.4, 6f9.3)') 'Processor ', rank, ':', vbar, uut, vvt, wwt, uvt, uwt, vwt



    call calc_turb_properties(yymesh, xvel, yvel, zvel, niter, MMy, MMz)

    write(77,'(e23.14)') xvelp
    call calc_length_scale(LSx_test, LSy_test, LSz_test, xvelp,niter,MMy,MMz, dt, Ly/real(MMy-1,dp), Lz/real(MMz-1,dp), yymesh, 1 )
    call calc_length_scale(LSx_test, LSy_test, LSz_test, yvelp,niter,MMy,MMz, dt, Ly/real(MMy-1,dp), Lz/real(MMz-1,dp), yymesh, 2 )
    call calc_length_scale(LSx_test, LSy_test, LSz_test, zvelp,niter,MMy,MMz, dt, Ly/real(MMy-1,dp), Lz/real(MMz-1,dp), yymesh, 3 )


    call closeDNSFilter()


    call MPI_Finalize(ierr)

contains

  function moment(u,ubar,v,vbar)
    use DNSbc, only : dp
    real(dp) :: u, ubar, v, vbar
    real(dp) :: moment

    moment = (u-ubar)*(v-vbar)

  end function

  !function interpolateToY(yin)!, yvar, var)
  !  real(dp) :: yin
  !  real(dp) :: interpolateToY
  !
  !  interpolateToY=0._dp
  !!Something is working right here...
  !  print*, yin
  !  return
  !
  !end function
end

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  ! btilde and bk Test
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  subroutine btildeTest
    use dnsbc, only : dp, btildek, bkfun

    print*, '-----------------------------------------------------------------------'
    print*, 'bk tilde test:'
    print*, '-----------------------------------------------------------------------'
    print*, 'Input: k=1, n=2.7'
    print*, 'Expected output: 0.8062'
    write(*,'(A,f6.4)') '          Output: ', btildek(1, 2.7_dp)
    print*, '-----------------------------------------------------------------------'
    print*,

    print*, '-----------------------------------------------------------------------'
    print*, 'bk test:'
    print*, '-----------------------------------------------------------------------'
    print*, 'Input: k=1, n=2.7, Nx=5'
    print*, 'Expected output: 0.4906'
    write(*,'(A,f6.4)') '          Output: ', bkfun(1, 2.7_dp, 5)
    print*, '-----------------------------------------------------------------------'
    print*,

  end subroutine btildeTest
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  ! interpolateToY Test
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  subroutine interpolateToYTest
    use dnsbc, only : dp, interpolateToY

    real(dp) :: interpolateToY2


    print*, '-----------------------------------------------------------------------'
    print*, 'interpolateToY test:'
    print*, '-----------------------------------------------------------------------'
    print*, 'Input: y=0.5, yvar=[0,1], var=[5, 10]'
    print*, 'Expected output: 7.5'
    interpolateToY2= interpolateToY(0.5_dp, [0._dp,1._dp],[5._dp,10._dp])
    write(*,'(A,f3.1)') '          Output: ', interpolateToY2
    print*, '-----------------------------------------------------------------------'
    print*,

  end subroutine interpolateToYTest

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  ! readTurbProperties and getTurbProperties test
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  subroutine readTurbPropertiesTest
    use dnsbc, only : dp, yuu, uu, vv, ww, uv, uw, vw, yuu_read, yuu_constant

    real(dp) :: uuOut, vvOut, wwOut, uvOut, uwOut, vwOut

    print*, '-----------------------------------------------------------------------'
    print*, 'readTurbProperties test:'
    print*, '-----------------------------------------------------------------------'
    open(21, file='turbFileTest.dat', status='unknown')
    write(21,*) 2
    write(21, '(7(E23.15))') 0._dp, 5._dp, 5._dp, 5._dp, 5._dp, 5._dp, 5._dp
    write(21, '(7(E23.15))') 1._dp, 9._dp, 9._dp, 9._dp, 9._dp, 9._dp, 9._dp
    close(21)

    write(*,'(A)') 'Expected output: |  y | uu | vv | ww | uv | uw | vw |'
    write(*, '(A17, 7(F5.1))') ' ', 0._dp, 5._dp, 5._dp, 5._dp, 5._dp, 5._dp, 5._dp
    write(*, '(A17, 7(F5.1))') ' ', 1._dp, 9._dp, 9._dp, 9._dp, 9._dp, 9._dp, 9._dp
    call readTurbProperties('turbFileTest.dat')

    print*,
    write(*,'(A)') '         Output: |  y | uu | vv | ww | uv | uw | vw |'
    write(*, '(A17, 7(F5.1))') ' ', yuu(1), uu(1), vv(1), ww(1), uv(1), uw(1), vw(1)
    write(*, '(A17, 7(F5.1))') ' ', yuu(2), uu(2), vv(2), ww(2), uv(2), uw(2), vw(2)

    print*, '-----------------------------------------------------------------------'
    print*,


    print*, '-----------------------------------------------------------------------'
    print*, 'getTurbProperties test:'
    print*, '-----------------------------------------------------------------------'
    call getTurbProperties(0.5_dp, uuOut, vvOut, wwOut, uvOut, uwOut, vwOut)


    write(*,'(A)') ' Expected Output: | uu | vv | ww | uv | uw | vw |'
    write(*,'(A18,6(f5.1))') ' ',7.0_dp,7.0_dp,7.0_dp,7.0_dp,7.0_dp,7.0_dp
    print*,
    write(*,'(A)') '          Output: | uu | vv | ww | uv | uw | vw |'
    write(*,'(A18,6(f5.1))') ' ',uuOut, vvOut, wwOut, uvOut, uwOut, vwOut
    print*, '-----------------------------------------------------------------------'
    print*,

    deallocate(yuu, uu, vv, ww, uv, uw, vw)
    yuu_read=.false.
    yuu_constant=.false.

  end subroutine

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  ! readVelProfile and getVelocity test
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  subroutine readVelProfileTest
    use dnsbc, only : dp, yvel, velprof, yvel_read, yvel_constant

    integer :: n=2
    real(dp), dimension(3) :: veltest

    print*, '-----------------------------------------------------------------------'
    print*, 'readVelProfile test:'
    print*, '-----------------------------------------------------------------------'

    open(21, file='velFileTest.dat', status='unknown')
    write(21,*) n
    write(21, '(4(E23.15))') 0._dp, 2._dp, 4._dp, 6._dp
    write(21, '(4(E23.15))') 1._dp, 7._dp, 8._dp, 9._dp
    close(21)

    write(*,'(A)') 'Expected output: | y | u | v | w |'
    write(*, '(A17, 7(F4.1))') ' ', 0._dp, 2._dp, 4._dp, 6._dp
    write(*, '(A17, 7(F4.1))') ' ', 1._dp, 7._dp, 8._dp, 9._dp
    call readVelProfile('velFileTest.dat')

    print*,
    write(*,'(A)') '         Output: | y | u | v | w |'
    write(*, '(A17, 7(F4.1))') ' ', yvel(1), velprof(1,1), velprof(1,2), velprof(1,3)
    write(*, '(A17, 7(F4.1))') ' ', yvel(2), velprof(2,1), velprof(2,2), velprof(2,3)

    print*, '-----------------------------------------------------------------------'
    print*,


    print*, '-----------------------------------------------------------------------'
    print*, 'getVelocity test:'
    print*, '-----------------------------------------------------------------------'
    call getVelocity(0.5_dp, veltest)


    write(*,'(A)') ' Expected Output: | u | v | w |'
    write(*,'(A18,3(f4.1))') ' ',4.5_dp, 6._dp, 7.5_dp
    print*,
    write(*,'(A)') '          Output: | u | v | w |'
    write(*,'(A18,3(f4.1))') ' ', veltest
    print*, '-----------------------------------------------------------------------'
    print*,

    deallocate( yvel, velprof)
    yvel_read=.false.
    yvel_constant=.false.

  end subroutine

  subroutine readLengthScaleTest
  use dnsbc, only : dp, yLx, yLy, yLz, Lxvec, Lyvec, Lzvec, yL_read, yL_constant

    integer :: n=2

    print*, '-----------------------------------------------------------------------'
    print*, 'readLengthScale test:'
    print*, '-----------------------------------------------------------------------'

    open(21, file='LengthScaletest.dat', status='unknown')
    write(21,*) n
    write(21, '(4(E23.15))') 0._dp, 2._dp, 4._dp, 6._dp
    write(21, '(4(E23.15))') 1._dp, 7._dp, 8._dp, 9._dp
    close(21)

    write(*,'(A)') 'Expected output: |   y   | Lx_uu | Lx_vv | Lx_ww |'
    write(*, '(A15, 4(F8.1))') ' ', 0._dp, 2._dp, 4._dp, 6._dp
    write(*, '(A15, 4(F8.1))') ' ', 1._dp, 7._dp, 8._dp, 9._dp

    call readLengthScale('LengthScaletest.dat', 1)
    print*,
    write(*,'(A)') '         Output: |   y   | Lx_uu | Lx_vv | Lx_ww |'
    write(*, '(A15, 4(F8.1))') ' ', yLx(1), Lxvec(1,1), Lxvec(1,2), Lxvec(1,3)
    write(*, '(A15, 4(F8.1))') ' ', yLx(2), Lxvec(2,1), Lxvec(2,2), Lxvec(2,3)
    write(*, *) 'yL_read for x: ',     yL_read(1)
    write(*, *) 'yL_constant for x: ', yL_constant(1)

    call readLengthScale('LengthScaletest.dat', 2)
    print*,
    write(*,'(A)') '         Output: |   y   | Ly_uu | Ly_vv | Ly_ww |'
    write(*, '(A15, 4(F8.1))') ' ', yLy(1), Lyvec(1,1), Lyvec(1,2), Lyvec(1,3)
    write(*, '(A15, 4(F8.1))') ' ', yLy(2), Lyvec(2,1), Lyvec(2,2), Lyvec(2,3)
    write(*, *) 'yL_read for y: ',     yL_read(2)
    write(*, *) 'yL_constant for y: ', yL_constant(2)

    call readLengthScale('LengthScaletest.dat', 3)
    print*,
    write(*,'(A)') '         Output: |   y   | Lz_uu | Lz_vv | Lz_ww |'
    write(*, '(A15, 4(F8.1))') ' ', yLz(1), Lzvec(1,1), Lzvec(1,2), Lzvec(1,3)
    write(*, '(A15, 4(F8.1))') ' ', yLz(2), Lzvec(2,1), Lzvec(2,2), Lzvec(2,3)
    write(*, *) 'yL_read for z: ',     yL_read(3)
    write(*, *) 'yL_constant for z: ', yL_constant(3)

    print*, '-----------------------------------------------------------------------'
    print*,

    deallocate( yLx, yLy, yLz, Lxvec, Lyvec, Lzvec)
    yL_read=.false.
    yL_constant=.false.

  end subroutine

  subroutine setLengthScaleTest
  use dnsbc, only : dp, yLx, yLy, yLz, Lxvec, Lyvec, Lzvec, yL_read, yL_constant

    print*, '-----------------------------------------------------------------------'
    print*, 'setLengthScale Components test:'
    print*, '-----------------------------------------------------------------------'

    write(*,'(A)') 'Expected output: |   y   | Lx_uu | Lx_vv | Lx_ww |'
    write(*, '(A15, 4(F8.1))') ' ', 1._dp, 2._dp, 4._dp, 6._dp

    call setLengthScale(1, 2._dp, 4._dp, 6._dp)
    print*,
    write(*,'(A)') '         Output: |   y   | Lx_uu | Lx_vv | Lx_ww |'
    write(*, '(A15, 4(F8.1))') ' ', yLx(1), Lxvec(1,1), Lxvec(1,2), Lxvec(1,3)

    call setLengthScale(1, 2._dp, 4._dp, 6._dp)
    print*,
    write(*,'(A)') '         Output: |   y   | Ly_uu | Ly_vv | Ly_ww |'
    write(*, '(A15, 4(F8.1))') ' ', yLy(1), Lyvec(1,1), Lyvec(1,2), Lyvec(1,3)

    call setLengthScale(1, 2._dp, 4._dp, 6._dp)
    print*,
    write(*,'(A)') '         Output: |   y   | Lz_uu | Lz_vv | Lz_ww |'
    write(*, '(A15, 4(F8.1))') ' ', yLz(1), Lzvec(1,1), Lzvec(1,2), Lzvec(1,3)

    write(*, *) 'yL_read for z: ',     yL_read(3)
    write(*, *) 'yL_constant for z: ', yL_constant(3)

    print*, '-----------------------------------------------------------------------'
    print*,

    deallocate( yLx, yLy, yLz, Lxvec, Lyvec, Lzvec)
    yL_read=.false.
    yL_constant=.false.

  end subroutine


  subroutine calc_turb_properties(ymesh, xvel, yvel, zvel, nx, ny, nz)
    use dnsbc, only : dp

    integer, intent(in) :: nx, ny, nz
    real(dp), dimension(ny),         intent(in) :: ymesh
    real(dp), dimension(ny, nz, nx), intent(in) :: xvel, yvel, zvel


    real(dp), dimension(ny) :: uavg, vavg, wavg, uu, vv, ww, uv, uw, vw
    real(dp), dimension(nz,nx) :: up, vp, wp

    real(dp) :: uuOut, vvOut, wwOut, uvOut, uwOut, vwOut
    real(dp), dimension(3) :: Vinf
    integer :: jj


    !Calc average velocities
    do jj=1,ny
      uavg(jj)=sum(xvel(jj,:,:))/real(nz*nx,dp)
      vavg(jj)=sum(yvel(jj,:,:))/real(nz*nx,dp)
      wavg(jj)=sum(zvel(jj,:,:))/real(nz*nx,dp)
    enddo

    do jj=1,ny
      up = xvel(jj,:,:)-uavg(jj)
      vp = yvel(jj,:,:)-vavg(jj)
      wp = zvel(jj,:,:)-wavg(jj)

      uu(jj) = sum( up*up)/real(nz*nx,dp)
      vv(jj) = sum( vp*vp)/real(nz*nx,dp)
      ww(jj) = sum( wp*wp)/real(nz*nx,dp)

      uv(jj) = sum( up*vp)/real(nz*nx,dp)
      uw(jj) = sum( up*wp)/real(nz*nx,dp)
      vw(jj) = sum( vp*wp)/real(nz*nx,dp)

    enddo


    print*, '--------------------------------------------------'
    print*, 'Turbulent statistics------------------------------'
    print*, '--------------------------------------------------'
    print*,
    print*, 'Average velocity'
    print*, '--------------------------------------------------'
    do jj=1,ny
      call getVelocity(ymesh(jj), Vinf)
      write(*,'(f10.6, 3(f10.6, A2, f9.6, A1))') ymesh(jj), uavg(jj), '[', Vinf(1), ']', &
                                                            vavg(jj), '[', Vinf(2), ']', &
                                                            wavg(jj), '[', Vinf(3), ']'
    enddo

    open(21,file='VelocityTest.txt',status='unknown')
    do jj=1,ny
      call getVelocity(ymesh(jj), Vinf)
      write(21,'(f10.6, 3(f10.6, f10.6))') ymesh(jj), uavg(jj), Vinf(1), &
                                                     vavg(jj), Vinf(2), &
                                                     wavg(jj), Vinf(3)
    enddo
    close(21)

    print*, '--------------------------------------------------'
    write(*,'(A10, 3(f10.6))') 'Average ', sum(uavg)/real(ny,dp), sum(vavg)/real(ny,dp), sum(wavg)/real(ny,dp)
    print*, '--------------------------------------------------'

    print*,
    print*, 'Second moments '
    print*, '--------------------------------------------------'
    do jj=1,ny
      call getTurbProperties(ymesh(jj), uuOut, vvOut, wwOut, uvOut, uwOut, vwOut)
      write(*,'(f10.6, 6(f10.6, A2, f9.6, A1))') ymesh(jj), uu(jj), ' [', uuOut,']', &
                                                            vv(jj), ' [', vvOut,']', &
                                                            ww(jj), ' [', wwOut,']', &
                                                            uv(jj), ' [', uvOut,']', &
                                                            uw(jj), ' [', uwOut,']', &
                                                            vw(jj), ' [', vwOut,']'
    enddo


    open(21,file='SecondMomentsTest.txt',status='unknown')
    do jj=1,ny
      call getTurbProperties(ymesh(jj), uuOut, vvOut, wwOut, uvOut, uwOut, vwOut)
      write(21,'(f10.6, 6(f10.6, f9.6))') ymesh(jj), uu(jj),  uuOut, &
                                                            vv(jj),  vvOut, &
                                                            ww(jj),  wwOut, &
                                                            uv(jj),  uvOut, &
                                                            uw(jj),  uwOut, &
                                                            vw(jj),  vwOut
    enddo
    close(21)

    print*, '--------------------------------------------------'
    write(*,'(A10, 6(f10.6))') 'Average ', sum(uu)/real(ny,dp), &
                                           sum(vv)/real(ny,dp), &
                                           sum(ww)/real(ny,dp), &
                                           sum(uv)/real(ny,dp), &
                                           sum(uw)/real(ny,dp), &
                                           sum(vw)/real(ny,dp)
    print*, '--------------------------------------------------'


  end subroutine calc_turb_properties

  subroutine calc_length_scale(Lx, Ly, Lz, vel,nx,ny,nz, dx, dy, dz, ymesh, selectVel)
  use dnsbc, only : dp, interpolateToY

    integer,                       intent(in) :: selectVel
    real(dp), dimension(ny),       intent(in) :: ymesh
    real(dp), dimension(ny,nz,nx), intent(in) :: vel
    real(dp),                      intent(in) :: dx, dy, dz

    real(dp), parameter :: pi = 4_dp*atan(1._dp)
    integer :: ii, jj, kk, j
    real(dp), dimension(ny-1) :: fy, fysum
    real(dp), dimension(nz-1) :: fz, fzsum
    real(dp), dimension(nx-1) :: fx, fxsum
    real(dp), dimension(ny)   :: vel_eq, ymesh_eq

    real(dp), dimension(ny), intent(out) :: Lx, Lz
    real(dp),                intent(out) :: Ly

    real(dp), dimension(3) :: LSx, LSy, LSz




    if (selectVel.eq.1) then
      open(22,file='Xvelocity_correlation_functionX.txt',status='unknown')
      open(23,file='Xvelocity_correlation_functionY.txt',status='unknown')
      open(24,file='Xvelocity_correlation_functionZ.txt',status='unknown')
    elseif (selectVel.eq.2) then
      open(22,file='Yvelocity_correlation_functionX.txt',status='unknown')
      open(23,file='Yvelocity_correlation_functionY.txt',status='unknown')
      open(24,file='Yvelocity_correlation_functionZ.txt',status='unknown')
    elseif (selectVel.eq.3) then
      open(22,file='Zvelocity_correlation_functionX.txt',status='unknown')
      open(23,file='Zvelocity_correlation_functionY.txt',status='unknown')
      open(24,file='Zvelocity_correlation_functionZ.txt',status='unknown')
    endif



    do jj=1,ny
      ymesh_eq(jj) = ymesh(1) + real(jj-1,dp)*dy
    enddo

    fysum=0._dp
    do ii=1,nx
      do kk=1,nz
        do jj=1,ny-1

          do j=1,ny
            vel_eq(j)=interpolateToY(ymesh_eq(j), ymesh, vel(:,kk,ii))
          enddo

          fy(jj)=sum(  vel_eq(1:ny+1-jj)*vel_eq(jj:ny) )

        enddo

        if (fy(1)<1.0e-12_dp) then
          fy = 1._dp
        endif

        fy = fy/fy(1)

        fysum = fysum + fy

      enddo
    enddo

    fysum = fysum/real( nx*nz, dp )

    do jj=1,ny
      write(23,'(f10.6)', advance='no') dy*real(jj-1,dp)
    enddo
    write(23,*) ' '

    do jj=1,ny
      write(23,'(f10.6)', advance='no') fysum(jj)
    enddo
    write(23,*) ' '

    do jj=1,ny-1
      if (fysum(jj)>exp(-pi/4) .and. fysum(jj+1)<exp(-pi/4)) then
        Ly = ( (exp(-pi/4)-fysum(jj))/(fysum(jj+1)-fysum(jj)) + (jj-1))*dy
      endif
    enddo


    !find z length scale
    write(24,'(f10.6)',advance='no') 0._dp
    do jj=1,nz-1
      write(24,'(f10.6)', advance='no') dz*real(jj-1,dp)
    enddo
    write(24,*) ' '


    do jj = 1,ny

      fzsum=0._dp
      do ii=1,nx
        do kk=1,nz-1
          fz(kk) = sum( vel(jj, 1:nz+1-kk, ii)*vel(jj,kk:nz,ii) )
        enddo

        if (fz(1)<1.0e-12_dp) then
          fz = 1._dp
        endif

        fz = fz/fz(1)
        fzsum = fzsum + fz
      enddo

      fzsum = fzsum/real( nx, dp)
      do kk=1,nz-2
        if (fzsum(kk)>exp(-pi/4) .and. fzsum(kk+1)<exp(-pi/4)) then
          Lz(jj) = ( (exp(-pi/4)-fzsum(kk))/(fzsum(kk+1)-fzsum(kk)) + (kk-1))*dz
        endif
      enddo

      write(24,'(f10.6)',advance='no') ymesh(jj)
      do j=1,ny-1
        write(24,'(f10.6)', advance='no') fzsum(j)
      enddo
      write(24,*) ' '


    enddo

    !find x length scale
    write(22,'(f10.6)',advance='no') 0._dp
    do jj=1,nx-1
      write(22,'(f10.6)', advance='no') dx*real(jj-1,dp)
    enddo
    write(22,*) ' '



    do jj = 1,ny

      fxsum=0._dp
      do kk=1,nz
        do ii=1,nx-1
          fx(ii) = sum( vel(jj, kk, 1:nx+1-ii)*vel(jj,kk,ii:nx) )
        enddo

        if (fx(1)<1.0e-12_dp) then
          fx = 1._dp
        endif

        fx = fx/fx(1)
        fxsum = fxsum + fx
      enddo

      fxsum = fxsum/real( nz, dp )
      do ii=1,nx-2
        if (fxsum(ii)>exp(-pi/4) .and. fxsum(ii+1)<exp(-pi/4)) then
          Lx(jj) = ( (exp(-pi/4)-fxsum(ii))/(fxsum(ii+1)-fxsum(ii)) + (ii-1))*dx
        endif
      enddo

      write(22,'(f10.6)',advance='no') ymesh(jj)
      do j=1,nx-1
        write(22,'(f10.6)', advance='no') fxsum(j)
      enddo
      write(22,*) ' '


    enddo


    close(22)
    close(23)
    close(24)


    print*, '--------------------------------------------------'

    if (selectVel.eq.1) then
      print*, 'X-velocity Length Scale --------------------------'
      open(21,file='Xvelocity_lengthscale.txt',status='unknown')
    elseif (selectVel.eq.2) then
      print*, 'Y-velocity Length Scale --------------------------'
      open(21,file='Yvelocity_lengthscale.txt',status='unknown')
    elseif (selectVel.eq.3) then
      print*, 'Z-velocity Length Scale --------------------------'
      open(21,file='Zvelocity_lengthscale.txt',status='unknown')
    endif
    print*, '--------------------------------------------------'
    print*,
    do jj=1,ny
      call getLengthScale(ymesh(jj), LSx, LSy, LSz)
      write(*,'(f10.6, 3(f10.6, A2, f9.6, A1))') ymesh(jj), Lx(jj), '[', LSx(selectVel), ']', &
                                                            Ly,     '[', LSy(selectVel), ']', &
                                                            Lz(jj), '[', LSz(selectVel), ']'

    enddo

    do jj=1,ny
      call getLengthScale(ymesh(jj), LSx, LSy, LSz)
      write(21,'(f10.6, 3(f10.6, f9.6))') ymesh(jj), Lx(jj), LSx(selectVel), &
                                                             Ly,     LSy(selectVel), &
                                                             Lz(jj), LSz(selectVel)
    enddo

    close(21)



    print*, '--------------------------------------------------'
    write(*,'(A10, 6(f10.6))') 'Average ', sum(Lx)/real(ny,dp), &
                                           Ly,                  &
                                           sum(Lz)/real(ny,dp)
    print*, '--------------------------------------------------'





  end subroutine calc_length_scale

