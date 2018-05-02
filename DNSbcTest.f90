program main
  use dnsbc, only : setupDNSFilter, closeDNSFilter, My, Mz, btildek, bk, &
                    interpolateToY
  use mpi

  implicit none

  integer,  parameter :: dp=selected_real_kind(15,307)
  integer :: niter=1

  integer :: n,ii,jj,kk
  real(dp), dimension(3) :: vel
  real(dp), dimension(3) :: vsum
  real(dp) :: uus, vvs, wws, uvs, uws, vws
  real(dp), dimension(3) :: vbar, ubar, vp
  real(dp) :: uut, vvt, wwt, uvt, uwt, vwt
  real(dp), dimension(3,3) :: Rij

  real(dp) :: uu, vv, ww, uv, uw, vw

  integer :: nave
  integer, parameter :: MMy=10
  integer, parameter :: MMz=10
  real(dp), dimension(MMy) :: YYmesh
  real(dp), dimension(MMz) :: ZZmesh
 
  integer :: rank, nproc, ierr
 
  real(dp) :: t0, t1, tstart, tfinal

   !MPI test components
   call mpi_init(ierr)
   call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr)
   call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr)

  if (rank==0) then
    call btildeTest
    call interpolateToYTest
    call readTurbPropertiesTest
    call readVelProfileTest
  endif


  ! Fast Test
  do ii=1,MMy
    YYmesh(ii) = real(ii,dp)*(1._dp/MMy-1)
  enddo 
  do ii=1,MMz
    ZZmesh(ii) = real(ii,dp)*(1._dp/MMz-1)
  enddo 
  call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .000001_dp, .005_dp/9._dp, 0.005_dp/9._dp, 10, 10, &
                      YYmesh, ZZmesh, .true., .true.)
  
  ! More realistic test but still fast-ish
  !call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .000001_dp, .005_dp/29._dp, 0.005_dp/29._dp, 30, 30, &
  !                    .true., .true.)

  ! Tests of scaling  ------------------------------------------------
  !call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .000001_dp, .005_dp/79._dp, 0.005_dp/79._dp, 80, 80, &
  !                    .true., .true.)

  !call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .0000005_dp, .005_dp/159._dp, 0.005_dp/159._dp, 160, 160, &
  !                    .false., .false.)

  !call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .00000025_dp, .005_dp/319._dp, 0.005_dp/319._dp, 320, 320, &
  !                    .false., .false.)

  !call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .000000125_dp, .005_dp/639._dp, 0.005_dp/639._dp, 640, 640, &
  !                    .false., .false.)
  

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
          call DNSVelocity(vel, vp, ubar, Rij)
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


      call DNSUpdate(0.000005*real(n,dp))

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
    use dnsbc, only : dp, btildek, bk

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
    write(*,'(A,f6.4)') '          Output: ', bk(1, 2.7_dp, 5)
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
    use dnsbc, only : dp, yuu, uu, vv, ww, uv, uw, vw

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
  end subroutine

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  ! readVelProfile and getVelocity test
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  subroutine readVelProfileTest
    use dnsbc, only : dp, yvel, velprof

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


  end subroutine
