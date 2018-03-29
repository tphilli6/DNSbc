program main
  use dnsbc
  use mpi
  implicit none

  integer :: niter=10

  integer :: n, jj,kk
  real(dp), dimension(3) :: vel
  real(dp), dimension(3) :: vsum
  real(dp) :: uus, vvs, wws, uvs, uws, vws
  real(dp), dimension(3) :: vbar, ubar, vp
  real(dp) :: uut, vvt, wwt, uvt, uwt, vwt
  real(dp), dimension(3,3) :: Rij

  real(dp) :: uu, vv, ww, uv, uw, vw

  integer :: nave
  integer :: rank, nproc, ierr
 
  real(dp) :: t0, t1, tstart, tfinal

   !MPI test components
   call mpi_init(ierr)
   call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr)
   call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr)

  if (rank==0) then
    print*, 'bk tilde test:'
    print*, 'Input: k=1, n=2.7'
    print*, 'Expected output: 0.8062'
    print*, '         Output: ', btildek(1, 2.7_dp)
    print*,

    print*, 'bk test:'
    print*, 'Input: k=1, n=2.7, Nx=5'
    print*, 'Expected output: 0.4906'
    print*, '         Output: ', bk(1, 2.7_dp, 5)
    print*,
  endif


  ! Fast Test
  call setupDNSFilter(0.00001_dp, .001_dp, .001_dp, .000001_dp, .005_dp/9._dp, 0.005_dp/9._dp, 10, 10, &
                      .true., .true.)
  
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

end

