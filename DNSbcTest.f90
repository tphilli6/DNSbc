program main
  use dnsbc
  use mpi
  implicit none

  integer :: niter = 667

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



  call setupDNSFilter(.5_dp, .5_dp, .5_dp, .10_dp, .10_dp, .10_dp, 10, 10)

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

  vsum=0.
  uus =0.
  vvs =0.
  wws =0.
  uvs =0.
  vws =0.
  vws =0.

    do n=1,niter

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

      !if (rank==0) then
      !  write(*,'(3f8.4, 6f9.3)') vbar, uut, vvt, wwt, uvt, uwt, vwt
      !endif
      call DNSUpdate(0.3_dp*real(n,dp))


    enddo

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


    call mpi_finalize(ierr)

contains

  function moment(u,ubar,v,vbar)
    use DNSbc, only : dp
    real(dp) :: u, ubar, v, vbar
    real(dp) :: moment

    moment = (u-ubar)*(v-vbar)

  end function

end

