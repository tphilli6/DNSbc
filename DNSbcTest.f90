program main
  use dnsbc

  integer :: niter = 500 

  integer :: jj,kk
  real(dp), dimension(3) :: vel
  real(dp), dimension(3) :: vsum
  real(dp) :: uus, vvs, wws, uvs, uws, vws 
  real(dp), dimension(3) :: vbar
  real(dp) :: uut, vvt, wwt, uvt, uwt, vwt 

  integer :: nave

!  allocate(vel(My,Mz))

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

!  call setupDNSFilter()
  call setupDNSFilter(.5_dp, .5_dp, .5_dp, .05_dp, .05_dp, .05_dp, 10, 10)

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
        call DNSVelocityPurturbation(vel, jj, kk)
        vsum = vsum + vel
        nave = nave + 1

        uus = uus + moment(vel(1),vsum(1)/nave, vel(1), vsum(1)/nave) 
        vvs = vvs + moment(vel(2),vsum(2)/nave, vel(2), vsum(2)/nave) 
        wws = wws + moment(vel(3),vsum(3)/nave, vel(3), vsum(3)/nave) 
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

    write(*,'(3f8.4, 6f9.3)') vbar, uut, vvt, wwt, uvt, uwt, vwt

    call updateRandomField()

  enddo

  print*,
  print*, 'Expected output:'
    write(*,'(3f8.4, 6f9.3)') ubar, uu, vv, ww, uv, uw, vw

contains

  function moment(u,ubar,v,vbar)
    use DNSbc, only : dp
    real(dp) :: u, ubar, v, vbar
    real(dp) :: moment
  
    moment = (u-ubar)*(v-vbar)
  
  end function

end

