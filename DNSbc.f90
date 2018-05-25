  module DNSbc

    implicit none

    integer,  parameter :: dp=selected_real_kind(15,307)
    real(dp), parameter :: pi=4._dp*atan(1._dp)

    ! User supplied input length scale and mesh spacing
    real(dp) :: Lx, Ly, Lz
    real(dp) :: dx, dy, dz

    integer :: My, Mz ! Mesh inflow plane dimensions
    integer :: Nx, Ny, Nz

    real(dp), allocatable, dimension(:,:,:) :: Rx, Ry, Rz, Ua
    real(dp), allocatable, dimension(:,:)   :: bjk !Filter coefficients
    real(dp), allocatable, dimension(:) :: bii, bjj, bkk

    real(dp), allocatable, dimension(:) :: Ymesh, Zmesh

    real(dp) :: currentTime

    logical, allocatable, dimension(:,:) :: updateR
    real(dp), allocatable, dimension(:,:) :: timeArray, dtArray !Filter coefficients

    ! Used for 3rd order interpolation for non equal time steps
    real(dp), allocatable, dimension(:,:,:) :: Uat1, Uat2, Uat3

    !if periodic, random field is also periodic
    logical :: pY=.false.
    logical :: pZ=.false.

    real(dp), allocatable, dimension(:)   :: y, z !Filter coefficients

    integer :: fidout = 101

    ! Variables for storage of turbulent quantities as a function of y
    real(dp), allocatable, dimension(:) :: yuu, uu, vv, ww, uv, uw, vw
    logical :: yuu_read=.false.
    logical :: yuu_constant=.false.

    ! variables for storage of velocity profile as a function of y
    real(dp), allocatable, dimension(:) :: yvel
    real(dp), allocatable, dimension(:,:) :: velprof
    logical :: yvel_read=.false.
    logical :: yvel_constant=.false.

    ! Variables for storage of turbulent lengthscale storage as a function of y
    real(dp), allocatable, dimension(:) :: yLx, yLy, yLz
    real(dp), allocatable, dimension(:,:) :: Lxvec
    real(dp), allocatable, dimension(:,:) :: Lyvec
    real(dp), allocatable, dimension(:,:) :: Lzvec

    real(dp), allocatable, dimension(:,:) :: Lxvec_Ymesh
    real(dp), allocatable, dimension(:,:) :: Lyvec_Ymesh
    real(dp), allocatable, dimension(:,:) :: Lzvec_Ymesh

    logical, dimension(3) :: yL_read=.false.
    logical, dimension(3) :: yL_constant=.false.


contains



    !Initialize the DNS inflow filter
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine setupDNSFilter(LLx, LLy, LLz, &
                              ddx, ddy, ddz, &
                              MMy, MMz, &
                              YYmesh, ZZmesh, &
                              ppY, ppZ)
    use mpi


      real(dp), intent(in) :: LLx, LLy, LLz, ddx, ddy, ddz
      integer, intent(in)  :: MMy, MMz
      real(dp), dimension(MMy), intent(in) :: YYmesh
      real(dp), dimension(MMz), intent(in) :: ZZmesh
      logical, intent(in)  :: ppY, ppZ


      integer  :: ii, jj, kk
      real(dp) :: nnx, nny, nnz
      real(dp) :: sqrt3=sqrt(3._dp)
      real(dp) :: ds

      integer :: rank, ierr
      real(dp) :: totalSize
      real(dp) :: start, finish, start0

      !integer, allocatable, dimension(:) :: seed


      call cpu_time(start)
      start0=start

      call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr)

      !if (rank==0) then
      !print*, 'Ymesh:', YYmesh
      !print*, 'Zmesh:', ZZmesh
      !endif

      Lx=LLx
      Ly=LLy
      Lz=LLz
      dx=ddx
      dy=ddy
      dz=ddz
      My=MMy
      Mz=MMz
      pY=ppY
      pZ=ppZ


      ! Compute filter width
      nnx = Lx/dx
      nny = Ly/dy
      nnz = Lz/dz
      Nx = int(ceiling(2*nnx))
      Ny = int(ceiling(2*nny))
      Nz = int(ceiling(2*nnz))



      allocate( Ymesh(-Ny+1:Ny+My), Zmesh(-Nz+1:Nz+Mz) )

      Ymesh(1:My)=YYmesh

      !Mirror spacing outside of the background mesh in the Y direction
      do jj = 1,Ny
        ds = Ymesh(jj+1)-Ymesh(1)
        Ymesh(1-jj) = Ymesh(1)-ds

        ds = Ymesh(My) - Ymesh(My-jj)
        Ymesh(My+jj) = Ymesh(My)+ds
      enddo

      Zmesh(1:Mz)=ZZmesh
      !Mirror spacing outside of the background mesh in the Z direction
      do kk = 1,Nz
        ds = Zmesh(kk+1)-Zmesh(1)
        Zmesh(1-kk) = Zmesh(1)-ds

        ds = Zmesh(Mz) - Zmesh(Mz-kk)
        Zmesh(Mz+kk) = Zmesh(Mz)+ds
      enddo


      if (rank==0) then
        allocate( Rx(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), &
                  Ry(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), &
                  Rz(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), &
                  updateR(-Ny+1:Ny+My, -Nz+1:Nz+Mz),    &
                  timeArray(-Ny+1:Ny+My, -Nz+1:Nz+Mz),  &
                  dtArray(-Ny+1:Ny+My, -Nz+1:Nz+Mz),    &
                  Ua(My, Mz, 3),                        &
                  Uat1(My, Mz, 3),                      &
                  Uat2(My, Mz, 3),                      &
                  Uat3(My, Mz, 3),                      &
                  bjk(-Ny:Ny, -Nz:Nz),                  &
                  bii(-Nx:Nx),                          &
                  bjj(-Ny:Ny),                          &
                  bkk(-Nz:Nz),                          &
                  Lxvec_Ymesh(My,3),                    &
                  Lyvec_Ymesh(My,3),                    &
                  Lzvec_Ymesh(My,3))

        call srand(20)

        write(*,'(A)') 'Initializing DNS bc module...'
        write(*,'(A)') '-----------------------------------'
        write(*,'(A)') 'Size of Padding:'
        write(*,'(A,e12.4,A,e12.4,A,e12.4)') 'Lx=',Lx, ', Ly=',Ly, ', Lz=',Lz
        write(*,'(A,e12.4,A,e12.4,A,e12.4)') 'Dx=',dx, ', Dy=',dy, ', Dz=',dz
        write(*,'(A,f5.0,A,f5.2,A,f5.0)') 'nx=',nnx, ', ny=',nny, ', nz=',nnz
        write(*,'(A,I5.0,A,I5.0,A,I5.0)') 'Nx=',Nx, ', Ny=',Ny, ', Nz=',Nz
        print*,
        write(*,'(A)') 'Size of background mesh:'
        write(*,'(A,I4.0,A,I4.0)') 'My=',My, ', Mz=',Mz
        print*,'Size of Rx = ', Size(Rx)
        print*,

        open(fidout,file='Velocity.dat')

        !Estimate size in bytes ------------------------------------------------------
        totalSize= real(3*size(Rx)+size(timeArray) &
                     +size(dtArray)+4*size(Ua)+size(bjk)+3*size(bii),dp)/2._dp**17 &
                   +real(size(updateR),dp)/2._dp**20

        if (totalSize<=1) then
          write(*,'(A,F5.1,A)') 'Estimated Memory Usage: ',real(totalSize,dp)*2**10, ' KB'
        elseif (totalSize<=1000) then
          write(*,'(A,F5.1,A)') 'Estimated Memory Usage: ',real(totalSize,dp), ' MB'
        else
          write(*,'(A,F4.2,A)') 'Estimated Memory Usage: ',real(totalSize,dp)/2**10, ' GB'
        endif
        !-----------------------------------------------------------------------------
        write(*,'(A)') '-----------------------------------'


        do kk=-Nz+1,Nz+Mz
          do jj=-Ny+1,Ny+My
            do ii=-Nx,Nx
              Rx(ii,jj,kk) = sqrt3*(2._dp*rand(0) - 1._dp)
              Ry(ii,jj,kk) = sqrt3*(2._dp*rand(0) - 1._dp)
              Rz(ii,jj,kk) = sqrt3*(2._dp*rand(0) - 1._dp)
            enddo
          enddo
        enddo

        updateR=.true.
        timeArray=0._dp
        currentTime=0._dp

        call periodic(Rx)
        call periodic(Ry)
        call periodic(Rz)

        !Compute filter coefficients
        do ii=-Nx,Nx
          bii(ii) = bkfun(ii,nnx,Nx)
        enddo
        do jj=-Ny,Ny
          bjj(jj) = bkfun(jj,nny,Ny)
        enddo
        do kk=-Nz,Nz
          bkk(kk) = bkfun(kk,nnz,Nz)
        enddo

        do kk=-Nz,Nz
          do jj=-Ny,Ny
              bjk(jj,kk) = bjj(jj)*bkk(kk)
          enddo
        enddo

        !----------------------------------------------------------------------------
        ! Set up lengthscale arrays as a function of y
        !if an array is given, then interpolate to the Ymesh
        if (yL_read(1)) then
          do ii=1,My
            do jj=1,3
              Lxvec_Ymesh(ii,jj) = interpolateToY(Ymesh(ii), yLx, Lxvec(:,jj))
            enddo
          enddo

        !if the length scales are constant then copy to appropriate array
        elseif (yL_constant(1)) then
          do jj=1,3
            Lxvec_Ymesh(:,jj) = Lxvec(1,jj)
          enddo

        !default is to use what's passed in
        else
          Lxvec_Ymesh=Lx

        endif

        !if an array is given, then interpolate to the Ymesh
        if (yL_read(2)) then
          do ii=1,My
            do jj=1,3
              Lyvec_Ymesh(ii,jj) = interpolateToY(Ymesh(ii), yLy, Lyvec(:,jj))
            enddo
          enddo

        !if the length scales are constant then copy to appropriate array
        elseif (yL_constant(2)) then
          do jj=1,3
            Lyvec_Ymesh(:,jj) = Lyvec(1,jj)
          enddo

        else
        !default is to use what's passed in
          Lyvec_Ymesh=Ly

        endif

        !if an array is given, then interpolate to the Ymesh
        if (yL_read(3)) then
          do ii=1,My
            do jj=1,3
              Lzvec_Ymesh(ii,jj) = interpolateToY(Ymesh(ii), yLz, Lzvec(:,jj))
            enddo
          enddo

        !if the length scales are constant then copy to appropriate array
        elseif (yL_constant(3)) then
          do jj=1,3
            Lzvec_Ymesh(:,jj) = Lzvec(1,jj)
          enddo

        !default is to use what's passed in
        else
          Lzvec_Ymesh=Lz

        endif
        !----------------------------------------------------------------------------



        ! Intialize Ua
        Ua=0._dp
        Uat1=0._dp
        Uat2=0._dp
        Uat3=0._dp

        ! t=0
        !call Ualpha(Uat3)
        call Ualpha_with_variable_LSyz(Uat3)
        call updateRandomField()

        !call Ualpha(Uat2)
        call Ualpha_with_variable_LSyz(Uat2)
        call updateRandomField()

        !call Ualpha(Uat1)
        call Ualpha_with_variable_LSyz(Uat1)

        Ua=Uat1


        call cpu_time(finish)
        write(*,*) 'Initialization time: ', finish-start0, ' seconds'
      else

        allocate( Rx(1,1,1), Ry(1,1,1), Rz(1,1,1), &
                  updateR(1,1), timeArray(1,1),    &
                  Ua(My, Mz, 3), Uat1(1,1,1), Uat2(1,1,1), Uat3(1,1,1), &
                  bjk(1,1), bii(1), bjj(1), bkk(1) )
        Rx=0._dp
        Ry=0._dp
        Rz=0._dp
        updateR=.false.
        timeArray=0._dp
        Ua = 0._dp
        Uat1 = 0._dp
        Uat2 = 0._dp
        Uat3 = 0._dp
        bjk  = 0._dp
        bii = 0._dp
        bjj = 0._dp
        bkk = 0._dp

      endif


!      allocate( y(My), &
!                z(Mz) )
!
      !Testing code for variable background mesh size


      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_Bcast( Ua, 3*My*Mz, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)


    end subroutine setupDNSFilter

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine closeDNSFilter()

      deallocate( Rx,   &
                  Ry,   &
                  Rz,   &
                  updateR,   &
                  timeArray, &
                  Ua,   &
                  Uat1, &
                  Uat2, &
                  Uat3, &
                  bjk,  &
                  bii,  &
                  bjj,  &
                  bkk,  &
                  Ymesh,&
                  Zmesh,&
                  Lxvec_Ymesh, &
                  Lyvec_Ymesh, &
                  Lzvec_Ymesh &
                  )

      if (yuu_read.or.yuu_constant) deallocate(yuu, uu, vv, ww, uv, uw, vw)
      if (yvel_read.or.yvel_constant) deallocate(yvel, velprof)

    endsubroutine closeDNSFilter

    subroutine Ualpha(Va)
      implicit none

      real(dp), dimension(My,Mz,3), intent(out) :: Va

      call bijk_times_R( Va(:,:,1), Rx(-Nx:Nx,:,:), bii, bjj, bkk )
      call bijk_times_R( Va(:,:,2), Ry(-Nx:Nx,:,:), bii, bjj, bkk )
      call bijk_times_R( Va(:,:,3), Rz(-Nx:Nx,:,:), bii, bjj, bkk )

    end subroutine Ualpha


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Evaluate the velocity perturbation for all velocity components
    subroutine Ualpha_with_variable_LSy(Va)
      implicit none

      real(dp), dimension(My,Mz,3), intent(out) :: Va
      integer :: ii

      ii=1
      call Velocity_wih_variable_LSy(Va(:,:,ii), Rx, Lxvec_Ymesh(:,ii), &
                                                    Lyvec_Ymesh(:,ii), &
                                                    Lzvec_Ymesh(:,ii))
      ii=2
      call Velocity_wih_variable_LSy(Va(:,:,ii), Ry, Lxvec_Ymesh(:,ii), &
                                                    Lyvec_Ymesh(:,ii), &
                                                    Lzvec_Ymesh(:,ii))
      ii=3
      call Velocity_wih_variable_LSy(Va(:,:,ii), Rz, Lxvec_Ymesh(:,ii), &
                                                    Lyvec_Ymesh(:,ii), &
                                                    Lzvec_Ymesh(:,ii))

      !stop
    end subroutine Ualpha_with_variable_LSy

    subroutine Ualpha_with_variable_LSyz(Va)
      implicit none

      real(dp), dimension(My,Mz,3), intent(out) :: Va
      integer :: ii

      ii=1
      call Velocity_wih_variable_LSyz(Va(:,:,ii), Rx, Lxvec_Ymesh(:,ii), &
                                                    Lyvec_Ymesh(:,ii), &
                                                    Lzvec_Ymesh(:,ii))
      ii=2
      call Velocity_wih_variable_LSyz(Va(:,:,ii), Ry, Lxvec_Ymesh(:,ii), &
                                                    Lyvec_Ymesh(:,ii), &
                                                    Lzvec_Ymesh(:,ii))
      ii=3
      call Velocity_wih_variable_LSyz(Va(:,:,ii), Rz, Lxvec_Ymesh(:,ii), &
                                                    Lyvec_Ymesh(:,ii), &
                                                    Lzvec_Ymesh(:,ii))

      !stop
    end subroutine Ualpha_with_variable_LSyz

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Evaluate the velocity perturbation for a single velocity component
    subroutine Velocity_wih_variable_LSyz(Va, Rijk, LSx, LSy, LSz)
      implicit none

      real(dp), dimension(My,Mz), intent(out) :: Va
      real(dp), dimension(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), intent(in) :: Rijk
      real(dp), dimension(My),    intent(in) :: LSx, LSy, LSz

      real(dp), dimension(-Nx:Nx) :: bbit, bsi
      real(dp), dimension(-Ny:Ny) :: bbjt, bsj
      real(dp), dimension(-Nz:Nz) :: bbkt, bsk

      real(dp), dimension(My, Mz, -Nx:Nx) :: bbi
      real(dp), dimension(My, Mz, -Ny:Ny) :: bbj
      real(dp), dimension(My, Mz, -Nz:Nz) :: bbk

      integer :: ii, jj, kk, mm, nn, i, j, k, jl, jh, jjl, jjh, kl, kh, kkl, kkh

      !real(dp), dimension(My, Mz, -Ny+1:Ny+My, -Nz+1:Nz+Mz) :: Rjk
      !real(dp), dimension(My, Mz, My,-Nz+1:Nz+Mz) :: Rk
      real(dp), dimension(-Nx:Nx) :: Ri


      bbi=0._dp
      bbj=0._dp
      bbk=0._dp

      do mm=1,My

        ! For a Random array
        ! These loops will generate the filter coefficients


        if (LSx(mm)>1.0e-12_dp) then
          do ii=-Nx,Nx
            bbit(ii) = btildek2(dx, ii, LSx(mm))
          enddo

          do nn=1,Mz
              bbi(mm,nn,:) = bbit/sqrt( sum( bbit*bbit)  )
          enddo

        endif

        if (LSy(mm)>1.0e-12_dp) then
          do jj=-Ny,Ny
            bbjt(jj) = btildek3(Ymesh(mm+jj)-Ymesh(mm), LSy(mm))
          enddo

          do nn=1,Mz
            bbj(mm,nn,:) = bbjt/sqrt( sum( bbjt*bbjt)  )
          enddo

        endif

        if (LSz(mm)>1.0e-12_dp) then
          do kk=-Nz,Nz
            bbkt(kk) = btildek3(Zmesh(mm+kk)-Zmesh(mm), LSz(mm))
          enddo

          do nn=1,Mz
            bbk(mm,nn,:) = bbkt/sqrt( sum( bbkt*bbkt)  )
          enddo

        endif

      enddo


          !write(*,'(I4, 1f9.6, f9.5, f6.1 )') mm
          !write(*,'(A4, 1f9.6, f9.5, f6.1 )') ' ',  dx, LSx(mm), LSx(mm)/dx
          !write(*,'(A4, 1f9.6, f9.5, f6.1 )') ' ',  DYMesh(mm), LSy(mm), LSy(mm)/DYmesh(mm)
          !write(*,'(A4, 1f9.6, f9.5, f6.1 )') ' ',  DZmesh(mm), LSz(mm), LSz(mm)/DZmesh(mm)
          !print*,


       !do nn = 1,Mz

       !  do kk = -Nz,Nz
       !     do jj = -Ny,Ny

       !       Ri = Rijk(:,1+jj:My+jj,kk)

       !       do mm = 1,My
       !         Rk(mm) = Rjk(jj,kk) + bbj(jj,mm,nn)*sum(bbi(:,mm,nn)*Ri(:,mm) )
       !       enddo

       !     enddo




       !   end do



       ! enddo

     Va = 0._dp
     do kk=-Nz+1,Mz+Nz
       do jj=-Ny+1,My+Ny

         Ri = Rijk(:,jj,kk)

         jl = max(jj-Ny, 1)
         jh = min(jj+Ny, My)
         kl = max(kk-Nz, 1)
         kh = min(kk+Nz, Mz)


         do k = kl,kh
         do j = jl,jh

           bsi = bbi(j,k,:)
           bsj = bbj(j,k,:)
           bsk = bbk(j,k,:)
         do i = -Nx,Nx
           jjl = jj-j
           jjh = j-jh
           kkl = kk-k
           kkh = k-kh

           Va(j,k) = Va(j,k) + bsk(kk-k)*bsj(jj-j)*bsi(i)*Ri(i)

         enddo
         enddo
         enddo


       enddo

     enddo





          !call bijk_times_Ryz( Va(mm,nn), Rijk, mm, nn, bbi(:,mm), bbj(:,mm), bbk(:,mm) )



    end subroutine Velocity_wih_variable_LSyz


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Evaluate the velocity perturbation
    subroutine bijk_times_Ryz(Vel, Rijk, j, k, bi, bj, bk)
      implicit none

      real(dp), dimension(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), intent(in) :: Rijk
      real(dp), dimension(-Nx:Nx), intent(in) :: bi
      real(dp), dimension(-Ny:Ny), intent(in) :: bj
      real(dp), dimension(-Nz:Nz), intent(in) :: bk
      integer,                     intent(in) :: j, k
      real(dp),                    intent(out) :: Vel

      real(dp), dimension(-Ny:Ny, -Nz:Nz) :: Rjk
      real(dp), dimension(-Nz:Nz) :: Rk

      integer :: jj, kk

      Vel=0._dp


        ! Reduce Rijk to a 2D array
        do kk=-Nz+k,Nz+k
          do jj=-Ny+j,Ny+j
            Rjk(jj-j,kk-k) = sum(bi*Rijk(:,jj,kk))
          enddo
        enddo

        ! Multiply by bjj
        do kk=-Nz,Nz
          Rk(kk) = sum(  bj*Rjk(:, kk ) )
        enddo

        ! Multiply by bkk
        Vel = sum( bk*Rk )


    end subroutine bijk_times_Ryz


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Evaluate the velocity perturbation for a single velocity component
    subroutine Velocity_wih_variable_LSy(Va, Rijk, LSx, LSy, LSz)
      implicit none

      real(dp), dimension(My,Mz), intent(out) :: Va
      real(dp), dimension(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), intent(in) :: Rijk
      real(dp), dimension(My),    intent(in) :: LSx, LSy, LSz

      real(dp), dimension(-Nx:Nx) :: bbit
      real(dp), dimension(-Ny:Ny) :: bbjt
      real(dp), dimension(-Nz:Nz) :: bbkt

      real(dp), dimension(-Nx:Nx, My) :: bbi
      real(dp), dimension(-Ny:Ny, My) :: bbj
      real(dp), dimension(-Nz:Nz, My) :: bbk

      integer :: ii, jj, kk, mm

      bbi=0._dp
      bbj=0._dp
      bbk=0._dp

      do mm=1,My

        ! For a Random array
        ! These loops will generate the filter coefficients


        if (LSx(mm)>1.0e-12_dp) then
          do ii=-Nx,Nx
            bbit(ii) = btildek2(dx, ii, LSx(mm))
          enddo
          bbi(:,mm) = bbit/sqrt( sum( bbit*bbit)  )

        endif

        if (LSy(mm)>1.0e-12_dp) then
          do jj=-Ny,Ny
            bbjt(jj) = btildek3(Ymesh(mm+jj)-Ymesh(mm), LSy(mm))
          enddo
          bbj(:,mm) = bbjt/sqrt( sum( bbjt*bbjt)  )
        endif

        if (LSz(mm)>1.0e-12_dp) then
          do kk=-Nz,Nz
            bbkt(kk) = btildek3(Zmesh(mm+kk)-Zmesh(mm), LSz(mm))
          enddo
          bbk(:,mm) = bbkt/sqrt( sum( bbkt*bbkt)  )
        endif

      enddo


          !write(*,'(I4, 1f9.6, f9.5, f6.1 )') mm
          !write(*,'(A4, 1f9.6, f9.5, f6.1 )') ' ',  dx, LSx(mm), LSx(mm)/dx
          !write(*,'(A4, 1f9.6, f9.5, f6.1 )') ' ',  DYMesh(mm), LSy(mm), LSy(mm)/DYmesh(mm)
          !write(*,'(A4, 1f9.6, f9.5, f6.1 )') ' ',  DZmesh(mm), LSz(mm), LSz(mm)/DZmesh(mm)
          !print*,

      do mm = 1,My
        call bijk_times_Ry( Va(mm,:), Rijk, mm, bbi(:,mm), bbj(:,mm), bbk(:,mm) )
      enddo


    !print*,
    !print*,
    !print*,

    end subroutine Velocity_wih_variable_LSy


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Evaluate the velocity perturbation
    subroutine bijk_times_Ry(Vel, Rijk, j, bi, bj, bk)
      implicit none

      real(dp), dimension(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), intent(in) :: Rijk
      real(dp), dimension(-Nx:Nx), intent(in) :: bi
      real(dp), dimension(-Ny:Ny), intent(in) :: bj
      real(dp), dimension(-Nz:Nz), intent(in) :: bk
      integer,                     intent(in) :: j
      real(dp), dimension(Mz), intent(out) :: Vel

      real(dp), dimension(-Ny:Ny, -Nz+1:Nz+Mz) :: Rjk
      real(dp), dimension(-Nz+1:Nz+Mz) :: Rk

      integer :: jj, kk

      Vel=0._dp


        ! Reduce Rijk to a 2D array
        do kk=-Nz+1,Nz+Mz
          do jj=-Ny+j,Ny+j
            Rjk(jj-j,kk) = sum(bi*Rijk(:,jj,kk))
          enddo
        enddo

        ! Multiply by bjj
        do kk=-Nz+1,Nz+Mz
          Rk(kk) = sum(  bj*Rjk(:, kk ) )
        enddo

        ! Multiply by bkk
        do kk=-Nz,Nz
          Vel(:) = Vel(:) + bk(kk)*Rk(1+kk:Mz+kk)
        enddo


    end subroutine bijk_times_Ry


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Evaluate the velocity perturbation
    subroutine bijk_times_R(Vel, Rijk, bi, bj, bk)
      implicit none

      real(dp), dimension(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), intent(in) :: Rijk
      real(dp), dimension(-Nx:Nx), intent(in) :: bi
      real(dp), dimension(-Ny:Ny),        intent(in) :: bj
      real(dp), dimension(-Nz:Nz),        intent(in) :: bk
      real(dp), dimension(My,Mz), intent(out) :: Vel

      real(dp), dimension(-Ny+1:Ny+My, -Nz+1:Nz+Mz) :: Rjk
      real(dp), dimension(My,-Nz+1:Nz+Mz) :: Rk

      integer :: jj, kk

      Vel=0._dp

      ! Reduce Rijk to a 2D array
      do kk=-Nz+1,Nz+Mz
        do jj=-Ny+1,Ny+My
          Rjk(jj,kk) = sum(bi*Rijk(-Nx:Nx,jj,kk))
        enddo
      enddo

      ! Multiply by bjj
      do kk=-Nz+1,Nz+Mz
        do jj=1,My
          Rk(jj,kk) = sum(  bj*Rjk(-Ny+jj:Ny+jj, kk ) )
        enddo
      enddo

      ! Multiply by bkk
      do kk=-Nz,Nz
        Vel(:,:) = Vel(:,:) + bk(kk)*Rk(:,kk+1:Mz+kk)
      enddo



    end subroutine bijk_times_R


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Evaluate the velocity perturbation (slowish algorithm)
    subroutine bijk_times_R2(Vel, Rijk, bi)
      implicit none

      real(dp), dimension(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), intent(in) :: Rijk
      real(dp), dimension(-Nx:Nx), intent(in) :: bi
      real(dp), dimension(My,Mz), intent(out) :: Vel

      real(dp), dimension(-Ny+1:Ny+My, -Nz+1:Nz+Mz) :: Rjk

      integer :: jj, kk

      Vel=0._dp

      ! Reduce Rijk to a 2D array
      do kk=-Nz+1,Nz+Mz
        do jj=-Ny+1,Ny+My
          Rjk(jj,kk) = sum(bi*Rijk(-Nx:Nx,jj,kk))
        enddo
      enddo

      ! Multiply by bjj*bkk
      do kk=1,Mz
        do jj=1,My
          Vel(jj,kk) = sum(  bjk*Rjk(-Ny+jj:Ny+jj, -Nz+kk:Nz+kk ) )
        enddo
      enddo

    end subroutine bijk_times_R2

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Evaluate the velocity purterbation (slowest algorithm. Go make some coffee...)
    subroutine bijk_times_R3(Vel, Rijk, bi)
      implicit none

      real(dp), dimension(-Nx:Nx, -Ny+1:Ny+My, -Nz+1:Nz+Mz), intent(in) :: Rijk
      real(dp), dimension(-Nx:Nx), intent(in) :: bi
      real(dp), dimension(My,Mz), intent(out) :: Vel

      integer :: ii, jj, kk

      Vel=0._dp
      ! Multiply by bii*bjj*bkk
      do kk=1,Mz
        do jj=1,My
          do ii=-Nx,Nx
            Vel(jj,kk) = Vel(jj,kk) + sum(  bjk*bi(ii)*Rijk(ii, -Ny+jj:Ny+jj, -Nz+kk:Nz+kk ) )
          enddo
        enddo
      enddo


    end subroutine bijk_times_R3

    !Interpolate to y
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    function interpolateToY(y, yvar, var)
      !use DNSbc, only : dp
      !implicit none

      real(dp) :: y
      real(dp), dimension(:) :: yvar, var
      real(dp) :: interpolateToY

      integer :: i, n
      real(dp) :: m

      interpolateToY=0._dp
      n = size(yvar)

      do i=1,n-1
        if ((y>=yvar(i)).and.(y<=yvar(i+1)))then
          m = ( var(i+1)-var(i) )/(yvar(i+1)-yvar(i))
          interpolateToY = m*(y-yvar(i)) + var(i)
          return
        endif
      enddo

    end function

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Update the velocity perturbation array
    subroutine updateUalpha(time)
      use mpi

      implicit none

      real(dp), intent(in) :: time

      integer  :: rank, ierr
      real(dp) :: dt
      !real(dp), dimension(3,3) :: Rij
      !real(dp), dimension(My,Mz,3) :: vel
      !real(dp), dimension(3) :: velave
      real(dp), dimension(My, Mz, 3) :: Uattest
      real(dp) :: t0, t1

      call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr)



      if (rank==0) then

        call cpu_time(t0)

        dt = dx ! just for clarity
        do while (time>currentTime)

          !currentTime=time
           currentTime= currentTime + dt


          !nUpdated=1
          !do while (nUpdated.gt.0)
          !  updateR=.false.

          !  where (time.gt.timeArray)
          !    updateR=.true.
          !    timeArray = timeArray+dtArray
          !  end where

          ! call updateRandomField(updateR, nUpdated)
          !end do

         call updateRandomField()

         !update if periodic
         call periodic(Rx)
         call periodic(Ry)
         call periodic(Rz)


          Uat3=Uat2
          Uat2=Uat1
          call Ualpha_with_variable_LSyz(Uat1)
          !call Ualpha_with_variable_LSy(Uat1)
          !call Ualpha(Uat1)
          !write(*,'(11(f12.8))') Uattest(:,:,1)
          !print*,
          !write(*,'(11(f12.8))') Uat1(:,:,1)
          !print*,
          !print*, maxval( abs( Uattest - Uat1) )
          !stop

          call interpolateVelocity(time)
          !Ua=Uat1


          !print*, 'shape Rx: ', shape(Rx)
          !print*, 'shape Ua: ', shape(Ua)
          !write(fidout,'(65536(E23.15))') Ua
          !write(88,'(10201(E23.15))') vel
          !write(89,'(324(E23.15))') Rx(1,:,:)
          !write(89,'(324(E23.15))') Ry(1,:,:)
          !write(89,'(324(E23.15))') Rz(1,:,:)
          !write(*,'(A)') 'DNSbc: Interpolated velocity to current time.'

        enddo

        call cpu_time(t1)
        write(*,*) 'DNS BC Update time: ', t1-t0, ' seconds'
      endif

      call MPI_Bcast( Ua, 3*My*Mz, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

    end subroutine updateUalpha

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Third-order polynomial interpolation to "time"
    subroutine interpolateVelocity(time)
      implicit none

      real(dp), intent(in) :: time
      real(dp) :: c1, c2, c3, dt1, dt2, dt3


      !do kk=1,Mz
      !  do jj=1,My

      !    dt1 = (timeArray(jj,kk) - time)/dtArray(jj,kk)
      !    dt2 = (timeArray(jj,kk) - dtArray(jj,kk) - time)/dtArray(jj,kk)
      !    dt3 = (timeArray(jj,kk) - 2._dp*dtArray(jj,kk) - time)/dtArray(jj,kk)

      !    c1 = dt2*dt3/( (dt3-dt1)*(dt1-dt2) )
      !    c2 = dt1*dt3/( (dt2-dt3)*(dt1-dt2) )
      !    c3 = dt1*dt2/( (dt2-dt3)*(dt3-dt1) )

      !    Ua(jj,kk,:) = -(Uat1(jj,kk,:)*c1 + Uat2(jj,kk,:)*c2 + Uat3(jj,kk,:)*c3)

      !  enddo
      !enddo

          dt1 = (currentTime - time)/dx
          dt2 = (currentTime - dx - time)/dx
          dt3 = (currentTime - 2._dp*dx - time)/dx

          c1 = dt2*dt3/( (dt3-dt1)*(dt1-dt2) )
          c2 = dt1*dt3/( (dt2-dt3)*(dt1-dt2) )
          c3 = dt1*dt2/( (dt2-dt3)*(dt3-dt1) )

          Ua = -(Uat1*c1 + Uat2*c2 + Uat3*c3)


    end subroutine interpolateVelocity


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !Update the random field (called after each time step)
    subroutine updateRandomField()
      implicit none

      !logical, dimension(-Ny+1:,-Nz+1:), intent(in) :: update
      !integer, intent(out)                :: nUpdated

      real(dp) :: sqrt3=sqrt(3._dp)
      integer  :: jj, kk

      !Update random numbers in last field
      !nUpdated=0
      !do kk=-Nz+1,Nz+Mz
      !  do jj=-Ny+1,Ny+My
      !    if (update(jj,kk)) then
      !      Rx(-Nx:Nx+Nextra-1,jj,kk) = Rx(-Nx+1:Nx+Nextra,jj,kk)
      !      Ry(-Nx:Nx+Nextra-1,jj,kk) = Ry(-Nx+1:Nx+Nextra,jj,kk)
      !      Rz(-Nx:Nx+Nextra-1,jj,kk) = Rz(-Nx+1:Nx+Nextra,jj,kk)

      !      Rx(Nx+Nextra,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
      !      Ry(Nx+Nextra,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
      !      Rz(Nx+Nextra,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)

      !      nUpdated=nUpdated+1
      !    endif

      !  enddo
      !enddo

      Rx(-Nx:Nx-1,:,:) = Rx(-Nx+1:Nx,:,:)
      Ry(-Nx:Nx-1,:,:) = Ry(-Nx+1:Nx,:,:)
      Rz(-Nx:Nx-1,:,:) = Rz(-Nx+1:Nx,:,:)

      do kk=-Nz+1,Nz+Mz
        do jj=-Ny+1,Ny+My
            Rx(Nx,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
            Ry(Nx,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
            Rz(Nx,jj,kk) =sqrt3*(2._dp*rand(0) - 1._dp)
        enddo
      enddo


    end subroutine updateRandomField


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Update random field for periodic boundary conditions
    subroutine periodic(R)
      implicit none

      real(dp), dimension(-Nx:, -Ny+1:, -Nz+1:), intent(inout) :: R

      !Periodic in y
      if (pY) then
        R(:, My:My+Ny,:) = R(:, 1:Ny, :)
        R(:, 1-Ny:0,:)    = R(:, My-Ny:My, :)
      endif

      !Periodic in z
      if (pZ) then
        R(: ,:, Mz:Mz+Nz) = R(: , :, 1:Nz)
        R(: ,:, 1-Nz:0)    = R(: , :, Mz-Nz:Mz)
      endif


    end subroutine periodic

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Support functions for the DNS filter
    function btildek(k,n)
      integer,  intent(in) :: k
      real(dp), intent(in) :: n

      real(dp) :: btildek

      btildek = exp( -pi*0.5_dp*real(k*k,dp)/(n*n) )

    end function btildek

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    function bkfun(k,n,NN)
      integer, intent(in) :: k, NN
      real(dp), intent(in) :: n
      real(dp) :: bkfun

      real(dp) :: bksum2
      integer  :: i

      bksum2=0._dp
      do i=-NN,NN
        bksum2 = bksum2 + btildek(i, n)**2
      enddo

      bkfun = btildek(k,n)/sqrt( bksum2 )

    end function bkfun


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Support functions for the DNS filter
    function btildek2(dx,k,L)
      integer,  intent(in) :: k
      real(dp), intent(in) :: dx, L

      real(dp) :: btildek2

      btildek2 = exp( -pi*0.5_dp*(dx*dx*real(k*k,dp))/(L*L) )

    end function btildek2

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    ! Support functions for the DNS filter
    function btildek3(r,L)
      real(dp), intent(in) :: r, L

      real(dp) :: btildek3

      btildek3 = exp( -pi*0.5_dp*(r*r)/(L*L) )

    end function btildek3

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    function bk2(dx,k,L,NN)
      integer, intent(in) :: k, NN
      real(dp), intent(in) :: dx, L
      real(dp) :: bk2

      real(dp) :: bksum2
      integer  :: i


      bksum2=0._dp
      do i=-NN,NN
        bksum2 = bksum2 + btildek2(dx, i, L)**2
      enddo

      bk2 = btildek2(dx, k, L)/sqrt( bksum2 )

    end function bk2

end



