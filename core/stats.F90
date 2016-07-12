module statistics
  use kinds, only : DP

  implicit none

  real(DP), allocatable :: stat(:,:,:)
  real(DP) :: start_time
  integer :: nstat, nrec

  contains

  ! allocate stat structures
  subroutine setup_stat(num)
    use size_m, only : nx1, ny1, nz1, nelv
    use tstep, only : time
    integer, intent(in) :: num

    nstat = num
    allocate(stat(nx1*ny1*nz1, nstat, nelv))
    stat = 0._dp
    start_time = time
    nrec = 0

  end subroutine setup_stat

  ! average statistics on the fly
  subroutine avg_stat_all(vx, vy, vz, p, t)
    use size_m, only : nx1, ny1, nz1, nelv
    real(DP), intent(in) :: vx(nx1*ny1*nz1, nelv)
    real(DP), intent(in) :: vy(nx1*ny1*nz1, nelv)
    real(DP), intent(in) :: vz(nx1*ny1*nz1, nelv)
    real(DP), intent(in) ::  p(nx1*ny1*nz1, nelv)
    real(DP), intent(in) ::  t(nx1*ny1*nz1, nelv)

    integer :: nxyz, i, j, k

    nxyz = nx1 * ny1 * nz1
   
    nrec = nrec + 1 
    do i = 1, nelv
      stat(:, 1,i) = stat(:, 1,i) + vx(:,i)
      stat(:, 2,i) = stat(:, 2,i) + vy(:,i)
      stat(:, 3,i) = stat(:, 3,i) + vz(:,i)
      stat(:, 4,i) = stat(:, 4,i) +  p(:,i)

      stat(:, 5,i) = stat(:, 5,i) + vx(:,i)*vx(:,i)
      stat(:, 6,i) = stat(:, 6,i) + vy(:,i)*vy(:,i)
      stat(:, 7,i) = stat(:, 7,i) + vz(:,i)*vz(:,i)
      stat(:, 8,i) = stat(:, 8,i) +  p(:,i)* p(:,i)

      stat(:, 9,i) = stat(:, 9,i) + vx(:,i)*vy(:,i)
      stat(:,10,i) = stat(:,10,i) + vy(:,i)*vz(:,i)
      stat(:,11,i) = stat(:,11,i) + vz(:,i)*vx(:,i)
    enddo

  end subroutine avg_stat_all

  subroutine output_stat(filename)
    use kinds, only : strlen
    use size_m, only : nx1, ny1, nz1
    use input, only : param
    use mesh, only : shape_x, start_x, end_x
    use tstep, only : time, dt
    
    character(strlen), intent(in) :: filename 

    character(strlen) val1, val2, val3, val4, val5, val6
    character(strlen) val7, val8, val9, val10

    real(DP), allocatable :: stat_xy(:,:), w1(:), w2(:)
    integer :: i

    allocate(stat_xy(nx1*ny1, nstat))
    allocate(w1(nx1*ny1*shape_x(1)*shape_x(2)))
    allocate(w2(nx1*ny1*shape_x(1)*shape_x(2)))

    do i = 1, nstat
      call z_average(stat_xy(1, i), stat(:,i,:), w1, w2)
    enddo

    open(unit=33,form='unformatted',file=filename)

    write(val1,'(1p15e17.9)') 1/param(2)           ! Reynolds number    
    write(val2,'(1p15e17.9)') end_x - start_x      ! domain size
    write(val3,'(9i9)') shape_x                    ! number of elements 
    write(val4,'(9i9)') nx1-1,ny1-1,nz1-1                ! polynomial order
    write(val5,'(9i9)')       nstat                     ! number of saved statistics 
    write(val6,'(1p15e17.9)') start_time                 ! start time
    write(val7,'(1p15e17.9)') time                       ! end time
    write(val8,'(1p15e17.9)') nrec*DT                    ! average time
    write(val9,'(1p15e17.9)') DT                         ! time step
    write(val10,'(9i9)')      nrec                       ! number of time records

    write(33) '(Re ='//trim(val1) &
      //') (Lx, Ly, Lz ='//trim(val2) &
      //') (nelx, nely, nelz ='//trim(val3) &
      //') (Polynomial order ='//trim(val4) &
      //') (Nstat ='//trim(val5) &
      //') (start time ='//trim(val6) &
      //') (end time ='//trim(val7) &
      //') (average time ='//trim(val8) &
      //') (time step ='//trim(val9) &
      //') (nrec ='//trim(val10) &
      //')'

    write(33) 1/param(2), &
      end_x - start_x, &
      shape_x, &
      nx1-1   , ny1-1   , nz1-1, &
      nstat, &
      start_time, &
      time, &
      nrec*DT, &
      DT, &
      nrec 

    do i = 1, nstat
      write(33) stat_xy(:,i)
    enddo

    close(33) 

    nrec = 0
    stat = 0._dp
    start_time = time

  end subroutine output_stat

end module statistics
