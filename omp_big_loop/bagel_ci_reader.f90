program bagel_ci_reader
   implicit none
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
   INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
   character(len=:), allocatable :: csfs(:)
   real(kind=dp), allocatable :: civs(:, :)
   real(kind=dp), allocatable :: coeffs(:)
   character(len=:), allocatable:: sds(:)
   character(len=:), allocatable:: finalsds(:)
   real(kind=dp), allocatable ::civecinsd(:, :)
   integer :: no_active_orbs
   character(len=100) :: filename

   filename = 'out.out'

   call get_civs_and_confs(filename, .false., finalsds, civecinsd)

contains

  subroutine get_civs_and_confs(logfile, caspt2, finalsds, civecinsd)

    ! main function
    implicit none
    character(len=*), intent(IN)  :: logfile
    logical, intent(IN) :: caspt2
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
   INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
   character(len=:), allocatable :: csfs(:)
   real(kind=dp), allocatable :: civs(:, :)
   real(kind=dp), allocatable :: coeffs(:)
   character(len=:), allocatable:: sds(:)
   character(len=:), allocatable, intent(out):: finalsds(:)
   character(len=100) :: fmtstring
   real(kind=dp), allocatable, intent(out) ::civecinsd(:, :)
   integer :: no_closed,i
   character(len=100) :: filename


   ! get the csfs  from the output file
   call get_civecs_in_sds(logfile, caspt2, csfs, civs, no_closed)

   ! format the output so that the it includes the number of closed
   call format_sd(csfs, no_closed, finalsds)

   fmtstring = '(A' // repeat(',2X, f12.9', size(civs(:,1))) // ')'

   do i = 1, size(civs(1,:))
   write(*,fmtstring ) finalsds(i), civs(:,i) ! this only writes four st
   end do
   end subroutine



   !make into elemental procedure?
   subroutine format_sd(sds, no_closed, finalsds)
    character(len=*), intent(in) :: sds(:)
    integer(kind=ikind), intent(in) :: no_closed
    character(len=:), allocatable, intent(out):: finalsds(:)
    character :: a

    integer :: total_number, i, j

    total_number = (len(sds(1)) + no_closed+1)*4
    

    allocate (character(len=total_number) :: finalsds(size(sds)))

    do i = 1, size(sds)
    finalsds(i) (1:4*no_closed) = REPEAT('1 2 ', no_closed)
    do j = 1, len(sds(1))
    a =sds(i)(j:j) 
    if (a .eq. '2') then
      finalsds(i)(4*(no_closed+j):4*(no_closed+j+1))= '1 2 '
    else if (a .eq. 'a') then
      finalsds(i)(4*(no_closed+j):4*(no_closed+j+1)) = '1 0 ' 
    else if (a .eq. 'b') then
      finalsds(i)(4*(no_closed+j):4*(no_closed+j+1)) = '0 2 ' 
    else if (a .eq. '.') then
      finalsds(i)(4*(no_closed+j):4*(no_closed+j+1)) = '0 0 ' 
    endif 
    end do
    end do
    
   end subroutine


   subroutine get_civecs_in_sds(filename, caspt2, sds, civs, no_closed)
    implicit none

     integer, intent(out) :: no_closed
    character(len=*), intent(in) :: filename
    logical, intent(in) :: caspt2
    logical :: cart = .false., found = .false.
    character(len=100) :: str, tempstr
    character(:), allocatable, intent(out) :: sds(:)
    real(kind=dp), allocatable, intent(out) :: civs(:, :)
    integer :: io, no_active, S, no_sds, no_roots, root, idx, q, j, i, no_electrons, no_a, no_b, k, running_total
    real(kind=dp) :: coeff
    real(kind=dp), allocatable :: rotmat(:,:)
    character(:), allocatable :: test_sd

    open (1, file=trim(filename))

    ! this is a bunch of similar things which look for things in the output file.

    ! They loop through the file, and when they find something (using index), then they read the correct thing
    do
    read (1, '(A)', iostat=io) str
    if (io .lt. 0) then
    write (*, *) 'Stopping'
    stop
  end if
  if (index(str, 'Cartesian basis functions are used') .gt. 0) then
  cart = .true.
  exit
end if
  end do
  if (.not. cart ) then
  write(*,*) 'Looks like this might be a calculation with spherical basis functions - this will fail!'
end if

  rewind(1) ! go back to beginning of file
  do
  read (1, '(A)', iostat=io) str
  if (io .lt. 0) then
  write (*, *) 'Stopping'
  stop
end if
  if (index(str, 'Number of electrons') .gt. 0) then
  read (str(30:35), *, iostat=io) no_electrons
  exit
end if
  end do

  do
  read (1, '(A)', iostat=io) str
  if (io .lt. 0) then
  write (*, *) 'Stopping'
  stop
end if
  if (index(str, 'nstate') .gt. 0) then
  read (str(19:25), *, iostat=io) no_roots
  exit
end if
  end do

  do
  read (1, '(A)', iostat=io) str
  if (io .lt. 0) then
  write (*, *) 'Stopping'
  stop
end if
  if (index(str, 'nclosed') .gt. 0) then
  read (str(19:25), *, iostat=io) no_closed
  exit
end if
  end do

  do
  read (1, '(A)', iostat=io) str
  if (io .lt. 0) then
  write (*, *) 'Stopping'
  stop
end if
  if (index(str, 'nact') .gt. 0) then
  read (str(19:25), *, iostat=io) no_active
  exit
end if
  end do

  ! need to find out the spin!

  allocate(character(no_active) :: test_sd)

  no_a = 0
  no_b = 0

  do
  read (1, '(A)', iostat=io) str
  if (io .lt. 0) then
  write (*, *) 'Stopping'
  stop
end if
  if (index(str, 'ci vector, state') .gt. 0) then
  read (1, '(A)', iostat=io) str
  read (str(9:9+no_active), '(A)', iostat=io) test_sd
  test_sd = str(8:8+no_active)
  do i = 1, len(test_sd)
  if (test_sd(i:i) .eq. '2') then
  no_a = no_a + 1
  no_b = no_b + 1
elseif (test_sd(i:i) .eq. 'a') then
  no_a = no_a + 1
elseif (test_sd(i:i) .eq. 'b') then
  no_b = no_b + 1
endif
  end do
  exit
endif
  enddo


  deallocate(test_sd)

  S = no_a - no_b

  no_sds = choose(no_active,no_a)*choose(no_active,no_b)


  allocate(character(no_active) :: sds(no_sds))
  allocate(civs(no_roots, no_sds))

  running_total=0
  rewind(1)

  do root = 1, no_roots
  do
  read (1, '(A)', iostat=io) str
  if (io .lt. 0) then
  write(*,*) 'stopping'
  exit
end if
  if (index(str, 'ci vector, state  ') .gt. 0) then
  do j = 1, no_sds
  found=.false.
  read (1, *,iostat=io) str, coeff
  test_sd = str(1:1+no_active)
  ! read(str(9+no_active:len(str)), *) coeff

  do k = 1, running_total
  if (test_sd .eq. sds(k)) then
  found = .true.
  civs(root, k) = coeff
endif
  end do
  if (.not. found) then
  running_total=running_total+1
  sds(k) = test_sd
  civs(root, k) = coeff
endif
  enddo

  exit
endif
  end do
  end do



  if (caspt2) then
  allocate(rotmat(no_roots,no_roots))
  do
  read (1, '(A)', iostat=io) str
  if (io .lt. 0) then
  write(*,*) 'stopping'
  exit
end if
  if (index(str, 'XMS-CASPT2 rotation matrix') .gt. 0) then
  do root = 1,no_roots
  read(1, *) rotmat(root,:)
  enddo
  exit
endif
  enddo


    do i = 1, no_sds
  civs(:,i) = matmul(rotmat,civs(:,i))
  enddo
  deallocate(rotmat)
    endif


  close (1)
   end subroutine
   function choose(n, k, err)
  integer :: choose
  integer, intent(in) :: n, k
  integer, optional, intent(out) :: err

  integer :: imax, i, imin, ie

  ie = 0
  if ((n < 0) .or. (k < 0)) then
  write (*, *) "negative in choose"
  choose = 0
  ie = 1
else
  if (n < k) then
  choose = 0
else if (n == k) then
  choose = 1
else
  imax = max(k, n - k)
  imin = min(k, n - k)
  choose = 1
  do i = imax + 1, n
  choose = choose*i
  end do
  do i = 2, imin
  choose = choose/i
  end do
end if
end if
  if (present(err)) err = ie
   end function choose

end program
