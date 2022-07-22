program molcas_ci_reader
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

   ! csfvec='22ud00'
   filename = 'NEON_scfvspt2vsfullci/molcas.log'

   call get_civecs_in_csfs(trim(filename), .false., csfs, civs)
   call csf_to_slater_basis_conversion(csfs, civs, sds, civecinsd, 0)
   call format_sd(sds, 1, finalsds)
   ! write(*,*) finalsds

contains

  subroutine get_civs_and_confs(logfile, caspt2, finalsds, civecinsd)
    implicit none
    character(len=:), intent(IN)  :: logfile
    logical, intent(IN) :: caspt2
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

   ! csfvec='22ud00'

   call get_civecs_in_csfs(logfile, caspt2, csfs, civs)
   call csf_to_slater_basis_conversion(csfs, civs, sds, civecinsd, 0)
   call format_sd(sds, 1, finalsds)

   end subroutine


   subroutine csf_to_slater_basis_conversion(csfs, civs, sds, civecinsd, S)
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
      INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
      integer(kind=ikind), intent(in) :: S
      character(len=*), intent(in) :: csfs(:)
      character(len=len(csfs(1))), allocatable :: sdvec(:)
      character(len=:), intent(out), allocatable:: sds(:)
      real(kind=dp), intent(in) ::civs(:, :)
      real(kind=dp), allocatable :: coeffs(:)
      logical :: found

      real(kind=dp), intent(out), allocatable ::civecinsd(:, :)
      integer(kind=ikind) :: no_roots, no_sds, no_u, no_d, no_mo, i, j, k, no_csfs, running_total
      character(len=len(csfs(1))) :: csfvec

      csfvec = csfs(1)
      write (*, *) csfvec
      no_mo = len(csfvec)
      no_roots = size(civs, dim=1)
      no_u = 0
      no_d = 0

      do i = 1, no_mo
      if (csfvec(i:i) .eq. 'u') then
         no_u = no_u + 1
      elseif (csfvec(i:i) .eq. 'd') then
         no_d = no_d + 1
      elseif (csfvec(i:i) .eq. '2') then
         no_u = no_u + 1
         no_d = no_d + 1
      end if
      end do
      write (*, *) 'U', no_u
      write (*, *) 'D', no_d

      no_sds = choose(no_mo, no_u)*choose(no_mo, no_d)
      no_csfs = size(civs, dim=2)

      allocate (character(len=no_mo) :: sds(no_sds))
      allocate (civecinsd(no_roots, no_sds))

      civecinsd = 0.0_dp

      running_total = 1
      do i = 1, no_csfs
         call csf2sd(csfs(i), S, coeffs, sdvec, no_mo)

         do j = 1, size(coeffs)
         if (coeffs(j) .ne. 0.0_dp) then

            found = .false.
            do k = 1, running_total
               if (sdvec(j) .eq. sds(k)) then
                  found = .true.
                  civecinsd(:, k) = civecinsd(:, k) + coeffs(j)*civs(:, i)
               end if
            end do
            if (.not. found) then
               sds(running_total) = sdvec(j)
               ! write(*,*) 'Putting sdvec ', sdvec(j),' to index ', running_total
               civecinsd(:, running_total) = coeffs(j)*civs(:, i)
               if (running_total .lt. no_sds) then
                  running_total = running_total + 1
               end if
            end if
         end if
         end do
         deallocate (coeffs, sdvec)
      end do

      write (*, *) 'Converted to ', no_sds, ' determinants'
   end subroutine

   !make into elemental procedure?
   subroutine format_sd(sds, no_closed, finalsds)
      character(len=*), intent(in) :: sds(:)
      integer(kind=ikind), intent(in) :: no_closed
      character(len=:), allocatable, intent(out):: finalsds(:)

      integer :: total_number, i

      total_number = len(sds(1)) + no_closed

      allocate (character(len=total_number) :: finalsds(size(sds)))

      do i = 1, size(sds)
         finalsds(i) (1:no_closed) = REPEAT('2', no_closed)
         finalsds(i) (no_closed:total_number) = sds(i)
      end do
   end subroutine

   subroutine get_number_of_active_orbitals(filename, no_active_orbs)
      character(len=*) :: filename
      character(len=60) :: str
      integer :: no_active_orbs
      integer :: io

      open (1, file=trim(filename))

      do
         read (1, '(A)', iostat=io) str
         if (index(str, 'RAS2 orbitals') .gt. 0) then
            read (str(34:len(str)), *, iostat=io) no_active_orbs
         end if
         if (io .lt. 0) then
            return
         end if
      end do
      close (1)
   end subroutine

   subroutine get_civecs_in_csfs(filename, caspt2, csfs, civs)
      implicit none

      integer :: no_active
      character(len=*) :: filename
      logical :: caspt2
      character(len=100) :: str, tempstr
      character(:), allocatable :: csfs(:)
      real(kind=dp), allocatable :: civs(:, :)
      integer :: io, no_closed, S, no_csfs, no_roots, root, idx, q, j, i

      open (1, file=trim(filename))

      ! this is a bunch of similar things which look for things in the output file.

      ! They loop through the file, and when they find something (using index), then they read the correct thing
      do
         read (1, '(A)', iostat=io) str
         if (io .lt. 0) then
            write (*, *) 'Stopping'
            stop
         end if
         if (index(str, 'Number of inactive orbitals') .gt. 0) then
            read (str(34:len(str)), *, iostat=io) no_closed
            exit
         end if
      end do
      do
         read (1, '(A)', iostat=io) str
         if (io .lt. 0) then
            write (*, *) 'Stopping'
            stop
         end if
         if (index(str, 'Number of active orbitals') .gt. 0) then
            read (str(34:len(str)), *, iostat=io) no_active
            exit
         end if
      end do
      do
         read (1, '(A)', iostat=io) str
         if (io .lt. 0) then
            write (*, *) 'Stopping'
            stop
         end if
         if (index(str, 'Spin quantum number') .gt. 0) then
            read (str(34:len(str)), *, iostat=io) S
            exit
         end if
      end do
      do
         read (1, '(A)', iostat=io) str
         if (io .lt. 0) then
            write (*, *) 'Stopping'
            stop
         end if
         if (index(str, 'Number of CSFs') .gt. 0) then
            read (str(34:len(str)), *, iostat=io) no_csfs
            exit
         end if
      end do

      allocate (character(no_active) :: csfs(no_csfs))

      if (caspt2) then

         do
            read (1, '(A)', iostat=io) str
            if (io .lt. 0) then
               return
            end if
            if (index(str, 'Number of CI roots used') .gt. 0) then
               read (str(34:len(str)), *, iostat=io) no_roots
               allocate (civs(no_roots, no_csfs))
               exit
            end if
         end do

         do root = 1, no_roots
         do
            read (1, '(A)', iostat=io) str
            if (io .lt. 0) then
               exit
            end if
            if (index(str, 'The CI coefficients for the MIXED state nr.') .gt. 0) then ! find this string
               do q = 1, 6 ! skip six lines
                  read (1, '(A)', iostat=io) str
               end do
               do j = 1, no_csfs
                  read (1, '(A)', iostat=io) str
                  read (str(1:15), *, iostat=io) idx ! first bit is index - can probably skip this and include a idx=idx+1
                  tempstr = trim(str(33:len(str))) ! second bit is the csf, the civ, and the coefficients
                  read (tempstr, *) csfs(idx), civs(root, idx)

               end do
               exit
            end if
         end do
         end do

      else
         do
            read (1, '(A)', iostat=io) str
            if (io .lt. 0) then
               return
            end if
            if (index(str, 'Number of root(s) required') .gt. 0) then
               read (str(34:len(str)), *, iostat=io) no_roots
               allocate (civs(no_roots, no_csfs))
               exit
            end if
         end do

         do root = 1, no_roots
         do
            read (1, '(A)', iostat=io) str
            if (io .lt. 0) then
               exit
            end if
            if (index(str, 'printout of CI-coefficients larger than') .gt. 0) then
               do q = 1, 2 ! skip two lines
                  read (1, '(A)', iostat=io) str
               end do
               do j = 1, no_csfs
                  read (1, '(A)', iostat=io) str
                  read (str, *, iostat=io) idx
                  read (str, *, iostat=io) idx, csfs(idx), civs(root, idx)  ! read twice to ensure that correct index is used.

               end do
               exit
            end if
         end do
         end do
      end if

      close (1)
      write (*, *) csfs(1)
      write (*, *) no_csfs, 'configurations read'
   end subroutine

   subroutine csf2sd(csfvec, m_s, coeffs, sdvec, no_mo)
      implicit none
      INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
      INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
      integer(kind=ikind), intent(in) :: m_s, no_mo
      character(len=no_mo), intent(IN) :: csfvec
      character(len=no_mo), allocatable, intent(out) :: sdvec(:)
      real(kind=dp), allocatable, intent(out) :: coeffs(:)
      integer(kind=ikind), allocatable :: comb(:)

      integer(kind=ikind) :: no_somo, no_domo, no_e, no_alpha, a(no_mo), b(no_mo), c(no_mo), i, &
                             no_a, i_s, i_a, i_b, idx, no_combs, i_domo, i_somo, i_alpha, i_beta, j
      integer(kind=ikind), allocatable :: combs(:, :)
      character :: curr_orb

      curr_orb = csfvec(1:1)
      no_somo = 0
      no_domo = 0
      a = 0
      b = 0
      c = 0

      do i = 1, no_mo
      if (csfvec(i:i) .eq. 'u') then
         no_somo = no_somo + 1
      elseif (csfvec(i:i) .eq. 'd') then
         no_somo = no_somo + 1
      elseif (csfvec(i:i) .eq. '2') then
         no_domo = no_domo + 1
      end if
      end do

      if (curr_orb .eq. '0') then
         c(1) = 1
      elseif (curr_orb .eq. 'u') then
         b(1) = 1
      elseif (curr_orb .eq. 'd') then
         a(1) = 1
         b(1) = -1
         c(1) = 1
      elseif (curr_orb .eq. '2') then
         a(1) = 1
      else
         write (*, *) csfvec
         write (*, *) 'ERROR - not a CSF'
         stop
      end if

      do i = 2, no_mo

         a(i) = a(i - 1)
         b(i) = b(i - 1)
         c(i) = c(i - 1)

         curr_orb = csfvec(i:i)

         if (curr_orb .eq. '0') then
            c(i) = c(i) + 1
         elseif (curr_orb .eq. 'u') then
            b(i) = b(i) + 1
         elseif (curr_orb .eq. 'd') then
            a(i) = a(i) + 1
            b(i) = b(i) - 1
            c(i) = c(i) + 1
         elseif (curr_orb .eq. '2') then
            a(i) = a(i) + 1
         else
            write (*, *) csfvec
            write (*, *) 'ERROR - not a CSF'
            stop
         end if

      end do

      no_e = 2*a(no_mo) + b(no_mo)
      ! spin = b(no_mo)
      no_alpha = (no_somo + m_s*2)/2
      no_combs = choose(no_somo, no_alpha)
      allocate (combs(no_alpha, no_combs), comb(no_alpha), sdvec(no_combs), coeffs(no_combs))
      idx = 1
      call gen(1, no_alpha, no_somo, idx, combs, comb)

      do j = 1, no_combs

         comb = combs(:, j)
         coeffs(j) = 1.0_dp
         i_domo = 1
         i_somo = 1
         i_alpha = 0
         i_beta = 0

         do i = 1, no_mo
            curr_orb = csfvec(i:i)
            if (curr_orb .eq. '0') then
               sdvec(j) (i:i) = '0'
            elseif (curr_orb .eq. 'u') then
               if (any(comb == i_somo)) then
                  sdvec(j) (i:i) = 'a'
                  coeffs(j) = coeffs(j)*(a(i) + b(i) - i_beta)
                  i_alpha = i_alpha + 1
               else
                  sdvec(j) (i:i) = 'b'
                  coeffs(j) = coeffs(j)*(a(i) + b(i) - i_alpha)
                  i_beta = i_beta + 1
               end if
               coeffs(j) = coeffs(j)/b(i)
               i_somo = i_somo + 1
            elseif (curr_orb .eq. 'd') then
               if (any(comb == i_somo)) then
                  sdvec(j) (i:i) = 'a'
                  coeffs(j) = coeffs(j)*(i_beta - a(i) + 1)
                  i_alpha = i_alpha + 1

                  if (iseven(b(i))) then
                     coeffs(j) = coeffs(j)*(-1)
                  end if
               else
                  sdvec(j) (i:i) = 'b'
                  coeffs(j) = coeffs(j)*(i_alpha - a(i) + 1)
                  i_beta = i_beta + 1
                  if (.not. iseven(b(i))) then
                     coeffs(j) = coeffs(j)*(-1)
                  end if
               end if
               coeffs(j) = coeffs(j)/(b(i) + 2)
               i_somo = i_somo + 1
            else if (curr_orb .eq. '2') then
               sdvec(j) (i:i) = '2'
               if (.not. iseven(b(i))) then
                  coeffs(j) = coeffs(j)*(-1)
               end if
               i_alpha = i_alpha + 1
               i_beta = i_beta + 1
            else
               write (*, *) 'error in step vector'
               stop
            end if
         end do
      end do

      ! sqrt the coeff as this is an amplitude
      do i = 1, size(coeffs)
         coeffs(i) = SIGN(sqrt(abs(coeffs(i))), coeffs(i))
      end do

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
   recursive subroutine gen(m, m_max, n_max, idx, output, comb)

      implicit none
      integer(kind=ikind), intent(in) :: m, m_max, n_max
      integer(kind=ikind) :: idx
      integer(kind=ikind) :: n
      integer(kind=ikind) :: output(:, :)
      integer(kind=ikind), dimension(m_max) :: comb

      if (m > m_max) then

         output(:, idx) = comb(:)
         idx = idx + 1
      else
         do n = 1, n_max
            if ((m == 1) .or. (n > comb(m - 1))) then
               comb(m) = n
               call gen(m + 1, m_max, n_max, idx, output, comb)
            end if
         end do
      end if

   end subroutine gen
   logical function iseven(a)
      integer, intent(in) :: a
      integer :: b
      b = mod(a, 2)

      if (b .eq. 0) then
         iseven = .true.
      elseif (b .eq. 1) then
         iseven = .false.
      else
         stop
      end if
   end function
end program
