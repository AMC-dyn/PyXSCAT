
   module onerdm
       use twordms
       contains

        function excitations(det1,det2,Nint) result(n_excitations)
           implicit none
           integer*8, intent(in) :: det1(Nint,2), det2(Nint,2)
           integer , intent(in) :: Nint
           integer :: n_excitations

           integer :: l

           n_excitations = &
                   popcnt(xor( det1(1,1), det2(1,1)) ) + &
                           popcnt(xor( det1(1,2), det2(1,2)) )

           do l=2,Nint
               n_excitations = n_excitations + &
                       popcnt(xor( det1(l,1), det2(l,1)) ) + &
                       popcnt(xor( det1(l,2), det2(l,2)) )
           end do
           n_excitations = ishft(n_excitations,-1)

       end function

     subroutine onerdm_creat(confs,civs,onerdm,maxnmo,state1,state2)
            INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
         integer(kind=ikind), intent(in), dimension(:,:) :: confs
         real(kind=dp), intent(in), dimension(:,:) :: civs
         real (kind=dp), intent(out), dimension(:,:), allocatable :: onerdm
         integer(kind=ikind),dimension(:,:), allocatable :: ep2, ndiff
         integer(kind=ikind), intent(out) :: maxnmo
            integer(kind=ikind),intent(in):: state1,state2



         call maxcoincidence(confs, ep2,ndiff)

         call createonerdm(confs,civs,ndiff,ep2,onerdm,maxnmo,state1,state2)



         end subroutine


!     subroutine maxcoincidence(confs, ep2,ndiff)
!        use types
!        implicit none
!        integer(kind=ikind), intent(in), dimension(:,:) :: confs
!        integer(kind=ikind), intent(out), dimension(:,:), allocatable :: ep2, ndiff
!        integer(kind=ikind), dimension(size(confs(:,1)), size(confs(1,:))):: matdum
!        integer(kind=ikind), dimension(:,:), allocatable :: mat1
!        integer(kind=ikind) :: i,j, c1,c2, count
!
!
!        allocate(ep2(size(confs(:,1)),size(confs(:,1))))
!        ep2=1
!        count=0
!        matdum=0
!        print*,'holaa'
!            do i=1,size(confs(:,1))
!                count=0
!                do j=1,size(confs(1,:))
!
!
!                    if (confs(i,j)/=0) then
!                        count=count+1
!                        matdum(i,count)=j
!                        end if
!                end do
!            end do
!           print*, 'matdum constructed'
!            allocate(mat1(size(confs(:,1)),count), ndiff(size(confs(:,1)),size(confs(:,1))))
!            ndiff=0
!            mat1=matdum(:,1:count)
!            print*, 'mat1 constructed'
!            do c1=1, size(confs(:,1))
!                do c2=c1+1,size(confs(:,1))
!                    do i=1,size(mat1(1,:))
!                        if  (mat1(c1,i) /= mat1(c2,i)) then
!
!
!                            do j=1,size(mat1(1,:))
!                                if (mat1(c1,i) /= mat1(c2,j)) then
!                                    ep2(c1,c2)=-ep2(c1,c2)
!                                    ep2(c2,c1)=ep2(c1,c2)
!
!                                end if
!                            end do
!                        end if
!                    end do
!
!                    do j=1,size(confs(1,:))
!                        if (confs(c1,j)/=confs(c2,j)) then
!                           ndiff(c1,c2)= ndiff(c1,c2) +1
!                           ndiff(c2,c1)= ndiff(c2,c1) +1
!                        end if
!                    end do
!
!                end do
!            end do
!
!    print*,ndiff(1,2)
!    end subroutine maxcoincidence


    subroutine createonerdm(confs,civs,ndiff,ep2,onerdm,maxnmo,state1,state2)

    implicit none

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

    integer(kind=ikind), intent(in),dimension(:,:) :: confs
    integer(kind=ikind), intent(in), dimension(:,:) :: ep2, ndiff
    real(kind=dp), intent(in), dimension(:,:) :: civs
    integer(kind=ikind), intent(in) :: state1,state2


    real(kind=dp), dimension(:,:), intent(out),allocatable :: onerdm
    integer(kind=ikind), intent(out) :: maxnmo

    integer(kind=ikind) :: ep,nc1,lconfs,norbs,nc2,sorb,rorb,qorb,porb,p,q,r,s,c1,c2,count1,count

    integer(kind=ikind) :: i,i1,n,count2,eg, ndiff1


    integer(kind=ikind):: diffs1,diffs2
    integer(kind=ikind) :: spin1,spin2

    lconfs=size(confs(1,:))
    maxnmo=lconfs/2
    nc1=1
    nc2=1
    ep=1
    allocate(onerdm(maxnmo,maxnmo))
    onerdm=0.0_dp
    do c1=1,size(confs(:,1))
        do c2=1,size(confs(:,1))

            ndiff1=ndiff(c1,c2)
            ep=ep2(c1,c2)


            if (ndiff1==2) then

                print*,'cross terms',c1,c2
                do n=1,size(confs(c1,:))
                    if (confs(c1,n) /= confs(c2,n)) then
                        if (confs(c1,n) /= 0) THEN
                            diffs1=nint((n) / 2.0_dp + 0.1)
                            spin1=confs(c1,n)



                        elseif (confs(c2,n) /= 0) THEN
                            diffs2=nint((n) / 2.0_dp + 0.1)
                            spin2=confs(c2,n)
                        end if
                    end if
                enddo


                    sorb = diffs1
                    qorb = diffs2


                    eg = 1.00


                        onerdm(sorb , qorb ) = onerdm(sorb , qorb ) + &
                                civs(c1,state1) * civs(c2,state2) * ep

!                        print*,'conf1', pack(confs(c1,:),confs(c1,:)/=0)
!                        print*,'conf2',pack(confs(c2,:),confs(c2,:)/=0)
!                        print*,'phase', c1,c2,ep

                        print*, 'orbitals involved', sorb,qorb, ep

            elseif (ndiff1== 0) then

                ep = 1

                do i1=1,size(confs(c1,:))
                    if (confs(c1,i1)/=0) then

                        sorb=nint((i1) / 2.0_dp + 0.1)

                        onerdm( sorb , sorb)  = onerdm(sorb , sorb) + civs(c1,state1) * civs(c1,state2) * ep


                    end if
                    enddo

            end if


        end do

    end do
    count=0.0
    open(file='onerdm.dat', unit=15)
    do i=1,size(onerdm(:,1))

         write(15,'(1000F14.7)')( onerdm(i,c1) ,c1=1,size(onerdm(i,:)))
    end do
    close(15)

     open(file='configurations.dat', unit=15)
    do i=1,size(confs(:,1))

         write(15,*)( confs(i,c1) ,c1=1,size(confs(i,:)))
    end do
    close(15)

    print*, 'onerdm calculated'



    end subroutine


     subroutine one_rdm_slow(file_read,onerdm,maxnmo)

          Use, intrinsic :: iso_fortran_env, Only : iostat_end
    implicit none

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

    character(len=20),intent(in) :: file_read



    real(kind=dp), dimension(:,:), intent(out),allocatable :: onerdm
    integer(kind=ikind), intent(inout) :: maxnmo

    integer(kind=ikind) :: ep,nc1,lconfs,norbs,nc2,sorb,rorb,qorb,porb,n1a,n1b,error_1,error_2
          integer(kind=ikind),allocatable, dimension(:) ::Nalpha,Nbeta
          real(kind=dp), dimension(:), allocatable :: civs
          real(kind=dp):: civ1,civ2

    integer(kind=ikind) :: i,i1,n,count2,eg, ndiff1,n2a,n2b,j,count,c1,c2
    integer(kind=ikind), dimension(:), allocatable :: c1a,c1b,c2a,c2b





    nc1=1
    nc2=1
    ep=1
    allocate(onerdm(maxnmo,maxnmo))
    onerdm=0.0_dp
    open(15,file=file_read)

    count=0
    do while (error_1 /= -1)
        read(15,*, iostat = error_1) j
        count=count+1

    enddo
    print*,count
    close(15)
      allocate(civs(count-1), Nalpha(count-1), Nbeta(count-1))
          civs=0
          Nalpha=0
          Nbeta=0
        open(15,file=file_read)
        do i=1,count-1
            read(15,*)j, civs(i), Nalpha(i), Nbeta(i)
        enddo
        close(15)

         do c1=1,count-1

             do c2=1,count-1

                 N1a=Nalpha(c1)
                 N2a=Nalpha(c2)

                 N1b=Nbeta(c1)
                 N2b=Nbeta(c2)
                 civ1=civs(c1)
                 civ2=civs(c2)


            ep=1

            if (N1a==N2a .and. N1b==N2b) then

            !case where everything is the same
                ep=1
                c1a=integer2binary_orbs(N1a)

                c1b=integer2binary_orbs(N1b)

                do i=1,size(c1a)
                   sorb=c1a(i)

                   onerdm( sorb , sorb)  = onerdm(sorb , sorb) + civ1 * civ2
                end do

                 do i=1,size(c1b)
                   sorb=c1b(i)

                   onerdm( sorb , sorb)  = onerdm(sorb , sorb) + civ1 * civ2

                end do




            elseif (N1a==N2a) then
                ndiff1=0
                c1b=integer2binary_orbs(N1b)
                c2b=integer2binary_orbs(N2b)


                do i=1,size(c1b)
                    if (c1b(i)/=c2b(i)) then
                        sorb = c1b(i)
                        qorb = c2b(i)

                        do j=1,size(c2b)
                            if (c1b(i)/=c2b(j)) then
                                ep=-ep

                            end if
                        end do

                        onerdm(sorb , qorb ) = onerdm(sorb , qorb ) + &
                                civ1*civ2*ep



                    end if
                end do


            elseif (N1b==N2b) then
                     ndiff1=0
                c1a=integer2binary_orbs(N1a)
                c2a=integer2binary_orbs(N2a)


                do i=1,size(c1a)
                    if (c1a(i)/=c2a(i)) then
                        sorb = c1a(i)
                        qorb = c2a(i)

                        do j=1,size(c2a)
                            if (c1a(i)/=c2a(j)) then
                                ep=-ep

                            end if
                        end do

                        onerdm(sorb , qorb ) = onerdm(sorb , qorb ) + &
                                civ1*civ2*ep



                    end if
                end do
                end if

            end do

enddo



    open(file='onerdm.dat', unit=15)
    do i=1,size(onerdm(:,1))

         write(15,'(1000F14.7)')( onerdm(i,c1) ,c1=1,size(onerdm(i,:)))
    end do
    close(15)


    print*, 'onerdm calculated'



    end subroutine


       subroutine one_rdm_bit(file_read,onerdm, maxnmo,numberlines)

           Use, intrinsic :: iso_fortran_env, Only : iostat_end
           implicit none

           INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
           INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)


           integer(kind=ikind):: numberlines
           character(len=20),intent(in) :: file_read
           real(kind=dp), dimension(:), allocatable :: civs
           integer*8,dimension(:,:,:), allocatable :: Nalphbet
           integer*8:: j, i, error_1,count,temp,sp
           real(kind=dp),dimension(:,:), allocatable :: onerdm
           integer(kind=ikind),intent(inout):: maxnmo
           integer*8:: exc(0:2,2,2),deg,buffer,c1,c2,ep,n1a,n2a,n1b,n2b,low,high,number,red_vec(maxnmo*2)
           double precision :: phase, c,civ1,civ2
           integer*8:: hole, particle, tmp,tz,final,max,prueba(maxnmo*2)
           real(kind=dp)::time1,time2


           allocate(onerdm(maxnmo,maxnmo))
           onerdm=0.0_dp
           count=numberlines


           allocate(civs(count), Nalphbet(1,2,count))
           civs=0
           Nalphbet=0
           open(15,file=file_read)
           do i=1,count
               read(15,*)j, civs(i), Nalphbet(1,1,i), Nalphbet(1,2,i)
           enddo
           close(15)

           maxnmo=maxval(integer2binary_orbs(maxval(Nalphbet(1,1,:))))

          ! call get_excitation(Nalphbet(1,1,1),Nalphbet(1,1,3),exc,deg,phase,1)
!           print*,exc(0,1,1),exc(1,1,1), exc(1,2,1),deg
!           print*,integer2binary_orbs(Nalphbet(1,1,1))
!           print*,integer2binary_orbs(Nalphbet(1,1,3))
!
!           print*,integer2binary_orbs(Nalphbet(1,2,1))
!           print*,integer2binary_orbs(Nalphbet(1,2,3))
!
!           print*,'xor', xor( Nalphbet(1,1,1), Nalphbet(1,1,3))
!           print*,'iand',iand(xor( Nalphbet(1,1,1), Nalphbet(1,1,3)), Nalphbet(1,1,3))
!           print*,'iand',iand(xor( Nalphbet(1,2,1), Nalphbet(1,2,3)), Nalphbet(1,2,3))
!           print*,'iand',iand(xor( Nalphbet(1,1,1), Nalphbet(1,1,3)), Nalphbet(1,1,3))/=0_8
!           print*,'iand',iand(xor( Nalphbet(1,2,1), Nalphbet(1,2,3)), Nalphbet(1,2,3))/=0_8
!           print*,'trail particle',trailz(iand(xor( Nalphbet(1,1,1), Nalphbet(1,1,3)), Nalphbet(1,1,3)))
!           print*,'trail hole',trailz(iand(xor( Nalphbet(1,1,1), Nalphbet(1,1,3)), Nalphbet(1,1,1)))
!           print*,'conf', integer2binary_orbs_bit(Nalphbet(1,:,1))
           onerdm=0.0_dp
!            call cpu_time(time1)
!           red_vec=integer2binary_orbs_bit(Nalphbet(1,:,3),maxnmo*2)
!           number=size(pack(red_vec(51:53),red_vec(51:53)/=0))
!           call cpu_time(time2)
!           phase=(-1)**number
!           print*,'previous number', number, time2-time1
!             call cpu_time(time1)
!            call combine_alpha_beta(Nalphbet(1,1,3), Nalphbet(1,2,3), final,max)
!             red_vec=integer2binary_orbs_bit(Nalphbet(1,:,3),maxnmo*2)
!
!            buffer=final
!           prueba=0
!            do while (buffer/= 0_8)
!                tz=trailz(buffer)
!                prueba(tz+1)=1
!                buffer = iand(buffer,buffer-1_8)
!
!           end do
!           print*,'suma vectores',sum(prueba-red_vec)
!           number=popcnt(ibits(final,53,2))
!            call cpu_time(time2)
!           print*,'new number', number, time2-time1


            do c1=1,count



                    civ1=civs(c1)



                        do i=1,2
                            buffer = Nalphbet(1,i,c1)
                            do while (buffer /= 0_8)
                                j = trailz(buffer)

                                onerdm(j+1,j+1) = onerdm(j+1,j+1 ) + &
                                civ1*civ1

                                buffer = iand(buffer,buffer-1_8)
                            end do


                        end do
                    do c2=c1+1,count

                        !number=excitations(,Nalphbet(1,:,c2),1)
                         number= popcnt(xor( Nalphbet(1,1,c1), Nalphbet(1,1,c2)) ) + &
                           popcnt(xor( Nalphbet(1,2,c1), Nalphbet(1,2,c2)) )
                        number=ishft(number,-1)

                    if  (number== 1) then

                        if (Nalphbet(1,1,c1)==Nalphbet(1,1,c2)) then
                            tmp = xor( Nalphbet(1,2,c1), Nalphbet(1,2,c2))
                            particle = iand(tmp, Nalphbet(1,2,c2))
                            hole = iand(tmp, Nalphbet(1,2,c1))
                            if (particle /= 0_8) then
                                tz = trailz(particle)
                                j=tz+1
                            end if

                            if (hole /= 0_8) then
                                tz = trailz(hole)
                                i=tz+1
                            end if
                            low=i*2
                            high=j*2
                        else
                            tmp = xor( Nalphbet(1,1,c1), Nalphbet(1,1,c2))
                            particle = iand(tmp, Nalphbet(1,1,c2))
                            hole = iand(tmp, Nalphbet(1,1,c1))
                            if (particle /= 0_8) then
                                tz = trailz(particle)
                                j=tz+1
                            end if

                            if (hole /= 0_8) then
                                tz = trailz(hole)
                                i=tz+1
                            end if
                             low=i*2-1
                             high=j*2-1



                        end if
                        if (low>high) then
                          temp=low
                          low=high
                          high=temp
                      end if
                        phase=1
                        civ2=civs(c2)




                       call combine_alpha_beta(Nalphbet(1,1,c2), Nalphbet(1,2,c2), final)

                       number=popcnt(ibits(final,low,(high-low-1)))
                      ! print*,'newnumber',number
!                       red_vec=integer2binary_orbs_bit(Nalphbet(1,:,c2),maxnmo*2)
!                       number=size(pack(red_vec(low+1:high-1),red_vec(low+1:high-1)/=0))
!                        print*,'oldnumber',number


                       phase=(-1)**number
                       c = phase*civ1*civ2

                       onerdm(i,j) = onerdm(i,j) + c
                       onerdm(j,i) = onerdm(j,i) + c


                    end if
                    end do

                end do
            print*,'finished', 1_8
           open(file='onerdm_mine.dat', unit=15)
             do i=1,size(onerdm(:,1))

                 write(15,'(1000F14.7)')( onerdm(i,c1) ,c1=1,size(onerdm(i,:)))
            end do
!           call compute_density_matrix(Nalphbet,count-1,civs,maxnmo, &
!               1,onerdm)
!
!
!               open(file='onerdm_2.dat', unit=15)
!             do i=1,size(onerdm(:,1))
!
!                 write(15,'(1000F14.7)')( onerdm(i,c1) ,c1=1,size(onerdm(i,:)))
!            end do
!    close(15)

       end subroutine


       subroutine compute_density_matrix(det,Ndet,coef,mo_num, &
               Nint,density_matrix)
           implicit none
           INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
           integer*8, intent(in) :: det(Nint,2,Ndet)
           integer, intent(in) :: Ndet, Nint, mo_num
           double precision, intent(in) :: coef(Ndet)
           double precision, intent(out) :: density_matrix(mo_num,mo_num)
           integer :: i,j,k,l,ispin,ishift
           integer*8:: buffer
           integer :: deg
           integer :: exc(0:2,2,2)
           double precision :: phase, c
           integer :: n_excitations
           density_matrix = 0.d0



           do k=1,Ndet
               do ispin=1,2
                   ishift = 1
                   do i=1,Nint
                       buffer = det(i,ispin,k)
                       do while (buffer /= 0_8)
                           j = trailz(buffer) + ishift
                           density_matrix(j,j) = density_matrix(j,j) &
                                   + coef(k)*coef(k)
                           buffer = iand(buffer,buffer-1_8)
                       end do
                       ishift = ishift+64
                   end do
               end do

               do l=1,k-1
                   if (excitations(det(1,1,k),det(1,1,l),Nint) /= 1) then
                       cycle
                   end if
                   call get_excitation(det(1,1,k),det(1,1,l),exc,deg,phase,Nint)
                   if (exc(0,1,1) == 1) then
                       i = exc(1,1,1)
                       j = exc(1,2,1)


                   else
                       i = exc(1,1,2)
                       j = exc(1,2,2)

                   end if
                   c = phase*coef(k)*coef(l)

                   c = c+c
                   density_matrix(j,i) = density_matrix(j,i) + c
                   density_matrix(i,j) = density_matrix(i,j) + c
               end do
           end do
         !   print*,'diag1',density_matrix(26,26)
       end subroutine

       subroutine get_double_excitation(det1,det2,exc,phase,Nint)
           implicit none
           integer, intent(in) :: Nint
           integer*8, intent(in) :: det1(Nint,2), det2(Nint,2)
           integer, intent(out) :: exc(0:2,2,2)
           double precision, intent(out) :: phase
           integer :: l, ispin, idx_hole, idx_particle, ishift
           integer :: i,j,k,m,n,high, low,a,b,c,d,nperm,tz,nexc
           integer*8 :: hole, particle, tmp
           double precision, parameter :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)
           exc(0,1,1) = 0
           exc(0,2,1) = 0
           exc(0,1,2) = 0
           exc(0,2,2) = 0
           nexc=0
           nperm=0
           do ispin = 1,2
               idx_particle = 0
               idx_hole = 0
               ishift = -63
               do l=1,Nint
                   ishift = ishift + 64
                   if (det1(l,ispin) == det2(l,ispin)) then
                       cycle
                   end if
                   tmp = xor( det1(l,ispin), det2(l,ispin) )
                   particle = iand(tmp, det2(l,ispin))
                   hole = iand(tmp, det1(l,ispin))
                   do while (particle /= 0_8)
                       tz = trailz(particle)
                       nexc = nexc+1
                       idx_particle = idx_particle + 1
                       exc(0,2,ispin) = exc(0,2,ispin) + 1
                       exc(idx_particle,2,ispin) = tz+ishift
                       particle = iand(particle,particle-1_8)
                   end do
                   do while (hole /= 0_8)
                       tz = trailz(hole)
                       nexc = nexc+1
                       idx_hole = idx_hole + 1
                       exc(0,1,ispin) = exc(0,1,ispin) + 1
                       exc(idx_hole,1,ispin) = tz+ishift
                       hole = iand(hole,hole-1_8)
                   end do
                   if (nexc == 4) exit
               end do

               do i=1,exc(0,1,ispin)
                   low = min(exc(i,1,ispin),exc(i,2,ispin))
                   high = max(exc(i,1,ispin),exc(i,2,ispin))
                   j = ishft(low-1,-6)+1
                   n = iand(low,63)
                   k = ishft(high-1,-6)+1
                   m = iand(high,63)
                   if (j==k) then
                       nperm = nperm + popcnt(iand(det1(j,ispin), &
                               iand( ibset(0_8,m-1)-1_8, ibclr(-1_8,n)+1_8 ) ))
                   else
                       nperm = nperm + popcnt(iand(det1(k,ispin), &
                               ibset(0_8,m-1)-1_8)) &
                               + popcnt(iand(det1(j,ispin), &
                                       ibclr(-1_8,n) +1_8))
                       do l=j+1,k-1
                           nperm = nperm + popcnt(det1(l,ispin))
                       end do
                   end if
               end do
               if (exc(0,1,ispin) == 2) then
                   a = min(exc(1,1,ispin), exc(1,2,ispin))
                   b = max(exc(1,1,ispin), exc(1,2,ispin))
                   c = min(exc(2,1,ispin), exc(2,2,ispin))
                   d = max(exc(2,1,ispin), exc(2,2,ispin))
                   if (c>a .and. c<b .and. d>b) nperm = nperm + 1
                   exit
               end if
           end do
           phase = phase_dble(iand(nperm,1))

       end

       subroutine get_single_excitation(det1,det2,exc,phase,Nint)
           implicit none
           integer, intent(in) :: Nint
           integer*8, intent(in) :: det1(Nint,2)
           integer*8, intent(in) :: det2(Nint,2)
           integer, intent(out) :: exc(0:2,2,2)
           double precision, intent(out) :: phase
           integer :: tz, l, ispin, ishift, nperm, i, j, k, m, n, high, low
           integer*8 :: hole, particle, tmp
           double precision, parameter :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)

           exc(0,1,1) = 0
           exc(0,2,1) = 0
           exc(0,1,2) = 0
           exc(0,2,2) = 0
           do ispin = 1,2
               ishift = -63
               do l=1,Nint
                   ishift = ishift + 64
                   if (det1(l,ispin) == det2(l,ispin)) cycle
                   tmp = xor( det1(l,ispin), det2(l,ispin) )
                   particle = iand(tmp, det2(l,ispin))
                   hole = iand(tmp, det1(l,ispin))
                   if (particle /= 0_8) then
                       tz = trailz(particle)
                       exc(0,2,ispin) = 1
                       exc(1,2,ispin) = tz+ishift
                   end if
                   if (hole /= 0_8) then
                       tz = trailz(hole)
                       exc(0,1,ispin) = 1
                       exc(1,1,ispin) = tz+ishift
                   end if
                    phase=1
!                   if ( iand(exc(0,1,ispin),exc(0,2,ispin)) == 1 ) then
!                       low = min(exc(1,1,ispin),exc(1,2,ispin))
!
!                       high = max(exc(1,1,ispin),exc(1,2,ispin))
!                       j = ishft(low-1,-6)+1
!                       n = iand(low,63)
!                       k = ishft(high-1,-6)+1
!                       m = iand(high,63)
!                       if (j==k) then
!                           nperm = popcnt(iand(det1(j,ispin), &
!                                   iand( ibset(0_8,m-1)-1_8, ibclr(-1_8,n)+1_8 ) ))
!                       else
!                           nperm = popcnt(iand(det1(k,ispin), ibset(0_8,m-1)-1_8)) + &
!                                   popcnt(iand(det1(j,ispin), ibclr(-1_8,n) +1_8))
!                           do i=j+1,k-1
!                               nperm = nperm + popcnt(det1(i,ispin))
!                           end do
!                       end if
!                       phase = phase_dble(iand(nperm,1))
!                       return
!                   end if
               end do
           end do
       end


       subroutine get_excitation(det1,det2,exc,degree,phase,Nint)
           implicit none
           integer, intent(in) :: Nint
           integer*8, intent(in) :: det1(Nint,2), det2(Nint,2)
           integer, intent(out) :: exc(0:2,2,2)
           integer, intent(out) :: degree
           double precision, intent(out) :: phase



           degree = excitations(det1,det2,Nint)

           select case (degree)

           case (3:)
               degree = -1
               return

           case (2)
               call get_double_excitation(det1,det2,exc,phase,Nint)
               return

           case (1)
               call get_single_excitation(det1,det2,exc,phase,Nint)
               return

           case(0)
               return

           end select
       end


   end module