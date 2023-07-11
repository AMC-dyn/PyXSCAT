Module variables

   use uniquemodule
   use integrals
    implicit none
    contains

SUBROUTINE fill_md_table(D, l, x, ga)

        implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        ! The MD table to be populated D(Ngto, Ngto, 2maxl+1,maxl+1,maxl+1)
        REAL(kind=dp), INTENT(OUT), DIMENSION(:,:,:,:,:)    :: D
        ! vectors that contains the x coordinates for all GTOs, and all gammas
        REAl(kind=dp), INTENT(IN), DIMENSION(:),allocatable             :: ga, x
        ! vector that contains all the angular momenta (lx) for the GTOs
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:),allocatable       :: l


        ! loop and temp vars
        INTEGER(kind=ikind) :: i, j, ii, jj, N, maxl, l1, l2
        REAL(kind=dp)       :: gaP, Px, PA, PB, a

        ! number of GTOs
        N = size(l)
        ! maximum angular momentum
        maxl = maxval(l)


        ! the offdiagonal elements
        do i = 1,N
            do j = i+1, N
                gaP=ga(i)+ga(j)
                Px=(ga(i)*x(i)+ga(j)*x(j))/gaP
                ! order the angular momenta so that l1>=l2
                if (l(i)<l(j)) then
                    ii=j
                    jj=i
                else
                    ii=i
                    jj=j
                end if
                PA=Px-x(ii)
                PB=Px-x(jj)

                l1=l(ii)
                l2=l(jj)

                if (l1==0 .and. l2==0) then
                    D(ii,jj,1,1,1)=1
                    D(jj,ii,1,1,1)=1

                elseif (l1==1 .and. l2==0) then
                   D(ii,jj,1,2,1)=PA
                   D(ii,jj,2,2,1)=0.5/gaP
                   D(jj,ii,1,1,2)=D(ii,jj,1,2,1)
                   D(jj,ii,2,1,2)=D(ii,jj,2,2,1)

               elseif (l1==1 .and. l2==1) then
                    a=0.5/gaP
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a

                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)

                    D(jj,ii,1,2,2)=D(ii,jj,1,2,2)
                    D(jj,ii,2,2,2)=D(ii,jj,2,2,2)
                    D(jj,ii,3,2,2)=D(ii,jj,3,2,2)

                elseif (l1==2 .and. l2==0) then
                    a=0.5/gaP

                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a

                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)

                    D(jj,ii,1,1,3)=D(ii,jj,1,3,1)
                    D(jj,ii,2,1,3)=D(ii,jj,2,3,1)
                    D(jj,ii,3,1,3)=D(ii,jj,3,3,1)

                elseif (l1==2 .and. l2==1 ) then
                    a=0.5/gaP
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)

                    D(ii,jj,1,3,2)=PB*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,3,2)=a*D(ii,jj,1,3,1)+PB*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,3,2)=a*D(ii,jj,2,3,1)+PB*D(ii,jj,3,3,1)
                    D(ii,jj,4,3,2)=a*D(ii,jj,3,3,1)

                    D(jj,ii,1,2,3)=D(ii,jj,1,3,2)
                    D(jj,ii,2,2,3)=D(ii,jj,2,3,2)
                    D(jj,ii,3,2,3)=D(ii,jj,3,3,2)
                    D(jj,ii,4,2,3)=D(ii,jj,4,3,2)

                elseif (l1==2 .and. l2==2) then
                    a=0.5/gaP

                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,1,3)=PB*D(ii,jj,1,1,2)+D(ii,jj,2,1,2)
                    D(ii,jj,2,1,3)=a*D(ii,jj,1,1,2)+PB*D(ii,jj,2,1,2)
                    D(ii,jj,3,1,3)=a*D(ii,jj,2,1,2)
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,2,3)=PB*D(ii,jj,1,2,2)+D(ii,jj,2,2,2)
                    D(ii,jj,2,2,3)=a*D(ii,jj,1,2,2)+PB*D(ii,jj,2,2,2)+2.*D(ii,jj,3,2,2)
                    D(ii,jj,3,2,3)=a*D(ii,jj,2,2,2)+PB*D(ii,jj,3,2,2)
                    D(ii,jj,4,2,3)=a*D(ii,jj,3,2,2)

                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,3,2)=PB*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,3,2)=a*D(ii,jj,1,3,1)+PB*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,3,2)=a*D(ii,jj,2,3,1)+PB*D(ii,jj,3,3,1)
                    D(ii,jj,4,3,2)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,3,3)=PB*D(ii,jj,1,3,2)+D(ii,jj,2,3,2)
                    D(ii,jj,2,3,3)=a*D(ii,jj,1,3,2)+PB*D(ii,jj,2,3,2)+2.*D(ii,jj,3,3,2)
                    D(ii,jj,3,3,3)=a*D(ii,jj,2,3,2)+PB*D(ii,jj,3,3,2)+3.*D(ii,jj,4,3,2)
                    D(ii,jj,4,3,3)=a*D(ii,jj,3,3,2)+PB*D(ii,jj,4,3,2)
                    D(ii,jj,5,3,3)=a*D(ii,jj,4,3,2)

                    D(jj,ii,1,3,3)=D(ii,jj,1,3,3)
                    D(jj,ii,2,3,3)=D(ii,jj,2,3,3)
                    D(jj,ii,3,3,3)=D(ii,jj,3,3,3)
                    D(jj,ii,4,3,3)=D(ii,jj,4,3,3)
                    D(jj,ii,5,3,3)=D(ii,jj,5,3,3)

                elseif (l1==3 .and. l2==0) then
                    a=0.5/gaP
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,4,1)=PA*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,4,1)=a*D(ii,jj,1,3,1)+PA*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,4,1)=a*D(ii,jj,2,3,1)+PA*D(ii,jj,3,3,1)
                    D(ii,jj,4,4,1)=a*D(ii,jj,3,3,1)

                    D(jj,ii,1,1,4)=D(ii,jj,1,4,1)
                    D(jj,ii,2,1,4)=D(ii,jj,2,4,1)
                    D(jj,ii,3,1,4)=D(ii,jj,3,4,1)
                    D(jj,ii,4,1,4)=D(ii,jj,4,4,1)

                elseif (l1==3 .and. l2==1) then
                    a=0.5/gaP
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)

                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,3,2)=PB*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,3,2)=a*D(ii,jj,1,3,1)+PB*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,3,2)=a*D(ii,jj,2,3,1)+PB*D(ii,jj,3,3,1)
                    D(ii,jj,4,3,2)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,4,1)=PA*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,4,1)=a*D(ii,jj,1,3,1)+PA*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,4,1)=a*D(ii,jj,2,3,1)+PA*D(ii,jj,3,3,1)
                    D(ii,jj,4,4,1)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,4,2)=PB*D(ii,jj,1,4,1)+D(ii,jj,2,4,1)
                    D(ii,jj,2,4,2)=a*D(ii,jj,1,4,1)+PB*D(ii,jj,2,4,1)+2.*D(ii,jj,3,4,1)
                    D(ii,jj,3,4,2)=a*D(ii,jj,2,4,1)+PB*D(ii,jj,3,4,1)+3.*D(ii,jj,4,4,1)
                    D(ii,jj,4,4,2)=a*D(ii,jj,3,4,1)+PB*D(ii,jj,4,4,1)
                    D(ii,jj,5,4,2)=a*D(ii,jj,4,4,1)


                    D(jj,ii,1,2,4)=D(ii,jj,1,4,2)
                    D(jj,ii,2,2,4)=D(ii,jj,2,4,2)
                    D(jj,ii,3,2,4)=D(ii,jj,3,4,2)
                    D(jj,ii,4,2,4)=D(ii,jj,4,4,2)
                    D(jj,ii,5,2,4)=D(ii,jj,5,4,2)
                    ! need to add symmetric combination

                elseif (l1==3 .and. l2==2) then
                    a=0.5/gaP
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,1,3)=PB*D(ii,jj,1,1,2)+D(ii,jj,2,1,2)
                    D(ii,jj,2,1,3)=a*D(ii,jj,1,1,2)+PB*D(ii,jj,2,1,2)
                    D(ii,jj,3,1,3)=a*D(ii,jj,2,1,2)
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,2,3)=PB*D(ii,jj,1,2,2)+D(ii,jj,2,2,2)
                    D(ii,jj,2,2,3)=a*D(ii,jj,1,2,2)+PB*D(ii,jj,2,2,2)+2.*D(ii,jj,3,2,2)
                    D(ii,jj,3,2,3)=a*D(ii,jj,2,2,2)+PB*D(ii,jj,3,2,2)
                    D(ii,jj,4,2,3)=a*D(ii,jj,3,2,2)
                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,3,2)=PB*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,3,2)=a*D(ii,jj,1,3,1)+PB*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,3,2)=a*D(ii,jj,2,3,1)+PB*D(ii,jj,3,3,1)
                    D(ii,jj,4,3,2)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,3,3)=PB*D(ii,jj,1,3,2)+D(ii,jj,2,3,2)
                    D(ii,jj,2,3,3)=a*D(ii,jj,1,3,2)+PB*D(ii,jj,2,3,2)+2.*D(ii,jj,3,3,2)
                    D(ii,jj,3,3,3)=a*D(ii,jj,2,3,2)+PB*D(ii,jj,3,3,2)+3.*D(ii,jj,4,3,2)
                    D(ii,jj,4,3,3)=a*D(ii,jj,3,3,2)+PB*D(ii,jj,4,3,2)
                    D(ii,jj,5,3,3)=a*D(ii,jj,4,3,2)
                    D(ii,jj,1,4,1)=PA*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,4,1)=a*D(ii,jj,1,3,1)+PA*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,4,1)=a*D(ii,jj,2,3,1)+PA*D(ii,jj,3,3,1)
                    D(ii,jj,4,4,1)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,4,2)=PB*D(ii,jj,1,4,1)+D(ii,jj,2,4,1)
                    D(ii,jj,2,4,2)=a*D(ii,jj,1,4,1)+PB*D(ii,jj,2,4,1)+2.*D(ii,jj,3,4,1)
                    D(ii,jj,3,4,2)=a*D(ii,jj,2,4,1)+PB*D(ii,jj,3,4,1)+3.*D(ii,jj,4,4,1)
                    D(ii,jj,4,4,2)=a*D(ii,jj,3,4,1)+PB*D(ii,jj,4,4,1)
                    D(ii,jj,5,4,2)=a*D(ii,jj,4,4,1)
                    D(ii,jj,1,4,3)=PB*D(ii,jj,1,4,2)+D(ii,jj,2,4,2)
                    D(ii,jj,2,4,3)=a*D(ii,jj,1,4,2)+PB*D(ii,jj,2,4,2)+2.*D(ii,jj,3,4,2)
                    D(ii,jj,3,4,3)=a*D(ii,jj,2,4,2)+PB*D(ii,jj,3,4,2)+3.*D(ii,jj,4,4,2)
                    D(ii,jj,4,4,3)=a*D(ii,jj,3,4,2)+PB*D(ii,jj,4,4,2)+4.*D(ii,jj,5,4,2)
                    D(ii,jj,5,4,3)=a*D(ii,jj,4,4,2)+PB*D(ii,jj,5,4,2)
                    D(ii,jj,6,4,3)=a*D(ii,jj,5,4,2)


                    D(jj,ii,1,3,4)=D(ii,jj,1,4,3)
                    D(jj,ii,2,3,4)=D(ii,jj,2,4,3)
                    D(jj,ii,3,3,4)=D(ii,jj,3,4,3)
                    D(jj,ii,4,3,4)=D(ii,jj,4,4,3)
                    D(jj,ii,5,3,4)=D(ii,jj,5,4,3)
                    D(jj,ii,6,3,4)=D(ii,jj,6,4,3)

                    ! need to add symmetric combination
                elseif (l1==3 .and. l2==3) then
                    a=0.5/gaP
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,1,1,2)=PB
                    D(ii,jj,2,1,2)=a
                    D(ii,jj,1,1,3)=PB*D(ii,jj,1,1,2)+D(ii,jj,2,1,2)
                    D(ii,jj,2,1,3)=a*D(ii,jj,1,1,2)+PB*D(ii,jj,2,1,2)
                    D(ii,jj,3,1,3)=a*D(ii,jj,2,1,2)
                    D(ii,jj,1,1,4)=PB*D(ii,jj,1,1,3)+D(ii,jj,2,1,3)
                    D(ii,jj,2,1,4)=a*D(ii,jj,1,1,3)+PB*D(ii,jj,2,1,3)+2.*D(ii,jj,3,1,3)
                    D(ii,jj,3,1,4)=a*D(ii,jj,2,1,3)+PB*D(ii,jj,3,1,3)
                    D(ii,jj,4,1,4)=a*D(ii,jj,3,1,3)
                    D(ii,jj,1,2,1)=PA
                    D(ii,jj,2,2,1)=a
                    D(ii,jj,1,2,2)=PB*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,2,2)=a*D(ii,jj,1,2,1)+PB*D(ii,jj,2,2,1)
                    D(ii,jj,3,2,2)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,2,3)=PB*D(ii,jj,1,2,2)+D(ii,jj,2,2,2)
                    D(ii,jj,2,2,3)=a*D(ii,jj,1,2,2)+PB*D(ii,jj,2,2,2)+2.*D(ii,jj,3,2,2)
                    D(ii,jj,3,2,3)=a*D(ii,jj,2,2,2)+PB*D(ii,jj,3,2,2)
                    D(ii,jj,4,2,3)=a*D(ii,jj,3,2,2)

                    D(ii,jj,1,2,4)=PB*D(ii,jj,1,2,3)+D(ii,jj,2,2,3)
                    D(ii,jj,2,2,4)=a*D(ii,jj,1,2,3)+PB*D(ii,jj,2,2,3)+2.*D(ii,jj,3,2,3)
                    D(ii,jj,3,2,4)=a*D(ii,jj,2,2,3)+PB*D(ii,jj,3,2,3)+3.*D(ii,jj,4,2,3)
                    D(ii,jj,4,2,4)=a*D(ii,jj,3,2,3)+PB*D(ii,jj,4,2,3)
                    D(ii,jj,5,2,4)=a*D(ii,jj,4,2,3)
                    D(ii,jj,1,3,1)=PA*D(ii,jj,1,2,1)+D(ii,jj,2,2,1)
                    D(ii,jj,2,3,1)=a*D(ii,jj,1,2,1)+PA*D(ii,jj,2,2,1)
                    D(ii,jj,3,3,1)=a*D(ii,jj,2,2,1)
                    D(ii,jj,1,3,2)=PB*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,3,2)=a*D(ii,jj,1,3,1)+PB*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,3,2)=a*D(ii,jj,2,3,1)+PB*D(ii,jj,3,3,1)
                    D(ii,jj,4,3,2)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,3,3)=PB*D(ii,jj,1,3,2)+D(ii,jj,2,3,2)
                    D(ii,jj,2,3,3)=a*D(ii,jj,1,3,2)+PB*D(ii,jj,2,3,2)+2.*D(ii,jj,3,3,2)
                    D(ii,jj,3,3,3)=a*D(ii,jj,2,3,2)+PB*D(ii,jj,3,3,2)+3.*D(ii,jj,4,3,2)
                    D(ii,jj,4,3,3)=a*D(ii,jj,3,3,2)+PB*D(ii,jj,4,3,2)
                    D(ii,jj,5,3,3)=a*D(ii,jj,4,3,2)

                    D(ii,jj,1,3,4)=PB*D(ii,jj,1,3,3)+D(ii,jj,2,3,3)
                    D(ii,jj,2,3,4)=a*D(ii,jj,1,3,3)+PB*D(ii,jj,2,3,3)+2.*D(ii,jj,3,3,3)
                    D(ii,jj,3,3,4)=a*D(ii,jj,2,3,3)+PB*D(ii,jj,3,3,3)+3.*D(ii,jj,4,3,3)
                    D(ii,jj,4,3,4)=a*D(ii,jj,3,3,3)+PB*D(ii,jj,4,3,3)+4.*D(ii,jj,5,3,3)
                    D(ii,jj,5,3,4)=a*D(ii,jj,4,3,3)+PB*D(ii,jj,5,3,3)
                    D(ii,jj,6,3,4)=a*D(ii,jj,5,3,3)
                    D(ii,jj,1,4,1)=PA*D(ii,jj,1,3,1)+D(ii,jj,2,3,1)
                    D(ii,jj,2,4,1)=a*D(ii,jj,1,3,1)+PA*D(ii,jj,2,3,1)+2.*D(ii,jj,3,3,1)
                    D(ii,jj,3,4,1)=a*D(ii,jj,2,3,1)+PA*D(ii,jj,3,3,1)
                    D(ii,jj,4,4,1)=a*D(ii,jj,3,3,1)
                    D(ii,jj,1,4,2)=PB*D(ii,jj,1,4,1)+D(ii,jj,2,4,1)
                    D(ii,jj,2,4,2)=a*D(ii,jj,1,4,1)+PB*D(ii,jj,2,4,1)+2.*D(ii,jj,3,4,1)
                    D(ii,jj,3,4,2)=a*D(ii,jj,2,4,1)+PB*D(ii,jj,3,4,1)+3.*D(ii,jj,4,4,1)
                    D(ii,jj,4,4,2)=a*D(ii,jj,3,4,1)+PB*D(ii,jj,4,4,1)
                    D(ii,jj,5,4,2)=a*D(ii,jj,4,4,1)

                    D(ii,jj,1,4,3)=PB*D(ii,jj,1,4,2)+D(ii,jj,2,4,2)
                    D(ii,jj,2,4,3)=a*D(ii,jj,1,4,2)+PB*D(ii,jj,2,4,2)+2.*D(ii,jj,3,4,2)
                    D(ii,jj,3,4,3)=a*D(ii,jj,2,4,2)+PB*D(ii,jj,3,4,2)+3.*D(ii,jj,4,4,2)
                    D(ii,jj,4,4,3)=a*D(ii,jj,3,4,2)+PB*D(ii,jj,4,4,2)+4.*D(ii,jj,5,4,2)
                    D(ii,jj,5,4,3)=a*D(ii,jj,4,4,2)+PB*D(ii,jj,5,4,2)
                    D(ii,jj,6,4,3)=a*D(ii,jj,5,4,2)


                    D(ii,jj,1,4,4)=PB*D(ii,jj,1,4,3)+D(ii,jj,2,4,3)
                    D(ii,jj,2,4,4)=a*D(ii,jj,1,4,3)+PB*D(ii,jj,2,4,3)+2.*D(ii,jj,3,4,3)
                    D(ii,jj,3,4,4)=a*D(ii,jj,2,4,3)+PB*D(ii,jj,3,4,3)+3.*D(ii,jj,4,4,3)
                    D(ii,jj,4,4,4)=a*D(ii,jj,3,4,3)+PB*D(ii,jj,4,4,3)+4.*D(ii,jj,5,4,3)
                    D(ii,jj,5,4,4)=a*D(ii,jj,4,4,3)+PB*D(ii,jj,5,4,3)+5.*D(ii,jj,6,4,3)
                    D(ii,jj,6,4,4)=a*D(ii,jj,5,4,3)+PB*D(ii,jj,6,4,3)
                    D(ii,jj,7,4,4)=a*D(ii,jj,6,4,3)

                    D(jj,ii,1,4,4)=D(ii,jj,1,4,4)
                    D(jj,ii,2,4,4)=D(ii,jj,2,4,4)
                    D(jj,ii,3,4,4)=D(ii,jj,3,4,4)
                    D(jj,ii,4,4,4)=D(ii,jj,4,4,4)
                    D(jj,ii,5,4,4)=D(ii,jj,5,4,4)
                    D(jj,ii,6,4,4)=D(ii,jj,6,4,4)
                    D(jj,ii,7,4,4)=D(ii,jj,7,4,4)
                    ! need to add symmetric combination
                else
                    print*, "case not programmed yet: l1/2= ", l1, l2
                    stop
                end if

            end do
        end do


        ! diagonal case

        do i = 1,N
            j=i

            gaP=ga(i)+ga(j)
            Px=(ga(i)*x(i)+ga(j)*x(j))/gaP
            PA=Px-x(i)
            PB=Px-x(j)

            l1=l(i)
            l2=l(j)

            if (l1==0 .and. l2==0) then
                D(i,j,1,1,1)=1

            elseif (l1==1) then
                a=0.5/gaP
                D(i,j,1,1,2)=PB
                D(i,j,2,1,2)=a
                D(i,j,1,2,1)=PA
                D(i,j,2,2,1)=a

                D(i,j,1,2,2)=PB*D(i,j,1,2,1)+D(i,j,2,2,1)
                D(i,j,2,2,2)=a*D(i,j,1,2,1)+PB*D(i,j,2,2,1)
                D(i,j,3,2,2)=a*D(i,j,2,2,1)

            elseif (l1==2 ) then
                a=0.5/gaP

                D(i,j,1,1,2)=PB
                D(i,j,2,1,2)=a
                D(i,j,1,1,3)=PB*D(i,j,1,1,2)+D(i,j,2,1,2)
                D(i,j,2,1,3)=a*D(i,j,1,1,2)+PB*D(i,j,2,1,2)
                D(i,j,3,1,3)=a*D(i,j,2,1,2)
                D(i,j,1,2,1)=PA
                D(i,j,2,2,1)=a
                D(i,j,1,2,2)=PB*D(i,j,1,2,1)+D(i,j,2,2,1)
                D(i,j,2,2,2)=a*D(i,j,1,2,1)+PB*D(i,j,2,2,1)
                D(i,j,3,2,2)=a*D(i,j,2,2,1)
                D(i,j,1,2,3)=PB*D(i,j,1,2,2)+D(i,j,2,2,2)
                D(i,j,2,2,3)=a*D(i,j,1,2,2)+PB*D(i,j,2,2,2)+2.*D(i,j,3,2,2)
                D(i,j,3,2,3)=a*D(i,j,2,2,2)+PB*D(i,j,3,2,2)
                D(i,j,4,2,3)=a*D(i,j,3,2,2)

                D(i,j,1,3,1)=PA*D(i,j,1,2,1)+D(i,j,2,2,1)
                D(i,j,2,3,1)=a*D(i,j,1,2,1)+PA*D(i,j,2,2,1)
                D(i,j,3,3,1)=a*D(i,j,2,2,1)
                D(i,j,1,3,2)=PB*D(i,j,1,3,1)+D(i,j,2,3,1)
                D(i,j,2,3,2)=a*D(i,j,1,3,1)+PB*D(i,j,2,3,1)+2.*D(i,j,3,3,1)
                D(i,j,3,3,2)=a*D(i,j,2,3,1)+PB*D(i,j,3,3,1)
                D(i,j,4,3,2)=a*D(i,j,3,3,1)
                D(i,j,1,3,3)=PB*D(i,j,1,3,2)+D(i,j,2,3,2)
                D(i,j,2,3,3)=a*D(i,j,1,3,2)+PB*D(i,j,2,3,2)+2.*D(i,j,3,3,2)
                D(i,j,3,3,3)=a*D(i,j,2,3,2)+PB*D(i,j,3,3,2)+3.*D(i,j,4,3,2)
                D(i,j,4,3,3)=a*D(i,j,3,3,2)+PB*D(i,j,4,3,2)
                D(i,j,5,3,3)=a*D(i,j,4,3,2)

            elseif (l1==3 ) then
                    a=0.5/gaP
                    D(i,j,1,1,2)=PB
                    D(i,j,1,1,2)=PB
                    D(i,j,2,1,2)=a
                    D(i,j,1,1,3)=PB*D(i,j,1,1,2)+D(i,j,2,1,2)
                    D(i,j,2,1,3)=a*D(i,j,1,1,2)+PB*D(i,j,2,1,2)
                    D(i,j,3,1,3)=a*D(i,j,2,1,2)
                    D(i,j,1,1,4)=PB*D(i,j,1,1,3)+D(i,j,2,1,3)
                    D(i,j,2,1,4)=a*D(i,j,1,1,3)+PB*D(i,j,2,1,3)+2.*D(i,j,3,1,3)
                    D(i,j,3,1,4)=a*D(i,j,2,1,3)+PB*D(i,j,3,1,3)
                    D(i,j,4,1,4)=a*D(i,j,3,1,3)
                    D(i,j,1,2,1)=PA
                    D(i,j,2,2,1)=a
                    D(i,j,1,2,2)=PB*D(i,j,1,2,1)+D(i,j,2,2,1)
                    D(i,j,2,2,2)=a*D(i,j,1,2,1)+PB*D(i,j,2,2,1)
                    D(i,j,3,2,2)=a*D(i,j,2,2,1)
                    D(i,j,1,2,3)=PB*D(i,j,1,2,2)+D(i,j,2,2,2)
                    D(i,j,2,2,3)=a*D(i,j,1,2,2)+PB*D(i,j,2,2,2)+2.*D(i,j,3,2,2)
                    D(i,j,3,2,3)=a*D(i,j,2,2,2)+PB*D(i,j,3,2,2)
                    D(i,j,4,2,3)=a*D(i,j,3,2,2)

                    D(i,j,1,2,4)=PB*D(i,j,1,2,3)+D(i,j,2,2,3)
                    D(i,j,2,2,4)=a*D(i,j,1,2,3)+PB*D(i,j,2,2,3)+2.*D(i,j,3,2,3)
                    D(i,j,3,2,4)=a*D(i,j,2,2,3)+PB*D(i,j,3,2,3)+3.*D(i,j,4,2,3)
                    D(i,j,4,2,4)=a*D(i,j,3,2,3)+PB*D(i,j,4,2,3)
                    D(i,j,5,2,4)=a*D(i,j,4,2,3)
                    D(i,j,1,3,1)=PA*D(i,j,1,2,1)+D(i,j,2,2,1)
                    D(i,j,2,3,1)=a*D(i,j,1,2,1)+PA*D(i,j,2,2,1)
                    D(i,j,3,3,1)=a*D(i,j,2,2,1)
                    D(i,j,1,3,2)=PB*D(i,j,1,3,1)+D(i,j,2,3,1)
                    D(i,j,2,3,2)=a*D(i,j,1,3,1)+PB*D(i,j,2,3,1)+2.*D(i,j,3,3,1)
                    D(i,j,3,3,2)=a*D(i,j,2,3,1)+PB*D(i,j,3,3,1)
                    D(i,j,4,3,2)=a*D(i,j,3,3,1)
                    D(i,j,1,3,3)=PB*D(i,j,1,3,2)+D(i,j,2,3,2)
                    D(i,j,2,3,3)=a*D(i,j,1,3,2)+PB*D(i,j,2,3,2)+2.*D(i,j,3,3,2)
                    D(i,j,3,3,3)=a*D(i,j,2,3,2)+PB*D(i,j,3,3,2)+3.*D(i,j,4,3,2)
                    D(i,j,4,3,3)=a*D(i,j,3,3,2)+PB*D(i,j,4,3,2)
                    D(i,j,5,3,3)=a*D(i,j,4,3,2)

                    D(i,j,1,3,4)=PB*D(i,j,1,3,3)+D(i,j,2,3,3)
                    D(i,j,2,3,4)=a*D(i,j,1,3,3)+PB*D(i,j,2,3,3)+2.*D(i,j,3,3,3)
                    D(i,j,3,3,4)=a*D(i,j,2,3,3)+PB*D(i,j,3,3,3)+3.*D(i,j,4,3,3)
                    D(i,j,4,3,4)=a*D(i,j,3,3,3)+PB*D(i,j,4,3,3)+4.*D(i,j,5,3,3)
                    D(i,j,5,3,4)=a*D(i,j,4,3,3)+PB*D(i,j,5,3,3)
                    D(i,j,6,3,4)=a*D(i,j,5,3,3)
                    D(i,j,1,4,1)=PA*D(i,j,1,3,1)+D(i,j,2,3,1)
                    D(i,j,2,4,1)=a*D(i,j,1,3,1)+PA*D(i,j,2,3,1)+2.*D(i,j,3,3,1)
                    D(i,j,3,4,1)=a*D(i,j,2,3,1)+PA*D(i,j,3,3,1)
                    D(i,j,4,4,1)=a*D(i,j,3,3,1)
                    D(i,j,1,4,2)=PB*D(i,j,1,4,1)+D(i,j,2,4,1)
                    D(i,j,2,4,2)=a*D(i,j,1,4,1)+PB*D(i,j,2,4,1)+2.*D(i,j,3,4,1)
                    D(i,j,3,4,2)=a*D(i,j,2,4,1)+PB*D(i,j,3,4,1)+3.*D(i,j,4,4,1)
                    D(i,j,4,4,2)=a*D(i,j,3,4,1)+PB*D(i,j,4,4,1)
                    D(i,j,5,4,2)=a*D(i,j,4,4,1)

                    D(i,j,1,4,3)=PB*D(i,j,1,4,2)+D(i,j,2,4,2)
                    D(i,j,2,4,3)=a*D(i,j,1,4,2)+PB*D(i,j,2,4,2)+2.*D(i,j,3,4,2)
                    D(i,j,3,4,3)=a*D(i,j,2,4,2)+PB*D(i,j,3,4,2)+3.*D(i,j,4,4,2)
                    D(i,j,4,4,3)=a*D(i,j,3,4,2)+PB*D(i,j,4,4,2)+4.*D(i,j,5,4,2)
                    D(i,j,5,4,3)=a*D(i,j,4,4,2)+PB*D(i,j,5,4,2)
                    D(i,j,6,4,3)=a*D(i,j,5,4,2)

                    D(i,j,1,4,4)=PB*D(i,j,1,4,3)+D(i,j,2,4,3)
                    D(i,j,2,4,4)=a*D(i,j,1,4,3)+PB*D(i,j,2,4,3)+2.*D(i,j,3,4,3)
                    D(i,j,3,4,4)=a*D(i,j,2,4,3)+PB*D(i,j,3,4,3)+3.*D(i,j,4,4,3)
                    D(i,j,4,4,4)=a*D(i,j,3,4,3)+PB*D(i,j,4,4,3)+4.*D(i,j,5,4,3)
                    D(i,j,5,4,4)=a*D(i,j,4,4,3)+PB*D(i,j,5,4,3)+5.*D(i,j,6,4,3)
                    D(i,j,6,4,4)=a*D(i,j,5,4,3)+PB*D(i,j,6,4,3)
                    D(i,j,7,4,4)=a*D(i,j,6,4,3)

            else

                    print*, "case not programmed yet: l1/2= " , l1, l2
                    stop
            end if
        end do

    END SUBROUTINE fill_md_table

subroutine variables_total(px,py,pz,ddx,ddy,ddz,z11,z22,z1t,z2t,&
        e12,maxl,ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq)


        implicit none


        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=ikind), intent(in)::ngto,ng, nmat, nq, maxl
        integer(kind=ikind), intent(in),dimension(:),allocatable :: l, m, n, group
        integer(kind=ikind), intent(in),dimension(:)::m1, m2, m3, m4,group_start,group_count

        real(kind=dp), intent(in),dimension(:),allocatable :: ga, xx, yy, zz, total,q
        real(kind=dp), dimension(ng):: ga2, xx2, yy2, zz2
        real(kind=dp), intent(in),dimension(:,:),allocatable :: mmod
        real(kind=dp), intent(out), dimension(ng,ng):: px,py,pz

        real(kind=dp), dimension(ngto,ngto) :: cmat
        !real(kind=dp), intent(out), dimension(ngto,ngto,ngto,ngto) :: zcontr

        real(kind=dp), intent(out), dimension(nq,ng,ng) :: e12

        real(kind=dp), intent(out), dimension(maxl*2+1,maxl+1,maxl+1,ng,ng) :: ddx,ddy,ddz


        integer(kind=ikind) :: max4l, max2l, i, j, k, N1, N2, nn1, nn2,ii,jj,ls, ms, ns,count
        integer(kind=ikind), dimension(ngto) :: ll
        integer(kind=ikind), dimension(nmat) :: indexing1,indexing2
        real(kind=dp), dimension(:,:,:),intent(out), allocatable :: z11, z22,z1t,z2t

        real(kind=dp) :: pi, gap,time1,time2,time3,time4

        integer(kind=ikind), dimension(:), allocatable :: iduplicates,jduplicates,m11,m22,m33,m44
        real(kind=dp),  dimension(ngto,ngto,maxl*2+1,maxl+1,maxl+1):: dx,dy,dz
        real(kind=dp),  dimension(nq):: preexp
        real(kind=dp),  dimension(:), allocatable:: temp1,temp2
        integer(kind=ikind),dimension(nmat,2) :: vec1, vec2
        integer(kind=ikind),dimension(nmat,4) :: mat1
        integer(kind=ikind),dimension(:,:),allocatable :: matfin
        real(kind=dp),dimension(:), allocatable :: totalfin
        real(kind=dp)   :: obwohl
        integer(kind=ikind) :: counter
        logical :: divided

        pi = acos(-1.0000)
        max4l=maxval(l)*4
        max2l=maxval(l)*2+1


       ! N1=size(ipos)
        do i=1,nmat

            vec1(i,:)=(/m1(i), m2(i)/)
            vec2(i,:)=(/m3(i), m4(i)/)
            call bubble_sort(vec1(i,:))
            call bubble_sort(vec2(i,:))
            mat1(i,1:2)=vec1(i,:)
            mat1(i,3:4)=vec2(i,:)

        end do
        print*,'wtf',nq,size(q)
        divided=.True.
        if (divided==.True.) then
            matfin=mat1
            totalfin=total
            else
            call unique_total(mat1,total,matfin,totalfin)
            print*,sum(total)
        end if

        !matfin=mat1
        !totalfin=total
        allocate(m11(size(totalfin)),m22(size(totalfin)),m33(size(totalfin)),m44(size(totalfin)))
        allocate(z11(size(totalfin),ngto,ngto), z22(size(totalfin),ngto,ngto), &
                z1t(ngto,size(totalfin),ngto), z2t(ngto,size(totalfin),ngto), &
               temp1(size(totalfin)),temp2(size(totalfin)))
        m11 = matfin(:,1)
        m22 = matfin(:,2)
        m33 = matfin(:,3)
        m44 = matfin(:,4)

        write(*,*)'shape of unique mat',shape(totalfin)
      !  print*,shape(z11), nmat
        z11=0.0_dp
        z22=0.0_dp
        z1t=0.0_dp
        z2t=0.0_dp
        print*,z11(1,1,1),m11(1)


        temp1=0.0
        temp2=0.0
        print*,ngto,shape(z11),shape(z22),shape(mmod)

        counter=0
         do  ii=1,ngto

            do jj=1,ngto


                        temp1 = totalfin * (mmod(m11, ii) * mmod(m22, jj) + mmod(m11, jj) * mmod(m22, ii))
                        temp2 = mmod(m33, ii) * mmod(m44, jj) + mmod(m33, jj) * mmod(m44, ii)

                        z11(:,ii, jj) =  temp1
                        z22(:, ii, jj) = temp2
                        z1t(ii,:,jj)=temp1
                        z2t(ii,:,jj)=temp2

                        if (sum(abs(temp1-temp2))<1E-10) then
                            counter=counter+1
                        end if



            enddo
        enddo


        obwohl = sum(abs(z11-z22))
        if (counter==ngto**2) then
            print*, 'Las gustones con las muchachas'
            else
            print*, 'Los gustones de los resultados, muchas',ngto,maxl,ng
        end if

        gap=0.0
        px=0.0
        py=0.0
        pz=0.0
        e12=0.0
        ddx=0.0
        ddy=0.0
        ddz=0.0
        dx=0.0
        dy=0.0
        dz=0.0
        print*,'hijos de la rata '



        call cpu_time(time3)

        print*,'hijos de la rata I'

        call fill_md_table(dx,l,xx,ga)
        call fill_md_table(dy,m,yy,ga)
        call fill_md_table(dz,n,zz,ga)

         print*,'hijos de la rata II'
       ! allocate( ga2(size(apos)), xx2(size(apos)), yy2(size(apos)), zz2(size(apos)) )

!        ll = l + m + n
!        ga2 = ga(apos)
!        xx2 = xx(apos)
!        yy2 = yy(apos)
!        zz2 = zz(apos)
!        ll2 = ll(apos)

     !   N2=size(apos)

       ! allocate(ddx(n2,n2,max2l,maxl,maxl),ddy(n2,n2,max2l,maxl,maxl),ddz(n2,n2,max2l,maxl,maxl))
        !allocate(preexp(nq))

       do jj = 1,Ng
            ! essentially taking the values for the first GTO in the group.
            ! All gtos in the group have the same prefactors as the prefactors
            ! do not't depend on l, m and n
            j = group_start(jj)
            do ii = 1, Ng
                i = group_start(ii)
                gaP=ga(i)+ga(j)
                Px(ii,jj)=(ga(i)*xx(i) + ga(j)*xx(j))/gaP
                Py(ii,jj)=(ga(i)*yy(i) + ga(j)*yy(j))/gaP
                Pz(ii,jj)=(ga(i)*zz(i) + ga(j)*zz(j))/gaP

                E12(:,ii,jj) = (pi/gaP)**1.5 * exp(-q*q*0.25/gaP) &
                    * exp(-ga(i)*ga(j)/gaP*((xx(i)-xx(j))**2. + (yy(i)-yy(j))**2. + (zz(i)-zz(j))**2.))
            end do
        end do
        print*,sum(Pz)


     do jj = 1, Ngto
            j = group(jj)
            do ii = 1, Ngto
                i = group(ii)
                do ls = 1, l(ii)+l(jj)+1
                    Ddx(ls, l(ii)+1,l(jj)+1,j,i)=Dx(ii,jj,ls, l(ii)+1,l(jj)+1)
                end do

                do ms = 1, m(ii)+m(jj)+1
                    Ddy(ms, m(ii)+1,m(jj)+1,j,i)=Dy(ii,jj,ms, m(ii)+1,m(jj)+1)
                end do

                do ns = 1, n(ii)+n(jj)+1
                    Ddz(ns, n(ii)+1,n(jj)+1,j,i)=Dz(ii,jj,ns, n(ii)+1,n(jj)+1)
                end do
            end do
        end do

    print*,ngto,size(E12(:,1,1)),size(E12(1,:,1)), size(E12(1,1,:))
    print*, 'leaving variables total'

    end subroutine variables_total





    subroutine variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl,ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,onerdm,maxnmo,q,nq)


    use uniquemodule
        implicit none


        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        integer(kind=ikind), intent(in)::ngto,ng,  nq, maxl,maxnmo
        integer(kind=ikind), intent(in),dimension(:),allocatable :: l, m, n,group
        integer(kind=ikind), intent(in),dimension(:):: group_count,group_start
        real(kind=dp), intent(in),dimension(:),allocatable :: ga, xx, yy, zz, q
        real(kind=dp), intent(in),dimension(:,:),allocatable :: mmod
        real(kind=dp), intent(in),dimension(:,:)::onerdm
        real(kind=dp), intent(out), dimension(ng,ng):: px,py,pz
        real(kind=dp), intent(out), dimension(ngto,ngto) :: z
        real(kind=dp), intent(out), dimension(nq,ng,ng) :: e12
        real(kind=dp), intent(out), dimension(maxl*2+1,maxl+1,maxl+1,ng,ng) :: ddx,ddy,ddz


        integer(kind=ikind) :: max4l, max2l, i, j, k, ii,jj,ls, ms, ns,count,cc,cc2
        real(kind=dp) :: pi, gap
        real(kind=dp),  dimension(ngto,ngto,maxl*2+1,maxl+1,maxl+1):: dx,dy,dz
        real(kind=dp),  dimension(nq):: preexp
        real(kind=dp):: temp



        integer(kind=ikind) :: counter

        pi = acos(-1.0000)
        max4l=maxval(l)*4
        max2l=maxval(l)*2+1


        z=0.0_dp
        print*,maxval(mmod), maxval(onerdm)
        counter=0
         do  ii=1,ngto

            do jj=1,ngto

                temp=0.0_dp

                do cc=1,maxnmo
                     do cc2=1,maxnmo


                      temp=temp+onerdm(cc,cc2)*mmod(cc,ii)*mmod(cc2,jj)
                      temp=temp+onerdm(cc,cc2)*mmod(cc,jj)*mmod(cc2,ii)

                    enddo
                enddo
                Z(ii,jj)=temp/2.0_dp;

            enddo
        enddo
        print*,'Z(1,2)',Z(1,2)


        gap=0.0
        dx=0.0
        dy=0.0
        dz=0.0
        px=0.0
        py=0.0
        pz=0.0
        e12=0.0
        ddx=0.0
        ddy=0.0
        ddz=0.0





        call fill_md_table(dx,l,xx,ga)
        call fill_md_table(dy,m,yy,ga)
        call fill_md_table(dz,n,zz,ga)


       do jj = 1,Ng
            ! essentially taking the values for the first GTO in the group.
            ! All gtos in the group have the same prefactors as the prefactors
            ! do not't depend on l, m and n
            j = group_start(jj)
            do ii = 1, Ng
                i = group_start(ii)
                gaP=ga(i)+ga(j)
                Px(ii,jj)=(ga(i)*xx(i) + ga(j)*xx(j))/gaP
                Py(ii,jj)=(ga(i)*yy(i) + ga(j)*yy(j))/gaP
                Pz(ii,jj)=(ga(i)*zz(i) + ga(j)*zz(j))/gaP
                E12(:,ii,jj) = (pi/gaP)**1.5 * exp(-q*q*0.25/gaP) &
                    * exp(-ga(i)*ga(j)/gaP*((xx(i)-xx(j))**2. + (yy(i)-yy(j))**2. + (zz(i)-zz(j))**2.))
            end do
        end do



     do jj = 1, Ngto
            j = group(jj)
            do ii = 1, Ngto
                i = group(ii)
                do ls = 1, l(ii)+l(jj)+1
                    Ddx(ls, l(ii)+1,l(jj)+1,j,i)=Dx(ii,jj,ls, l(ii)+1,l(jj)+1)
                end do

                do ms = 1, m(ii)+m(jj)+1
                    Ddy(ms, m(ii)+1,m(jj)+1,j,i)=Dy(ii,jj,ms, m(ii)+1,m(jj)+1)
                end do

                do ns = 1, n(ii)+n(jj)+1
                    Ddz(ns, n(ii)+1,n(jj)+1,j,i)=Dz(ii,jj,ns, n(ii)+1,n(jj)+1)
                end do
            end do
        end do

       print*,maxval(e12)
    end subroutine variables_elastic
End Module variables
