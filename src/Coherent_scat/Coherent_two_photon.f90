Program Coherent_scattering
    use calculate_form_factors
    implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER,PARAMETER :: nq=100

        integer(kind=ikind) ::natoms,i,j,count,ab,ang,cd

        real(kind=dp),dimension(:),allocatable:: bond

        real(kind=dp),dimension(:,:),allocatable::vbond,multaff1,multaff2,multiff1,multiff2
        real(kind=dp), dimension(nq):: q1,q2
        character(len=10),dimension(:), allocatable ::atoms
        real(kind=dp),dimension(:),allocatable :: x,y,z
        real(kind=dp),dimension(:,:),allocatable:: aff1,aff2,iff1,iff2
        real(kind=dp),dimension(:),allocatable ::ff
        character(len=10),dimension(151):: atoms_table
        character(len=10),dimension(55) :: atoms_table_2
        real(kind=dp),dimension(13):: q_table_2
        real(kind=dp),dimension(55,13):: iff
        real(kind=dp),dimension(nq) :: iff_fin
        real(kind=dp),dimension(:,:), allocatable :: iff_interp
        real(kind=dp), dimension(4,151) :: a,b
        real(kind=dp),  dimension(151):: c
        real(kind=dp), dimension(nq) :: f1,f2
        real(kind=dp),dimension(180) :: angleq
        real(kind=dp),dimension(nq,nq,180) :: sigma
        real(kind=dp) :: anglebond,pi,q12p,q12m,tt


        pi=dacos(-1.00_dp)
        call linspace(0.000001_dp,4.0_dp,nq,q1)
        call linspace(0.000001_dp,4.0_dp,nq,q2)
        call linspace(0.00001_dp,pi,180,angleq)

        open(unit=16,file="mol.xyz")
        read(16,*)natoms

        allocate(atoms(natoms), x(natoms), y(natoms), z(natoms),aff1(nq,natoms),aff2(nq,natoms),ff(nq))
        allocate(bond((natoms**2-natoms)/2),vbond(3,(natoms**2-natoms)/2),multaff1(nq,(natoms**2-natoms)/2))
        allocate(multaff2(nq,(natoms**2-natoms)/2),iff1(nq,natoms),iff2(nq,natoms), multiff1(nq,(natoms**2-natoms)/2))
        allocate(multiff2(nq,(natoms**2-natoms)/2))
        call table_of_ff("form_factors.dat",atoms_table,a,b,c)
        call table_of_iff("incoherent_factors.dat",atoms_table_2,iff,q_table_2)


        do i=1,natoms
            read(16,*)atoms(i),x(i),y(i),z(i)
            call obtain_form_factors(ff,atoms(i),q1,atoms_table,a,b,c)
            aff1(:,i)=ff
            call obtain_form_factors(ff,atoms(i),q2,atoms_table,a,b,c)
            aff2(:,i)=ff
            call obtain_i_form_factors(atoms_table_2,atoms(i),iff, q1,q_table_2,iff_fin)
            iff1(:,i)=iff_fin
            call obtain_i_form_factors(atoms_table_2,atoms(i),iff, q2,q_table_2,iff_fin)
            iff2(:,i)=iff_fin


        end do
        close(16)
        f1=0.0_dp
        f2=0.0_dp
        count=0
        do i=1,natoms
            f1=f1+abs(aff1(:,i))**2+iff1(:,i)
            f2=f2+abs(aff2(:,i))**2+iff2(:,i)
            do j=i+1,natoms
                count=count+1
                bond(count)=norm2((/x(i),y(i),z(i)/)-(/x(j),y(j),z(j)/))
                vbond(:,count)=(/x(i),y(i),z(i)/)-(/x(j),y(j),z(j)/)
                multaff1(:,count)= aff1(:,i)*aff1(:,j)
                multaff2(:,count)= aff2(:,i)*aff2(:,j)
                multiff1(:,count)= iff1(:,i)*iff1(:,j)
                multiff2(:,count)= iff2(:,i)*iff2(:,j)

            end do
        end do
        sigma=0.0_dp



        print*,"TTRA result"

       ! print*,'form_factors squared',f1(1),f2(1)
        do i=1,nq
               print*,100.00_dp-(100.00_dp*float(nq-i)/float(nq)), '%'
            do j=1,nq
                do ang=1,180
                    q12p=sqrt(q1(i)**2+q2(j)**2+2*q1(i)*q2(j)*cos(angleq(ang)))
                    q12m=sqrt(q1(i)**2+q2(j)**2-2*q1(i)*q2(j)*cos(angleq(ang)))

                    sigma(i,j,ang)=f1(i)*f2(j)
                   if (i==1 .and. j==1 .and. ang==1) then


                    end if

                    do ab=1,count
                        sigma(i,j,ang)=sigma(i,j,ang)+2*multaff1(i,ab)*multaff2(j,ab)*(sinc(q12p*bond(ab)) &
                        +sinc(q12m*bond(ab)))
                         if (i==1 .and. j==1 .and. ang==1) then


                        end if
                        sigma(i,j,ang)=sigma(i,j,ang)+2*(f2(j)*multaff1(i,ab)*sinc(q1(i)*bond(ab)) &
                        + f1(i) * multaff2(j,ab) * sinc(q2(j)*bond(ab)))
                        if (i==1 .and. j==1 .and. ang==1) then


                        end if

                        do cd=ab+1,count

                            anglebond=dot_product(vbond(:,ab),vbond(:,cd))/(bond(ab)*bond(cd))
                            if (isnan(anglebond)) then

                            end if
                            tt=TTRA(bond(ab),bond(cd),q1(i),q2(j),anglebond,angleq(ang))


                            sigma(i,j,ang)=sigma(i,j,ang)+(4*multaff1(i,ab)*multaff2(j,cd)*tt)
                             if (i==1 .and. j==1 .and. ang==1) then


                            end if


                            anglebond=dot_product(vbond(:,cd),vbond(:,ab))/(bond(ab)*bond(cd))
                            if (isnan(anglebond)) then

                            end if
                            tt=TTRA(bond(cd),bond(ab),q1(i),q2(j),anglebond,angleq(ang))


                            sigma(i,j,ang)=sigma(i,j,ang)+(4*multaff1(i,cd)*multaff2(j,ab)*tt)
                            if (i==1 .and. j==1 .and. ang==1) then


                            end if
                        end do
                    end do

                end do
            end do
    end do


print*,sigma(10,10,10)






End Program Coherent_scattering




