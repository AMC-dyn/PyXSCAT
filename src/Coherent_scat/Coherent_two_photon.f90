Program Coherent_scattering
    use calculate_form_factors
    implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER,PARAMETER :: nq=25

        integer(kind=ikind) ::natoms,i,j,count

        real(kind=dp),dimension(:),allocatable:: bond

        real(kind=dp),dimension(:,:),allocatable::vbond,multaff1,multaff2
        real(kind=dp), dimension(nq):: q1,q2
        character(len=10),dimension(:), allocatable ::atoms
        real(kind=dp),dimension(:),allocatable :: x,y,z
        real(kind=dp),dimension(:,:),allocatable:: aff1,aff2,iff1,iff2
        real(kind=dp),dimension(:),allocatable ::ff
        character(len=10),dimension(151):: atoms_table
        character(len=10),dimension(55) :: atoms_table_2
        real(kind=dp),dimension(12):: q_table_2
        real(kind=dp),dimension(55,12):: iff
        real(kind=dp),dimension(nq) :: iff_fin
        real(kind=dp),dimension(:,:), allocatable :: iff_interp
        real(kind=dp), dimension(4,151) :: a,b
        real(kind=dp),  dimension(151):: c
        real(kind=dp), dimension(nq) :: q12p, q12m,f1,f2,s1,s2
        real(kind=dp),dimension(180) :: angleq
        real(kind=dp),dimension(nq,nq,180) :: sigma
        real(kind=dp) :: anglebond



        call linspace(0.00001_dp,4.0_dp,nq,q1)
        call linspace(0.00001_dp,4.0_dp,nq,q2)
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
        f1=sum(abs(aff1)**2.0_dp,dim=2)
        f2=sum(abs(aff2)**2.0_dp,dim=2)
        s1=sum(iff1,dim=2)
        s2=sum(iff2,dim=2)
        count=0
        do i=1,natoms
            do j=i+1,natoms
                bond(count)=norm2((/x(i),y(i),z(i)/)-(/x(j),y(j),z(j)/))
                vbond(:,count)=(/x(i),y(i),z(i)/)-(/x(j),y(j),z(j)/)
                multaff1(:,count)= aff1(:,i)*aff1(:,j)
                multaff2(:,count)= aff2(:,i)*aff2(:,j)
                multiffi1(:,count)= iff1(:,i)*iff1(:,j)
                multiff2(:,count)= iff2(:,i)*iff2(:,j)
                count=count+1
            end do
        end do











End Program Coherent_scattering




