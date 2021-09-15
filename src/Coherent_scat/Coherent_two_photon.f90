Program Coherent_scattering
    use calculate_form_factors
    implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        INTEGER,PARAMETER :: nq=25

        integer(kind=ikind) ::natoms,i,j,count

        real(kind=dp),dimension(:),allocatable:: bond,bondabs
        real(kind=dp), dimension(nq):: q1,q2
        character(len=10),dimension(:), allocatable ::atoms
        real(kind=dp),dimension(:),allocatable :: x,y,z
        real(kind=dp),dimension(:,:),allocatable:: aff1,aff2
        real(kind=dp),dimension(:),allocatable ::ff
        character(len=10),dimension(151):: atoms_table
        real(kind=dp), dimension(4,151) :: a,b
        real(kind=dp),  dimension(151):: c



        call linspace(0.00001_dp,4.0_dp,nq,q1)
        call linspace(0.00001_dp,4.0_dp,nq,q2)

        open(unit=16,file="mol.xyz")
        read(16,*)natoms

        allocate(atoms(natoms), x(natoms), y(natoms), z(natoms),aff1(nq,natoms),aff2(nq,natoms),ff(nq))
        allocate(bond((natoms**2-natoms)/2))
        call table_of_ff("form_factors.dat",atoms_table,a,b,c)
        do i=1,natoms
            read(16,*)atoms(i),x(i),y(i),z(i)
            call obtain_form_factors(ff,atoms(i),q1,atoms_table,a,b,c)
            aff1(:,i)=ff
            call obtain_form_factors(ff,atoms(i),q2,atoms_table,a,b,c)
            aff2(:,i)=ff

        end do
        count=0
        do i=1,natoms
            do j=i+1,natoms
                bond(count)=
            end do
        end do
        close(16)








End Program Coherent_scattering




