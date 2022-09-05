module molden_molpro_reader
    implicit none

    contains

    subroutine read_data()
        implicit none
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        character*60, intent(in)::file_read
        integer(kind=ikind), dimension(:), allocatable, intent(out)::nat
        real(kind=dp), dimension(:), allocatable, intent(out):: x,y,z
        integer(kind=ikind):: i,j,k
        logical:: atoms,gtos

        atoms=.False.
        gtos=.true.


        open(file=file_read, unit=15)

        do while(gtos)

        end do






    end subroutine read_data


end module molden_molpro_reader