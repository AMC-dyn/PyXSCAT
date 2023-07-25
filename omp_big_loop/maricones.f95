!-----------------------------------------------------------------------
! 2RDM Fortran Modules for PyXSCAT code
! Andres Moreno Carrascosa and Mats Simmermacher, (2021)
!-----------------------------------------------------------------------










MODULE twordmreader

    implicit none

    contains




    SUBROUTINE reduce_density(matst,totalst,m1,m2,m3,m4,total2)

        use types
        use membermodule
        use uniquemodule

        real(kind=dp), intent(in), dimension(:)                    :: totalst
        integer(kind=ikind), intent(in), dimension(:,:) :: matst
        integer(kind=ikind), dimension(:), intent(out), allocatable      :: m1, m2, m3, m4
        real(kind=dp), dimension(:), intent(out),allocatable             :: total2
        integer(kind=ikind)                                              :: i, j, k, l, cc
        integer(kind=ikind)                                              :: n, sdr
        integer(kind=ikind)                                              :: cnt, num
        logical                                                          :: memb

        integer(kind=ikind), dimension(:,:), allocatable                 :: newmat
        real(kind=dp), dimension(:), allocatable      :: stotal,newtotal
        integer(kind=ikind), dimension(:,:), allocatable    :: smat
        integer(kind=ikind)                                              :: mo1, mo2, mo3, mo4
        integer(kind=ikind), dimension(4)                                :: b

        integer(kind=ikind), dimension(:,:), allocatable                 :: mat2

        integer(kind=ikind), dimension(:), allocatable :: num1




        call unique_total(matst, totalst, mat2, total2 )
        
        sdr = size(mat2(:,1))

        allocate(m1(sdr), m2(sdr), m3(sdr), m4(sdr))


        m1 = mat2(:,1)
        m2 = mat2(:,2)
        m3 = mat2(:,3)
        m4 = mat2(:,4)
        


   
    END SUBROUTINE reduce_density



END MODULE twordmreader


