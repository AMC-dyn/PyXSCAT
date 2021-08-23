 function dyadic_11(T1, T2)
        implicit none

        type(Tensor1), intent(in) :: T1, T2
        type(Tensor2) :: dyadic_11
        integer i, j

        forall(i=1:3,j=1:3) dyadic_11%ab(i,j) = T1%a(i) * T2%a(j)

       end function dyadic_11
 