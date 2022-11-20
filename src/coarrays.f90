program cob
implicit none
integer:: x[*]
integer,allocatable:: a(:)[*]
integer:: i,val,ncap,newcap

print*, this_image(),num_images()

sync all

    x[this_image()]=this_image()
    call multiplicando(x,a)



end program cob

subroutine multiplicando(x,a)
    implicit none
    integer, intent(in)::x
    integer,intent(out)::a

    a=x

end subroutine