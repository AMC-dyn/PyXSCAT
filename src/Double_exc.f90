program double_exc
    implicit none
    integer*8 nel, orbs, initial,fullocc,diff
    integer*8 i,j,k,l,m,n
    integer*8 fin_number,exit_flag,exit_flag_2
    integer*8 buffer,buffer_2,buffer_3,buffer_4,count,len_rm
    integer*8, dimension(:), allocatable :: occ,unocc

    nel=3
    count=0
    initial=11
    orbs=6
    fullocc=0
     allocate(occ(nel),unocc(orbs-nel))

    count=1
    buffer=initial
    do while(buffer/=0_8)
        i=trailz(buffer)
        occ(count)=i
        count=count+1
        buffer=iand(buffer,buffer-1_8)
    end do


    do i=1,orbs
        fullocc=fullocc+2**(i-1)
    end do

    diff=xor(initial,fullocc)

    buffer=diff
    count=1
    do while(buffer/=0_8)
        i=trailz(buffer)
        unocc(count)=i
        count=count+1
        buffer=iand(buffer,buffer-1_8)
    end do
print*,fullocc
    count=0
   do i=1,size(occ)-1
       do j=i+1,size(occ)
           do k=1,size(unocc)-1
               do l=k+1,size(unocc)
                  buffer=ibclr(initial,occ(i))

                  buffer=ibclr(buffer,occ(j))

                  buffer=ibset(buffer,unocc(k))

                  buffer=ibset(buffer,unocc(l))

                  print*,buffer,count

                  count=count+1



               end do
           end do
       end do
    end do












end program double_exc