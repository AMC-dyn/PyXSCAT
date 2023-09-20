program readpython
    implicit none
    integer, parameter :: DPR = selected_real_kind(p=15)
    character(len=*), parameter :: filename = 'Zcotr.dat'
    real(dpr), allocatable:: rata2(:,:,:,:)
    allocate(rata2(44,44,44,44))
        open(40, file=filename, status='old', access='stream', form='unformatted')
        read(40) rata2
        close(40)

        print*,rata2(1,1:10,1,1)




end program readpython