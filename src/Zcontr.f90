
Program Zcontr
    use Reader
    use onerdm


    implicit none
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
    character(len=60):: conf1, file_out,var,var2rdm,file_in
    character(len=100)::file_bit
    real(kind=dp), dimension(:), allocatable:: atoms
    integer(kind=ikind):: ngtos,norbs,nconfs,nstates,state1,state2,natoms,ncontr,numberlines,norbs2,ngtos2
    integer(kind=ikind),dimension(:,:), allocatable:: confs

    real(kind=dp),dimension(:), allocatable:: q
    real(kind=dp),dimension(:), allocatable:: coeffs,total
    real(kind=dp),dimension(:,:), allocatable:: civs
    real(kind=dp),dimension(:,:), allocatable:: mmod
    real(kind=dp),dimension(:,:), allocatable:: geom
    real(kind=dp),dimension(:), allocatable:: ga, xx,yy,zz
    real(kind=dp), dimension(:), allocatable:: result,q_abs
    integer*8, dimension(1):: nnn,start1,end1, nnn2, start_2, end_2,ordering1,ordering2
    integer*8, dimension(1,1):: newdat
    real(kind=dp):: cutoffcentre,cutoffz,cutoffmd,qmin,qmax
    integer(kind=ikind),dimension(:), allocatable:: l,m,n,group,gs,gf,gc,contrvec
    integer(kind=ikind):: typec, i, j,k, npoints,ncivs,lconfs,maxl,ng,nq,count
    logical:: jeremyR, mcci, hf,molpro,molcas,bagel,bitwise,fci
    integer(kind = ikind), dimension(:, :), allocatable :: mat
    real(kind=dp), dimension(:,:,:,:), allocatable :: twordm,zcontract,prueba1
    real(kind=dp), dimension(:,:), allocatable:: interm,interm2,mos



    call read_files(nconfs,ngtos,norbs,ng,ncontr,state1,state2,natoms,typec,maxl,npoints, &
               cutoffcentre,cutoffz,cutoffmd,atoms,coeffs,xx,yy,zz,ga,l,m,n,group,mmod,geom, &
               jeremyR,mcci,hf,molpro,molcas,bagel,qmin,qmax,&
               file_out,file_bit,var,gs,gc,gf,confs,civs,lconfs,ncivs,contrvec,bitwise)

            open(unit=15, file='bitwise.dat')
              numberlines=0
              do while(.true.)
                  read (15, *, end=999) i
                  numberlines=numberlines+1
              end do


999 continue
        close(15)

    file_in='bitwise.dat'
    file_out='es.dat'
    call mcci_to_bit(file_in,file_out,numberlines)

    call createtwordm_bit(file_out,numberlines,mat,total)
    norbs=maxval(mat(:,1))
    allocate(twordm(norbs,norbs,norbs,norbs), zcontract(ncontr,ncontr,ncontr,ncontr))
    open(unit=15, file='MOs2.dat')
    read(15,*)norbs2,ngtos2
    allocate(mos(norbs2,ngtos2))

    do i=1,norbs2
        read(15,*) (mos(i,j), j=1, ngtos2)
    end do
    close(15)
    print*,shape(mos),norbs
    do i=1,size(total)


        twordm(mat(i,1), mat(i,2), mat(i,3),mat(i,4))=total(i)

    end do


    allocate(interm(norbs**3,norbs),interm2(norbs**3,ncontr))
    interm=reshape(twordm, (/norbs**3,norbs/))
    print*,shape(interm)

    call dgemm('n', 't', norbs**3, ncontr, norbs, 1.0_dp ,interm , &
                                        &           norbs**3,mos , ncontr, 0.0_dp, interm2, norbs**3)


    allocate(prueba1(norbs,norbs,norbs,ncontr))
     prueba1=reshape(interm2,(/norbs,norbs,norbs,ncontr/),order=[2,1])
    print*,prueba1(3,1,1,1)
    stop
    deallocate(interm)
    allocate(interm(norbs**2*ncontr,norbs))
    interm=reshape(interm2, (/norbs**2*ncontr,norbs/))
    deallocate(interm2)
    allocate(interm2(norbs**2*ncontr,ncontr))
    print*,shape(interm), shape(interm2)

    call dgemm('n', 't', norbs**2*ncontr, ncontr, norbs, 1.0_dp ,interm , &
                                        &           ncontr*norbs**2,mos , ncontr, 0.0_dp, interm2, norbs**2*ncontr)




    print*,shape(interm2)

    deallocate(interm)
    allocate(interm(norbs*ncontr**2,norbs))
    interm=reshape(interm2, (/norbs*ncontr**2,norbs/))
    deallocate(interm2)
    allocate(interm2(norbs*ncontr**2,ncontr))
    print*,shape(interm), shape(interm2)

    call dgemm('n', 't', norbs*ncontr**2, ncontr, norbs, 1.0_dp ,interm , &
                                        &           ncontr**2*norbs,mos , ncontr, 0.0_dp, interm2, norbs*ncontr**2)

    print*,shape(interm2)
      print*,shape(interm2)

    deallocate(interm)
    allocate(interm(ncontr**3,norbs))
    interm=reshape(interm2, (/ncontr**3,norbs/))
    deallocate(interm2)
    allocate(interm2(ncontr**3,ncontr))
    print*,shape(interm), shape(interm2)

    call dgemm('n', 't', ncontr**3, ncontr, norbs, 1.0_dp ,interm , &
                                        &           ncontr**3,mos , ncontr, 0.0_dp, interm2, ncontr**3)

    print*,shape(interm2)
    zcontract=reshape(interm2,(/ncontr,ncontr,ncontr,ncontr/))
    print*,zcontract(1,1,1,1)
End Program Zcontr