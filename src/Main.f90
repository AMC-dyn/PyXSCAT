program main
    use omp_lib
    use main_calculation_mod
    use linspace
    implicit none
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
    character(len=60):: conf1, file_out,var,var2rdm
    character(len=100)::file_bit
    real(kind=dp), dimension(:), allocatable:: atoms
    integer(kind=ikind):: ngtos,norbs,nconfs,nstates,state1,state2,natoms
    integer(kind=ikind),dimension(:,:), allocatable:: confs
    real(kind=dp),dimension(:), allocatable:: q
    real(kind=dp),dimension(:,:), allocatable:: civs
    real(kind=dp),dimension(:,:), allocatable:: mmod
    real(kind=dp),dimension(:,:), allocatable:: geom
    real(kind=dp),dimension(:), allocatable:: ga, xx,yy,zz
    real(kind=dp), dimension(:), allocatable:: result,q_abs
    integer*8, dimension(1):: nnn,start1,end1, nnn2, start_2, end_2,ordering1,ordering2
    integer*8, dimension(1,1):: newdat
    real(kind=dp):: cutoffcentre,cutoffz,cutoffmd,qmin,qmax
    integer(kind=ikind),dimension(:), allocatable:: l,m,n,group
    integer(kind=ikind):: typec, i, j,k, npoints,ncivs,lconfs,maxl,ng,nq
    logical:: jeremyR, mcci, hf,molpro,molcas,bagel,bitwise,fci
     
     call OMP_set_num_threads(20) 
     print*,OMP_get_num_threads()    
    open(unit=15, file='basis.dat')
    read(15,*)ngtos
    allocate(xx(ngtos), yy(ngtos), zz(ngtos), ga(ngtos), l(ngtos), m(ngtos), n(ngtos), group(ngtos))
    do i=1,ngtos
        read(15,*)xx(i), yy(i),zz(i), ga(i), l(i), m(i), n(i), group(i)
    end do
    close(15)
    maxl=maxval(l)
    ng=maxval(group)

    open(unit=15, file='MOs.dat')
    read(15,*)norbs,ngtos
    allocate(mmod(norbs,ngtos))
    print*,norbs,ngtos
    do i=1,norbs
        read(15,*) (mmod(i,j), j=1, ngtos)
    end do
    close(15)

    open(unit=15, file='options.dat')
    read(15,*)natoms
    allocate(atoms(natoms))
    read(15,*)(atoms(i), i=1,natoms)
    allocate(geom(natoms,3))
    do i=1,natoms
        read(15,*)(geom(i,j), j=1,3)
    end do
    read(15,*)cutoffcentre
    read(15,*)cutoffz
    read(15,*)cutoffmd
    read(15,*)jeremyR
    print*,'jeremyR is ', jeremyR

    read(15,*)mcci
    read(15,*)hf
    read(15,*)qmin,qmax,npoints
    read(15,*)Typec
    read(15,*)state1,state2
    read(15,*)file_out
    read(15,*)molpro
    print*,'molpro is ', molpro
    if (molpro.eqv..True.) then
        read(15,*)nconfs
        read(15,*)lconfs
        read(15,*)ncivs
        if (nconfs<1E5) then

            allocate(confs(nconfs,lconfs),civs(nconfs,ncivs))

            do i=1,nconfs

                read(15,*)(confs(i,j),j=1,lconfs),(civs(i,j), j=1,ncivs)
            end do
        else
          bitwise=.True.
          file_bit='bitwise.dat'
        end if


    else
            read(15,*) var
            print*,var
            var2rdm='readtwordm'
            print*,var2rdm

            if (var==var2rdm) then
              print*,'We have a 2RDM to read'

              read(15,*) file_bit
              print*,file_bit
              close(15)

              call linspace_1(qmin,qmax,npoints,q)
              nq=npoints
              print*,nq, npoints
              newdat = 0
              start1 = 0
              end1 = 0
              start_2 = 0
              fci = .False.
              nnn = 0
              nnn2 = 0
              ordering1 = 0
              ordering2 = 0
              end_2 = 0
              nconfs=52236
              print*,mmod
              call  total_scattering_calculation_2(Typec, atoms, geom, state1, state2, maxl, &
      Ngtos, ng, ga, l, m, n, xx, yy, zz, mmod, q, nq, group, cutoffz, cutoffmd, cutoffcentre,file_bit, nconfs, newdat, &
      start1,end1, start_2, end_2, fci, nnn,nnn2, ordering1, ordering2, .True., .False.,q_abs,result)
              open(unit=15,file=file_out)

              do i=1,npoints
                  write(15,*)q(i),result(i)
              end do
              close(15)
              print*, file_out,' created'
              stop
            end if
        end if


    print*,'out of loop'
    close(15)
    call linspace_1(qmin,qmax,npoints,q)
    print*,nq, npoints
    nq=npoints
    if (bitwise.eqv..False.) then
    call total_scattering_calculation(Typec, atoms, geom, state1, &
            state2, maxl, Ngtos, ng,ga, l, m, n, xx, yy, zz, mmod, q, nq, group, &
            cutoffz, cutoffmd, cutoffcentre, confs, civs,q_abs, result)
    else
        newdat = 0
        start1 = 0
        end1 = 0
        start_2 = 0
        fci = .False.
        nnn = 0
        nnn2 = 0
        ordering1 = 0
        ordering2 = 0
        end_2 = 0
      call  total_scattering_calculation_2(Typec, atoms, geom, state1, state2, maxl, &
      Ngtos, ng, ga, l, m, n, xx, yy, zz, mmod, q, nq, group, cutoffz, cutoffmd, cutoffcentre,file_bit, nconfs, newdat, &
      start1,end1, start_2, end_2, fci, nnn,nnn2, ordering1, ordering2, .False., .False.,q_abs,result)


    end if

    open(unit=15,file=file_out)

    do i=1,npoints
        write(15,*)q(i),result(i)
    end do
    close(15)
    print*, file_out,' created'

end program main

