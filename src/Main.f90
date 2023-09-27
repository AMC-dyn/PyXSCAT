program main
    use omp_lib
    use main_calculation_mod
    use Reader
    use linspace
    use TSj0groups
    use TSj0contr
    implicit none
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
    character(len=60):: conf1, file_out,var,var2rdm
    character(len=100)::file_bit
    real(kind=dp), dimension(:), allocatable:: atoms
    integer(kind=ikind):: ngtos,norbs,nconfs,nstates,state1,state2,natoms,ncontr
    integer(kind=ikind),dimension(:,:), allocatable:: confs
    real(kind=dp),dimension(:), allocatable:: q
    real(kind=dp),dimension(:), allocatable:: coeffs
    real(kind=dp),dimension(:,:), allocatable:: civs
    real(kind=dp),dimension(:,:), allocatable:: mmod
    real(kind=dp),dimension(:,:), allocatable:: geom
    real(kind=dp),dimension(:), allocatable:: ga, xx,yy,zz
    real(kind=dp), dimension(:), allocatable:: result,q_abs
    integer*8, dimension(1):: nnn,start1,end1, nnn2, start_2, end_2,ordering1,ordering2
    integer*8, dimension(1,1):: newdat
    real(kind=dp):: cutoffcentre,cutoffz,cutoffmd,qmin,qmax
    integer(kind=ikind),dimension(:), allocatable:: l,m,n,group,gs,gf,gc
    integer(kind=ikind):: typec, i, j,k, npoints,ncivs,lconfs,maxl,ng,nq,count
    logical:: jeremyR, mcci, hf,molpro,molcas,bagel,bitwise,fci
     

    call read_files(nconfs,ngtos,norbs,ng,ncontr,state1,state2,natoms,typec,maxl,npoints, &
               cutoffcentre,cutoffz,cutoffmd,atoms,coeffs,xx,yy,zz,ga,l,m,n,group,mmod,geom, &
               jeremyR,mcci,hf,molpro,molcas,bagel,qmin,qmax,&
               file_out,file_bit,var,gs,gc,gf,confs,civs,lconfs,ncivs,bitwise)



    if (molpro.eqv..False. .and. molcas.eqv..False. .and. bagel.eqv..False.)  then
            !read(15,*) var
            print*,var
            var2rdm='readtwordm'
            print*,var2rdm

            if (var==var2rdm) then
              print*,'We have a 2RDM to read'

              !read(15,*) file_bit
              print*,file_bit
              open(unit=15, file=file_bit)
              count=0
              do while(.true.)
                  read (10, *, end=999) i
                  count=count+1
              end do


999           continue
              nconfs=count

              !close(15)

              call linspace_1(qmin,qmax,npoints,q)

              nq=npoints
!              open(unit=15,file='Exp_data/q_exp.dat')
!              do i=1,nq
!                 read(15,*) q(i)
!
!
!                  end do
!              close(15)
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
    print*,nq, npoints,maxl,ngtos
    nq=npoints
    if (bitwise.eqv..False.) then
    print*,'calling total_scattering'
    if (typec/=1 .and. typec/=11) then
    call total_scattering_calculation(Typec, atoms, geom, state1, &
            state2, maxl, Ngtos, ng,ga, l, m, n, xx, yy, zz, mmod, coeffs,q, nq, group, &
            ncontr,gs,gf,gc,cutoffz, cutoffmd, cutoffcentre, confs, civs,q_abs, result)
    elseif (typec==1) then
        call total_scattering_j0_groups(q, l, m, n, ngtos, ng, nq, maxl, typec, state1, &
            state2, ncontr, group, gs, gf, gc, confs, ga, xx, yy, zz, coeffs, mmod, civs, geom, &
            cutoffmd, cutoffz, cutoffcentre, result)
    elseif (typec==11) then
        call total_scattering_j0_ncontr(q, l, m, n, ngtos, ng, nq, maxl, typec, state1, &
            state2, ncontr, group, gs, gf, gc, confs, ga, xx, yy, zz, coeffs, mmod, civs, geom, &
            cutoffmd, cutoffz, cutoffcentre, result)

    end if
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

