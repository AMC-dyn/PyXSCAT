! Created by  on 26/09/2023.

module Reader

    implicit none
    contains

    subroutine read_files(nconfs,ngtos,norbs,ng,ncontr,state1,state2,natoms,typec,maxl,npoints, &
               cutoffcentre,cutoffz,cutoffmd,atoms,coeffs,xx,yy,zz,ga,l,m,n,group,mmod,geom, &
               jeremyR,mcci,hf,molpro,molcas,bagel,qmin,qmax,&
               file_out,file_bit,var,gs,gc,gf,confs,civs,lconfs,ncivs,bitwise)
    implicit none

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
    integer(kind=ikind),intent(out):: ngtos,norbs,nconfs,state1,typec
    integer(kind=ikind),intent(out):: state2,natoms,ncontr,maxl,npoints,ng
    integer(kind=ikind),intent(out):: lconfs,ncivs
    real(kind=dp), intent(out) :: cutoffcentre,cutoffz,cutoffmd, qmin,qmax
    real(kind=dp),dimension(:), allocatable,intent(out)::atoms,coeffs,ga,xx,yy,zz
    integer(kind=ikind), dimension(:), allocatable, intent(out):: l,m,n,group,gs,gc, gf
    integer(kind=ikind),dimension(:,:), allocatable,intent(out):: confs
    real(kind=dp),dimension(:,:),intent(out), allocatable:: civs,mmod,geom
    logical, intent(out) :: jeremyR,mcci,hf,molpro,molcas,bagel,bitwise
    character(len=60), intent(out):: file_out,file_bit,var
    character(len=60) :: var2rdm
    integer(kind=ikind)::i,j




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
     open(unit=15, file='coeffs.dat')
     read(15,*)
     allocate(coeffs(ngtos))
     do i=1,ngtos
      read(15,*) coeffs(i)
    end do
     read(15,*)ncontr
     allocate(gs(ncontr), gf(ncontr), gc(ncontr))
     do i=1,ncontr
         read(15,*) gs(i), gf(i),gc(i)
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
  read(15,*)molcas
  read(15,*)bagel
  print*,'molpro is ', molpro
  print*,'molcas is ', molcas
  print*,'bagel is ', bagel
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

  else if (molcas.eqv..True.) then
    read(15,*)nconfs
    read(15,*)lconfs
    read(15,*)ncivs
    allocate(confs(nconfs,lconfs),civs(nconfs,ncivs))

    do i=1,nconfs

    read(15,*)(confs(i,j),j=1,lconfs),(civs(i,j), j=1,ncivs)
    end do
  else if (bagel.eqv..True.) then
    read(15,*)nconfs
    read(15,*)lconfs
    read(15,*)ncivs
    allocate(confs(nconfs,lconfs),civs(nconfs,ncivs))

    do i=1,nconfs

    read(15,*)(confs(i,j),j=1,lconfs),(civs(i,j), j=1,ncivs)
    end do
    else
            read(15,*) var

            var2rdm='readtwordm'


            if (var==var2rdm) then
              print*,'We have a 2RDM to read'

              read(15,*) file_bit
              print*,file_bit
              close(15)

              end if
              end if
    end subroutine

end module Reader