program trying
    interface

    subroutine TwoRDM_SD(length)
        integer :: length
        end subroutine TwoRDM_SD
    end interface
    call TwoRDM_SD(149850)
end program trying



subroutine TwoRDM_SD(length) !length of wavefunction, energy from 2RDM (without nuclear contribution)
 ! basis functions, electrons, coefficients, list of orbitals in SD
implicit none
integer ::  length
double precision, allocatable :: SpinFree2RDM(:,:,:,:)
integer porb,rorb,sorb,qorb
integer spins,spinq,spinp,spinr
double precision dtemp,ep,eg
integer ici,jci,ndiff,idiff1,idiff2
integer i,i1,i2,l,l2,k,k2
double precision Eone_e,Etwo_e,TwoRDM_e
integer newdiff,mytemp,n
  integer hole,part,tz,tz2,buffer,nword
  integer myhigh,mylow,nperm
  integer idx_part,idx_hole
  integer exc(0:2,2,2),mya,myb,myc,myd
  integer myspin,samespin,ispin,nbft
  integer buffer2,ispin2, dummy,length2
  double precision, parameter :: phase_dbl(0:1)=(/1.d0,-1.d0/)
  double precision c_term
  double precision :: c(length)
  integer icij(2,1,length)

open(15,file='O3/civ_out_0.0001')
nword=1
 do i=1,length
     read(15,*)dummy, c(i), icij(1,1,i), icij(2,1,i)
 enddo
close(15)


 nbft=27

print*,length 
length2=length

ALLOCATE(SpinFree2RDM(nbft,nbft,nbft,nbft))

print*,'allocated matrix'
  !set to zero (a while ago Array=0.0D0 did not work for some compilers)
 do porb=1,nbft
 do rorb=1,nbft
 do sorb=1,nbft
 do qorb=1,nbft
 SpinFree2RDM(porb,rorb,sorb,qorb)=0.0D0
 end do
 end do
 end do
 end do

print*, 'spinfreee filled'
dtemp=0.0D0
!ici=jci first





do ici=1,length2


!only zero differences by construction  ! could try to make further improvements but there are only length terms for no differences not O(length**2) as for 1 and 2
c_term=c(ici)**2 !calculate once as  common to all zero differences for this ici
do ispin=1,2
buffer=icij(ispin,1,ici)

do while(buffer.ne.0)
sorb=trailz(buffer)+1
buffer=IAND(buffer,buffer-1)

do ispin2=1,2
buffer2=icij(ispin2,1,ici)
do while(buffer2.ne.0)
qorb=trailz(buffer2)+1
buffer2=IAND(buffer2,buffer2-1)

if((qorb.eq.sorb).AND.(ispin.eq.ispin2)) cycle ! all possible choices of two so can't pick the same twice

  !Case 1 p is same as s and r is q
  if(ispin.eq.ispin2) THEN

SpinFree2RDM(sorb,qorb,sorb,qorb)=SpinFree2RDM(sorb,qorb,sorb,qorb)-c_term !porb=sorb rorb=qorb
  END IF

!Case 2 p is same as q and r is same as s
!All s and q spins valid for this contribution

SpinFree2RDM(qorb,sorb,sorb,qorb)=SpinFree2RDM(qorb,sorb,sorb,qorb)+c_term !p=q r=s

end do
end do
end do
end do ! end of zero differences


end do !end of ici loop

print*,'starting big loop'

!now jci>ici and double values so we don't need to do jci<ici
do ici=1,length2-1
do jci=ici+1,length2


!!!!!!!!!!!!!
newdiff=0
do n=1,nword
mytemp=IEOR(icij(1,n,ici),icij(1,n,jci))

newdiff=newdiff+POPCNT(mytemp) ! calcs number of bits set to 1 as we used xor (IEOR) bitwise that must be - note one difference adds two to newdiff as there are two places where the bits will be set to 1

mytemp=IEOR(icij(2,n,ici),icij(2,n,jci))

newdiff=newdiff+POPCNT(mytemp)

end do
 if (newdiff.gt.4) cycle !more than two differences so matrix element is zero



        if(newdiff.eq.4) THEN ! get differences and phase

exc(0,1,1)=0
exc(0,2,1) =0
exc(0,1,2)=0
exc(0,2,2)=0
nperm=0
samespin=0
do myspin=1,2
idx_part=0
idx_hole=0
mytemp=(IEOR(icij(myspin,1,ici),icij(myspin,1,jci)))
hole= IAND(mytemp,icij(myspin,1,ici))
part= IAND(mytemp,icij(myspin,1,jci))

do while(part.ne.0)
tz=trailz(part)
idx_part=idx_part+1
exc(0,2,myspin)=exc(0,2,myspin)+1
exc(idx_part,2,myspin)=tz+1
part=iand(part,part-1)
end do

do while(hole.ne.0)
tz=trailz(hole)
idx_hole=idx_hole+1
exc(0,1,myspin)=exc(0,1,myspin)+1
exc(idx_hole,1,myspin)=tz+1
hole=iand(hole,hole-1)
end do

do i=1,exc(0,1,myspin)
mylow=min(exc(i,1,myspin),exc(i,2,myspin))
myhigh=max(exc(i,1,myspin),exc(i,2,myspin))
nperm=nperm+POPCNT(IAND(icij(myspin,1,ici),&
IAND(ibset(0,myhigh-1)-1,ibclr(-1,mylow)+1)))


end do
!both holes have same spin
if(exc(0,1,myspin).eq.2) THEN !this seems less efficient...
samespin=myspin
mya=min(exc(1,1,myspin),exc(1,2,myspin))
myb=max(exc(1,1,myspin),exc(1,2,myspin))
myc=min(exc(2,1,myspin),exc(2,2,myspin))
myd=max(exc(2,1,myspin),exc(2,2,myspin))

if((myc>mya).AND.(myc<myb).AND.(myd>myb)) nperm=nperm+1
exit
END IF ! end of check if both are same spin

end do !loop over myspin

  c_term=c(ici)*c(jci)*phase_dbl(iand(nperm,1)) ! calculate once as common to all in this loop

   if(samespin.eq.0) THEN

      !Case 1
   sorb=exc(1,2,1) !korb
     qorb=exc(1,2,2)!lorb
     porb=exc(1,1,2) !jorb
     rorb=exc(1,1,1) !iorb
 ! print*,ici,jci,c_term
 ! print*, porb,rorb,sorb,qorb
      !spins.eq.spinr and spinq.eq.spinp by construction
SpinFree2RDM(porb,rorb,sorb,qorb)=SpinFree2RDM(porb,rorb,sorb,qorb)&
+c_term

SpinFree2RDM(qorb,sorb,rorb,porb)=SpinFree2RDM(qorb,sorb,rorb,porb)&
+c_term
   !Case 2 Swap p and r introduces negative sign but means spins.ne.spinr so no contribution
   !Case 3 from Case 1 swap s and q to give negative sign but spins.ne.spinr so no contribution
   !Case 4 from Case 1 swap s and q then swap p and r so no sign change

         !spins.eq.spinr and spinq.eq.spinp by construction
SpinFree2RDM(rorb,porb,qorb,sorb)=SpinFree2RDM(rorb,porb,qorb,sorb)&
+c_term


SpinFree2RDM(sorb,qorb,porb,rorb)=SpinFree2RDM(sorb,qorb,porb,rorb)&
+c_term
    ELSE
  !samespin.eq.1 or   2
  !Case 1
     sorb=exc(1,2,samespin) !korb
     qorb=exc(2,2,samespin)!lorb
     porb=exc(2,1,samespin) !jorb
     rorb=exc(1,1,samespin) !iorb
SpinFree2RDM(porb,rorb,sorb,qorb)=SpinFree2RDM(porb,rorb,sorb,qorb)&
+c_term

SpinFree2RDM(qorb,sorb,rorb,porb)=SpinFree2RDM(qorb,sorb,rorb,porb)&
+c_term
!all same spin so all swaps are allowed
   !Case 2 Swap p and r introduces negative sign

SpinFree2RDM(rorb,porb,sorb,qorb)=SpinFree2RDM(rorb,porb,sorb,qorb)&
-c_term

SpinFree2RDM(sorb,qorb,rorb,porb)=SpinFree2RDM(sorb,qorb,rorb,porb)&
-c_term
   !Case 3 from Case 1 swap s and q to give negative sign
SpinFree2RDM(porb,rorb,qorb,sorb)=SpinFree2RDM(porb,rorb,qorb,sorb)&
-c_term

SpinFree2RDM(qorb,sorb,porb,rorb)=SpinFree2RDM(qorb,sorb,porb,rorb)&
-c_term
   !Case 4 from Case 1 swap s and q then swap p and r so no sign change

SpinFree2RDM(rorb,porb,qorb,sorb)=SpinFree2RDM(rorb,porb,qorb,sorb)&
+c_term

       SpinFree2RDM(sorb,qorb,porb,rorb)=SpinFree2RDM(sorb,qorb,porb,rorb)&
+c_term


   end if ! end of samespin check



cycle


end if ! end of two differences (newdiff.eq.4)

!one difference
if(newdiff.eq.2) THEN !  newdiff.eq.2 is one difference

do myspin=1,2
mytemp=(IEOR(icij(myspin,1,ici),icij(myspin,1,jci)))
hole= IAND(mytemp,icij(myspin,1,ici))
part= IAND(mytemp,icij(myspin,1,jci))
if(hole.ne.0) THEN
  tz=1+trailz(hole)
  tz2=1+trailz(part)
  EXIT
  END IF
  end do


!get sign for single
!myspin shows if alpha or beta at this point
mylow=min(tz,tz2)
myhigh=max(tz,tz2)
!PRINT *,mylow,myhigh
nperm=POPCNT(IAND(icij(myspin,1,ici),IAND(ibset(0,myhigh-1)-1,ibclr(-1,mylow)+1)))



  c_term=c(ici)*c(jci)*phase_dbl(iand(nperm,1)) ! calculate once as common to all in this loop

!case 1 and 4 and 2 and 3
porb=tz
qorb=tz2

!ispin=myspin
buffer=IBCLR(icij(myspin,1,ici),porb-1) ! remove porb of this spin so rorb cannot be porb


do while(buffer.ne.0)
rorb=trailz(buffer)+1
buffer=IAND(buffer,buffer-1)


!case 1
SpinFree2RDM(porb,rorb,rorb,qorb)=SpinFree2RDM(porb,rorb,rorb,qorb)+c_term  !sorb=rorb
!case 4 s and r are the differences so signs of moving through occupied will cancel
SpinFree2RDM(rorb,porb,qorb,rorb)=SpinFree2RDM(rorb,porb,qorb,rorb)+c_term !sorb=rorb
!case 2 spins and spinr are the same by construction as are spinp and spinq
SpinFree2RDM(rorb,porb,rorb,qorb)=SpinFree2RDM(rorb,porb,rorb,qorb)-c_term !sorb=rorb
! case 3 s and p are the different ones - sign change so both operators to come before occupied when acting on their det (left det for p and r, right det for s and q)
SpinFree2RDM(porb,rorb,qorb,rorb)=SpinFree2RDM(porb,rorb,qorb,rorb)-c_term !sorb=rorb


SpinFree2RDM(qorb,rorb,rorb,porb)=SpinFree2RDM(qorb,rorb,rorb,porb)+c_term  !sorb=rorb
!case 4 s and r are the differences so signs of moving through occupied will cancel
SpinFree2RDM(rorb,qorb,porb,rorb)=SpinFree2RDM(rorb,qorb,porb,rorb)+c_term !sorb=rorb
!case 2 spins and spinr are the same by construction as are spinp and spinq
SpinFree2RDM(rorb,qorb,rorb,porb)=SpinFree2RDM(rorb,qorb,rorb,porb)-c_term !sorb=rorb
! case 3 s and p are the different ones - sign change so both operators to come before occupied when acting on their det (left det for p and r, right det for s and q)
SpinFree2RDM(qorb,rorb,porb,rorb)=SpinFree2RDM(qorb,rorb,porb,rorb)-c_term !sorb=rorb


!  if (rorb==1 .and. porb==4 .and. qorb==3) then
!      print*, SpinFree2RDM(rorb,porb,qorb,rorb),ici,jci
!  end if

end do

buffer=icij(IEOR(myspin,3),1,ici)  !map spin 1 to spin 2 and 2 to 1

do while(buffer.ne.0)
rorb=trailz(buffer)+1
buffer=IAND(buffer,buffer-1)
!only have case 1 and 2 as spins are different
!case 1
SpinFree2RDM(porb,rorb,rorb,qorb)=SpinFree2RDM(porb,rorb,rorb,qorb)+c_term  !sorb=rorb
!case 4 s and r are the differences so signs of moving through occupied will cancel
SpinFree2RDM(rorb,porb,qorb,rorb)=SpinFree2RDM(rorb,porb,qorb,rorb)+c_term !sorb=rorb

SpinFree2RDM(qorb,rorb,rorb,porb)=SpinFree2RDM(qorb,rorb,rorb,porb)+c_term  !sorb=rorb
!case 4 s and r are the differences so signs of moving through occupied will cancel
SpinFree2RDM(rorb,qorb,porb,rorb)=SpinFree2RDM(rorb,qorb,porb,rorb)+c_term !sorb=rorb

end do

!end of case 1 and case 4 and 2 and 3




cycle
end if  ! end of  single difference, newdiff.eq.2





end do !end of ici loop
end do !end of jci loop

print*,'here we are'



!Write out 2RDM

OPEN(UNIT=15,FILE='NonZero2RDM_MO_Coeffs')
do l=1,nbft
do l2=1,nbft
do k=1,nbft
do k2=1,nbft


if (ABS(SpinFree2RDM(l,l2,k,k2)).eq.0.0D0) cycle
!WRITE(15,*) l,l2,k,k2,SpinFree2RDM(l,l2,k,k2)

WRITE(15,*) l,k2,l2,k,SpinFree2RDM(l,l2,k,k2) ! notation from arxiv 1809.09058
end do
end do
end do
end do
CLOSE(15)

end subroutine

