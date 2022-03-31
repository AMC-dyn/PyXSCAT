!-----------------------------------------------------------------------
! Fortran Modules for PyXSCAT code
! Andres Moreno Carrascosa and Mats Simmermacher, 2021
!-----------------------------------------------------------------------

module main_calculation_mod


    implicit none 

    contains

subroutine total_scattering_calculation(type,Zn,geom,state1,state2,maxl,ngto,ng,ga,l,m,n,xx,yy,zz, &
        mmod,q,nq, group,&
        cutoffz,cutoffmd,cutoffcentre,confs,civecs,q_abs,result)




    use onerdm

    use twordms

    use variables

    use integrals

    implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)

        integer(kind=ikind), intent(in):: ngto, ng,  nq, maxl,type, state1,state2
        integer(kind=ikind), intent(in),dimension(:) :: l, m, n,group
        integer(kind=ikind), dimension(:,:), intent(in):: confs

        real(kind=dp), intent(in),dimension(:) :: ga, xx, yy, zz, q, Zn
        real(kind=dp), intent(in),dimension(:,:) :: mmod, civecs,geom
        real(kind=dp),dimension(:,:,:),allocatable :: z1,z2
        real(kind=dp),dimension(ngto,ngto):: z

        real(kind=dp), intent(in) :: cutoffmd, cutoffz,cutoffcentre

        real(kind=dp), intent(out), dimension(nq):: result,q_abs
        REAL(kind=dp), DIMENSION(size(q),4*maxval(l)+1,4*maxval(l)+1,4*maxval(l)+1) :: P0matrix
         real(kind=dp),  dimension(maxl*2+1,maxl+1,maxl+1,ng,ng) :: ddx,ddy,ddz
        real(kind=dp), dimension(ng,ng) :: px,py,pz
        real(kind=dp),  dimension(:,:,:,:), allocatable :: zcontr
        real(kind=dp),  dimension(:,:), allocatable :: onerdm_matrix
        real(kind=dp),  dimension(:), allocatable :: total,newtotal
        real(kind=dp),  dimension(nq,ngto,ngto) :: e12
        INTEGER(kind=ikind), DIMENSION(maxval(group))   :: group_start, group_count
        integer(kind=ikind), dimension(:), allocatable :: m1, m2, m3, m4
        integer(kind=ikind), dimension(:,:), allocatable :: mat,ep3,ndiff2
        integer(kind=ikind):: nmat,i,j,nmomax,LL,MM,NN,k
        real(kind=dp) :: start,time1,time2,time3,time4,co,wl,rr,k0,rij
        complex(kind=dp), dimension(:,:,:,:), allocatable :: exponent1, exponent2
        complex(kind=dp), dimension(size(q)):: resultaligned
        real(kind=dp),dimension(3,size(q)) :: q_al
        real(KIND=dp), dimension(size(q)) :: ss,ss2,Iee,Inn,Ine

        print*, 'here we are '

       ! allocate(ep3(size(confs(:,1)),size(confs(:,1))),ndiff2(size(confs(:,1)),size(confs(:,1))) )

        do i = 1, Ng
            do j = 1, Ngto
                if (group(j) == i) then
                    group_start(i) = j
                    group_count(i) = count(group==i)
                    exit
                end if
            end do
        end do
        q_abs=0
       ! call onerdm_creat(confs,civecs,onerdm_matrix,nmomax)
        if (type==1) then
            call cpu_time(time2)
            P0matrix = 0
            CALL set_P0(P0matrix, 4*maxval(l), q)
            write(*,*)
            call maxcoincidence(confs,ep3,ndiff2)
            call createtwordm(confs,civecs,ndiff2,ep3,mat,total,state1,state2)
            allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
             m1 = mat(:,1)
             m2 = mat(:,2)
             m3 = mat(:,3)
             m4 = mat(:,4)
             nmat=size(m1)


            call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq)

            call cpu_time(time3)
            print*,'Time variables', time3-time2

            call tot_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z1,z2,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result)
            print*,'here'
        else if (type==2) then
               P0matrix = 0
               CALL set_P0(P0matrix, 4*maxval(l), q)
            call maxcoincidence(confs,ep3,ndiff2)
            call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)
               print*,'maxnmo',nmomax
               co=0
               do i=1,nmomax
                   co=co+onerdm_matrix(i,i)
               end do
               print*,'trace',co
            call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,onerdm_matrix,nmomax,q,nq)
            call elastic_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result)


        else if (type==3) then
            call cpu_time(time2)


            call maxcoincidence(confs,ep3,ndiff2)
            call createtwordm(confs,civecs,ndiff2,ep3,mat,total,state1,state2)
            allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
             m1 = mat(:,1)
             m2 = mat(:,2)
             m3 = mat(:,3)
             m4 = mat(:,4)
             nmat=size(m1)

              allocate(exponent1(nq,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1), &
                    exponent2(nq,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1))

             wl=4.0_dp*dacos(-1.0_dp)/maxval(q)
             wl=3.00
               call linspace(0.0_dp,dacos(-1.0_dp)/2,nq,ss);
            ss=(ss-dacos(-1.0_dp))/2.0_dp
            print*,'linspace calculated', wl
            do mm=1,nq
                k0 = wl*1.889726125*0.529
                rr=2*k0*cos(ss(mm))
                q_al(1,mm)=rr*sin(ss(mm))*cos(0.0_dp)
                q_al(2,mm)=rr*sin(ss(mm))*sin(0.0_dp)
                q_al(3,mm)=rr*cos(ss(mm))

                q_abs(mm)=sqrt(q_al(1,mm)**2.0_dp+q_al(2,mm)**2.0_dp+q_al(3,mm)**2.0_dp)
                end do
            call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q_abs,nq)



            print*,'q_al calculated', q_al(1,:)


            do LL=0,maxval(l)*4
                do MM=0,maxval(m)*4
                    do NN=0,maxval(n)*4

                    exponent2(:,LL+1,MM+1,NN+1)=((0,-1.00)*q_al(1,:))**(LL)* &
                            ((0,-1.00)*q_al(2,:))**(MM)*((0,-1.00)*q_al(3,:))**(NN)

                    exponent1(:,LL+1,MM+1,NN+1)=((0,1.00)*q_al(1,:))**(LL) &
                            *((0,1.00)*q_al(2,:))**(MM)*((0,1.00)*q_al(3,:))**(NN)
                     enddo
                enddo
            enddo
            print*,'exponents calculated', (exponent1(1,1,1,1))

            call tot_integration_aligned(ng,px,py,pz,l,m,n,ddx,ddy,ddz,z1,z2,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q_al,e12,exponent1,exponent2,resultaligned)

            result=real(resultaligned)



         else if (type==4) then


             call maxcoincidence(confs,ep3,ndiff2)
            call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)
               print*,'maxnmo',nmomax
             wl=4.0_dp*dacos(-1.0_dp)/maxval(q)
             wl=3.00
             call linspace(0.0_dp,dacos(-1.0_dp)/2,nq,ss);
            ss=(ss-dacos(-1.0_dp))/2.0_dp
            print*,'linspace calculated', wl
            do mm=1,nq
                k0 = wl*1.889726125*0.529
                rr=2*k0*cos(ss(mm))
                q_al(1,mm)=rr*sin(ss(mm))*cos(0.0_dp)
                q_al(2,mm)=rr*sin(ss(mm))*sin(0.0_dp)
                q_al(3,mm)=rr*cos(ss(mm))

                q_abs(mm)=sqrt(q_al(1,mm)**2.0_dp+q_al(2,mm)**2.0_dp+q_al(3,mm)**2.0_dp)
            end do
             co=0
             print*,'qvector',size(q_abs),maxval(q_abs),maxval(ss)
               do i=1,nmomax
                   co=co+onerdm_matrix(i,i)
               end do
               print*,'trace',co
             call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,onerdm_matrix,nmomax,q_abs,nq)
            print*, 'variables calculated'
            allocate(exponent1(nq,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1))
             do LL=0,maxval(l)*4
                do MM=0,maxval(m)*4
                    do NN=0,maxval(n)*4


                    exponent1(:,LL+1,MM+1,NN+1)=((0,1.00)*q_al(1,:))**(LL) &
                            *((0,1.00)*q_al(2,:))**(MM)*((0,1.00)*q_al(3,:))**(NN)
                     enddo
                enddo
            enddo
            print*, 'calling elastic'
             call elastic_integration_alig(ng,px,py,pz,l,m,n,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q_al,e12,exponent1,resultaligned)
                 result=abs(resultaligned)**2



            !Starting the rotationally averaged elastic electron scattering, in the end 4 different versions must be programmed

            else if (type==5) then

                P0matrix = 0
                CALL set_P0(P0matrix, 4*maxval(l), q)
                call maxcoincidence(confs,ep3,ndiff2)
                call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)

                call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,onerdm_matrix,nmomax,q,nq)

                print*,'first Z', Z(1,1)

                call elastic_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,Iee)



                print*,Iee(1),geom

                call nuclei_electron_integration(Zn,geom,ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,Ine)

                print*,Ine(1)

                Inn=0.0_dp
                do i=1,size(Zn)
                    do j=1,size(Zn)
                        do k=1,nq
                         Rij=sqrt(sum((geom(i,:)-geom(j,:))**2) )
                         Inn(k)=Inn(k)+Zn(i)*Zn(j)* sinc(q(k)*abs(Rij))

                         end do
                    enddo
                end do

                print*,Inn(1)

                print*,sum(Ine-Iee)
                result=Iee-2*Ine+Inn






                !we need three terms here, Inn, Iee and Ine



        end if
        end subroutine total_scattering_calculation


        subroutine total_scattering_calculation_2(type,Zn,geom,state1,state2,maxl,ngto,ng,ga,l,m,n,xx,yy,zz, &
        mmod,q,nq, group,&
        cutoffz,cutoffmd,cutoffcentre,fileJ,numberlines,newdat,start,end,start_2, end_2, fci,irep1,irep2,ordering1 &
                ,ordering2, read2rdm,mcci,q_abs,result)




    use onerdm

    use twordms

    use variables

    use integrals

    implicit none

        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
        INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
        logical, intent(in) :: fci,read2rdm,mcci
        integer(kind=ikind), intent(in):: ngto, ng,  nq, maxl,type, state1,state2,numberlines
        integer(kind=ikind), intent(in),dimension(:) :: l, m, n,group
        integer*8, intent(in), dimension(:,:):: newdat
        integer*8,intent(in), dimension(:):: irep1,start,end, irep2, start_2, end_2,ordering1,ordering2


        real(kind=dp), intent(in),dimension(:) :: ga, xx, yy, zz, q, Zn
        real(kind=dp), intent(in),dimension(:,:) :: mmod, geom
        real(kind=dp),dimension(:,:,:),allocatable :: z1,z2
        real(kind=dp),dimension(ngto,ngto):: z
        character (len=60), intent(in) :: fileJ
        character (len=60):: fileout_8
        real(kind=dp), intent(in) :: cutoffmd, cutoffz,cutoffcentre

        real(kind=dp), intent(out), dimension(nq):: result,q_abs
        REAL(kind=dp), DIMENSION(size(q),4*maxval(l)+1,4*maxval(l)+1,4*maxval(l)+1) :: P0matrix
         real(kind=dp),  dimension(maxl*2+1,maxl+1,maxl+1,ng,ng) :: ddx,ddy,ddz
        real(kind=dp), dimension(ng,ng) :: px,py,pz
        real(kind=dp),  dimension(:,:,:,:), allocatable :: zcontr
        real(kind=dp),  dimension(:,:), allocatable :: onerdm_matrix,onerdm_matrix_2
        real(kind=dp),  dimension(:), allocatable :: total,newtotal
        real(kind=dp),  dimension(nq,ngto,ngto) :: e12
        INTEGER(kind=ikind), DIMENSION(maxval(group))   :: group_start, group_count
        integer(kind=ikind), dimension(:), allocatable :: m1, m2, m3, m4
        integer(kind=ikind), dimension(:,:), allocatable :: mat,ep3,ndiff2
        integer(kind=ikind):: nmat,i,j,nmomax,LL,MM,NN,k,c1,sizenmat
        real(kind=dp) ::time1,time2,time3,time4,co,wl,rr,k0,rij
        complex(kind=dp), dimension(:,:,:,:), allocatable :: exponent1, exponent2
        complex(kind=dp), dimension(size(q)):: resultaligned
        real(kind=dp),dimension(3,size(q)) :: q_al
        real(KIND=dp), dimension(size(q)) :: ss,ss2,Iee,Inn,Ine

        print*, 'here we are '

       ! allocate(ep3(size(confs(:,1)),size(confs(:,1))),ndiff2(size(confs(:,1)),size(confs(:,1))) )

        do i = 1, Ng
            do j = 1, Ngto
                if (group(j) == i) then
                    group_start(i) = j
                    group_count(i) = count(group==i)
                    exit
                end if
            end do
        end do
        q_abs=0
       ! call onerdm_creat(confs,civecs,onerdm_matrix,nmomax)
        if (type==1) then
            call cpu_time(time2)
            P0matrix = 0
            CALL set_P0(P0matrix, 4*maxval(l), q)
          !  call maxcoincidence(confs,ep3,ndiff2)
          !  call createtwordm(confs,civecs,ndiff2,ep3,mat,total)
            nmomax=15

                if (read2rdm) then
                sizenmat=numberlines
                open(file=fileJ,unit=15)
                allocate(mat(sizenmat,4), total(sizenmat))
                do i=1,sizenmat
                    read(15,*)mat(i,1), mat(i,2), mat(i,3), mat(i,4), total(i)
                end do

              !  call mcci_to_bit(fileJ,fileout,numberlines)
              !  print*,'Created new file'
                !call one_rdm_two_rdm(mat,total,onerdm_matrix_2)
                close(15)
                else
                    if (fci) then
                call createtwordm_bit_fci(fileJ,numberlines,newdat,irep1,start,end,mat,total)
                elseif (mcci) then
                 fileout_8='newdat.dat'
                call mcci_to_bit(fileJ,fileout_8,numberlines)
                 print*, 'created new file'
                call createtwordm_bit(fileout_8,numberlines,&
                        mat,total)
                else
                    fileout_8='es.dat'
                    print*,fileJ,numberlines

                    call mcci_to_bit(fileJ,fileout_8,numberlines)
                    print*,fileJ,numberlines
                  call createtwordm_bit(fileout_8,numberlines,&
                        mat,total)
            end if
                end if
            allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
             m1 = mat(:,1)
             m2 = mat(:,2)
             m3 = mat(:,3)
             m4 = mat(:,4)
             print*,m1(1), total(1)
             nmat=size(m1)
             print*,'size of nmat',nmat
             nmomax=maxval(mat)
             print*,nmomax
!            call one_rdm_two_rdm(mat,total,onerdm_matrix_2)
!            co=0
!            do i=1,nmomax
!                   co=co+onerdm_matrix_2(i,i)
!               end do
!
!            call one_rdm_bit(fileJ,onerdm_matrix,nmomax,numberlines,newdat,irep1)
!            print*,'diagonal part of the onerdm',co
!            open(file='onerdm_from_twordm.dat', unit=15)
!
!             do i=1,size(onerdm_matrix_2(:,1))
!
!                 write(15,'(1000F14.7)')( onerdm_matrix_2(i,c1),c1=1,size(onerdm_matrix_2(i,:)))
!            end do
!            close(15)


           call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq)

            call cpu_time(time3)
            print*,'Time variables', time3-time2

            call tot_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z1,z2,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result)

            print*,'here'



        else if (type==2) then
            call cpu_time(time2)

               P0matrix = 0

               CALL set_P0(P0matrix, 4*maxval(l), q)
           ! call maxcoincidence(confs,ep3,ndiff2)
            !call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)
               nmomax=27
             !  call one_rdm_slow(fileJ,onerdm_matrix,nmomax)

               call cpu_time(time3)
                print*,'time slow one_rdm',time3-time2

                nmomax=27

            print*, 'number of lines', numberlines
            if (read2rdm) then
                print*, fileJ
                sizenmat=numberlines
                open(file=fileJ,unit=15)
                allocate(mat(sizenmat,4), total(sizenmat))
                do i=1,sizenmat
                    read(15,*)mat(i,1), mat(i,2), mat(i,3), mat(i,4), total(i)
                end do
                   close(15)
                nmomax=maxval(mat)
                call one_rdm_two_rdm(mat,total,onerdm_matrix_2,24)

            else
               call one_rdm_bit(fileJ,onerdm_matrix_2,nmomax,numberlines,newdat,irep1)

            end if

                 call cpu_time(time2)
                print*,'time bit one_rdm', time2-time3
               print*,'maxnmo',nmomax
               co=0
               do i=1,nmomax
                   co=co+onerdm_matrix_2(i,i)
               end do
               print*,'trace',co
            call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,onerdm_matrix_2,nmomax,q,nq)
            call elastic_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result)


!        else if (type==3) then
!            call cpu_time(time2)
!
!
!            call maxcoincidence(confs,ep3,ndiff2)
!            call createtwordm(confs,civecs,ndiff2,ep3,mat,total)
!            allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
!             m1 = mat(:,1)
!             m2 = mat(:,2)
!             m3 = mat(:,3)
!             m4 = mat(:,4)
!             nmat=size(m1)
!
!              allocate(exponent1(nq,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1), &
!                    exponent2(nq,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1))
!
!             wl=4.0_dp*dacos(-1.0_dp)/maxval(q)
!             wl=3.00
!               call linspace(0.0_dp,dacos(-1.0_dp)/2,nq,ss);
!            ss=(ss-dacos(-1.0_dp))/2.0_dp
!            print*,'linspace calculated', wl
!            do mm=1,nq
!                k0 = wl*1.889726125*0.529
!                rr=2*k0*cos(ss(mm))
!                q_al(1,mm)=rr*sin(ss(mm))*cos(0.0_dp)
!                q_al(2,mm)=rr*sin(ss(mm))*sin(0.0_dp)
!                q_al(3,mm)=rr*cos(ss(mm))
!
!                q_abs(mm)=sqrt(q_al(1,mm)**2.0_dp+q_al(2,mm)**2.0_dp+q_al(3,mm)**2.0_dp)
!                end do
!            call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
!        mmod,m1,m2,m3,m4,nmat, total,q_abs,nq)
!
!
!
!            print*,'q_al calculated', q_al(1,:)
!
!
!            do LL=0,maxval(l)*4
!                do MM=0,maxval(m)*4
!                    do NN=0,maxval(n)*4
!
!                    exponent2(:,LL+1,MM+1,NN+1)=((0,-1.00)*q_al(1,:))**(LL)* &
!                            ((0,-1.00)*q_al(2,:))**(MM)*((0,-1.00)*q_al(3,:))**(NN)
!
!                    exponent1(:,LL+1,MM+1,NN+1)=((0,1.00)*q_al(1,:))**(LL) &
!                            *((0,1.00)*q_al(2,:))**(MM)*((0,1.00)*q_al(3,:))**(NN)
!                     enddo
!                enddo
!            enddo
!            print*,'exponents calculated', (exponent1(1,1,1,1))
!
!            call tot_integration_aligned(ng,px,py,pz,l,m,n,ddx,ddy,ddz,z1,z2,group_start,group_count,group, &
!                cutoffz,cutoffmd, cutoffcentre,q_al,e12,exponent1,exponent2,resultaligned)
!
!            result=real(resultaligned)
!
!
!
!         else if (type==4) then
!
!
!             call maxcoincidence(confs,ep3,ndiff2)
!            call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)
!               print*,'maxnmo',nmomax
!             wl=4.0_dp*dacos(-1.0_dp)/maxval(q)
!             wl=3.00
!             call linspace(0.0_dp,dacos(-1.0_dp)/2,nq,ss);
!            ss=(ss-dacos(-1.0_dp))/2.0_dp
!            print*,'linspace calculated', wl
!            do mm=1,nq
!                k0 = wl*1.889726125*0.529
!                rr=2*k0*cos(ss(mm))
!                q_al(1,mm)=rr*sin(ss(mm))*cos(0.0_dp)
!                q_al(2,mm)=rr*sin(ss(mm))*sin(0.0_dp)
!                q_al(3,mm)=rr*cos(ss(mm))
!
!                q_abs(mm)=sqrt(q_al(1,mm)**2.0_dp+q_al(2,mm)**2.0_dp+q_al(3,mm)**2.0_dp)
!            end do
!             co=0
!             print*,'qvector',size(q_abs),maxval(q_abs),maxval(ss)
!               do i=1,nmomax
!                   co=co+onerdm_matrix(i,i)
!               end do
!               print*,'trace',co
!             call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
!        mmod,onerdm_matrix,nmomax,q_abs,nq)
!            print*, 'variables calculated'
!            allocate(exponent1(nq,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1))
!             do LL=0,maxval(l)*4
!                do MM=0,maxval(m)*4
!                    do NN=0,maxval(n)*4
!
!
!                    exponent1(:,LL+1,MM+1,NN+1)=((0,1.00)*q_al(1,:))**(LL) &
!                            *((0,1.00)*q_al(2,:))**(MM)*((0,1.00)*q_al(3,:))**(NN)
!                     enddo
!                enddo
!            enddo
!            print*, 'calling elastic'
!             call elastic_integration_alig(ng,px,py,pz,l,m,n,ddx,ddy,ddz,z,group_start,group_count,group, &
!                cutoffz,cutoffmd, cutoffcentre,q_al,e12,exponent1,resultaligned)
!                 result=abs(resultaligned)**2
!
!
!
!            !Starting the rotationally averaged elastic electron scattering, in the end 4 different versions must be programmed
!
!            else if (type==5) then
!
!                P0matrix = 0
!                CALL set_P0(P0matrix, 4*maxval(l), q)
!                call maxcoincidence(confs,ep3,ndiff2)
!                call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)
!
!                call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
!        mmod,onerdm_matrix,nmomax,q,nq)
!
!                print*,'first Z', Z(1,1)
!
!                call elastic_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
!                cutoffz,cutoffmd, cutoffcentre,q,e12,Iee)
!
!
!
!                print*,Iee(1),geom
!
!                call nuclei_electron_integration(Zn,geom,ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
!                cutoffz,cutoffmd, cutoffcentre,q,e12,Ine)
!
!                print*,Ine(1)
!
!                Inn=0.0_dp
!                do i=1,size(Zn)
!                    do j=1,size(Zn)
!                        do k=1,nq
!                         Rij=sqrt(sum((geom(i,:)-geom(j,:))**2) )
!                         Inn(k)=Inn(k)+Zn(i)*Zn(j)* sinc(q(k)*abs(Rij))
!
!                         end do
!                    enddo
!                end do
!
!                print*,Inn(1)
!
!                print*,sum(Ine-Iee)
!                result=Iee-2*Ine+Inn
!
!
!
!
!
!
!                !we need three terms here, Inn, Iee and Ine
!
!
!
        end if
        end subroutine total_scattering_calculation_2

        subroutine linspace(from, to, n,array)
                INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
         INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
         integer,intent(in) :: n
         real(kind=dp), intent(in) :: from, to
         real(kind=dp), intent(out),dimension(n) :: array
        real(kind=dp) :: range

    integer :: i

    range = to - from

    do i=1, n
        array(i) = from + range * (i - 1) / (n - 1)
    end do
    end subroutine

     function sinc (a)
             INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
            INTEGER, PARAMETER :: ikind = SELECTED_INT_KIND(8)
             real(kind=dp):: sinc, a
             if (abs(a) < 1.0d-10) then
                sinc = 1
             else
                sinc = sin(a) / (a)
             end if
    end function

        end module main_calculation_mod
