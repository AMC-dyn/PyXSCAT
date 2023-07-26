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



      use mod_types
    use omp_lib
    use onerdm

    use twordms

    use variables

    use integrals
!    use p0_cases

 !   use bessel_calcs

 !   use total_alig

!    use total_j0

 !   use total_j1

 !   use total_j2

 !   use elastic_alig

 !   use elastic_j0

 !   use elastic_j2

    implicit none

        
        

        integer(kind=ikind), intent(in):: ngto, ng,  nq, maxl,type, state1,state2
        integer(kind=ikind), intent(in),dimension(:),allocatable :: l, m, n,group
        integer(kind=ikind), dimension(:,:), allocatable,intent(in):: confs

        real(kind=dp), intent(in),dimension(:),allocatable :: ga, xx, yy, zz, q
        real(kind=dp), intent(in), dimension(:), allocatable:: zn
        real(kind=dp), intent(in),dimension(:,:), allocatable :: mmod, civecs,geom
        real(kind=dp),dimension(:,:,:),allocatable :: z1,z2,z1t,z2t
        real(kind=dp),dimension(ngto,ngto):: z

        real(kind=dp), intent(in) :: cutoffmd, cutoffz,cutoffcentre
        real(kind=dp), dimension(nq):: result2
        real(kind=dp), intent(out), dimension(:), allocatable:: result,q_abs
        REAL(kind=dp), DIMENSION(size(q),4*maxval(l)+1,4*maxval(l)+1,4*maxval(l)+1) :: P0matrix
         real(kind=dp),  dimension(maxl*2+1,maxl+1,maxl+1,ng,ng) :: ddx,ddy,ddz
        real(kind=dp), dimension(ng,ng) :: px,py,pz
        real(kind=dp),  dimension(:,:,:,:), allocatable :: zcontr
        real(kind=dp),  dimension(:,:), allocatable :: onerdm_matrix
        real(kind=dp),  dimension(:), allocatable :: total,newtotal
        real(kind=dp),  dimension(nq,ng,ng) :: e12
        real(kind=dp),  dimension(3131,ngto,ngto):: E123
        INTEGER(kind=ikind), DIMENSION(maxval(group))   :: group_start, group_count
        integer(kind=ikind), dimension(:), allocatable :: m1, m2, m3, m4
        integer(kind=ikind), dimension(:,:), allocatable :: mat,ep3,ndiff2
        integer(kind=ikind):: nmat,i,j,nmomax,LL,MM,NN,k
        real(kind=dp) :: start,time1,time2,time3,time4,co,wl,rr,k0,rij
        complex(kind=dp), dimension(:,:,:,:), allocatable :: exponent1, exponent2
        complex(kind=dp), dimension(3131):: resultaligned
        real(kind=dp),dimension(3,3131) :: q_al
        real(kind=dp),dimension(3131):: q_abs_2
        real(KIND=dp), dimension(nq) :: ss,ss2,Iee,Inn,Ine

         print*,nq
         print*,ngto,size(E123(:,1,1)),size(E123(1,:,1)), size(E123(1,1,:))
         print*,ngto,size(E12(:,1,1)),size(E12(1,:,1)), size(E12(1,1,:))
        do i = 1, Ng
            do j = 1, Ngto
                if (group(j) == i) then
                    group_start(i) = j
                    group_count(i) = count(group==i)
                    exit
                end if
            end do
        end do


        !Now we explore the different cases
        allocate(q_abs(nq), result(nq))
              q_abs=0.0_dp
        select case(type)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TOTAL J0 SCATTERING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CASE(1)

            P0matrix = 0.0_dp
            CALL set_P0(P0matrix, 4*maxval(l), q)
           ! CALL set_P0(P0matrix, 4*maxval(l), q)
            call maxcoincidence(confs,ep3,ndiff2)

            call createtwordm(confs,civecs,ndiff2,ep3,mat,total,state1,state2)

            allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
             m1 = mat(:,1)
             m2 = mat(:,2)
             m3 = mat(:,3)
             m4 = mat(:,4)
             nmat=size(m1)

            print*,'size nmat', nmat
            call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,z1t,z2t,&
                    e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq)

            call tot_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z1,z2,z1t,z2t,&
                    group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result2)

            result=result2

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ELASTIC J0 SCATTERING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       CASE(2)
               P0matrix = 0.0_dp

            CALL set_P0(P0matrix, 4*maxval(l), q)
            call maxcoincidence(confs,ep3,ndiff2)
            call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)
            print*,'onerdm ayyy', maxval(onerdm_matrix)
            call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,onerdm_matrix,nmomax,q,nq)
            print*,maxval(e12), 'and this is the maximum value'

            call elastic_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result2)

        result=result2
               print*,'CALCULATION FINISHED', q(1),result(1)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TOTAL ALIGNED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        CASE(3)



            call maxcoincidence(confs,ep3,ndiff2)
            call createtwordm(confs,civecs,ndiff2,ep3,mat,total,state1,state2)
            allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
             m1 = mat(:,1)
             m2 = mat(:,2)
             m3 = mat(:,3)
             m4 = mat(:,4)
             nmat=size(m1)
            allocate(q_abs(3131))
             call QCOORD3(10.583544175_dp ,q_al(1,:),q_al(2,:),q_al(3,:))
              allocate(exponent1(3131,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1), &
                    exponent2(3131,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1))

                q_abs=sqrt(q_al(1,:)**2.0_dp+q_al(2,:)**2.0_dp+q_al(3,:)**2.0_dp)
            do LL=1,3131
                write(*,'(4F15.8)') q_al(1,LL), q_al(2,LL), q_al(3,LL), dsqrt(q_al(1,LL)**2.0_dp+ q_al(2,LL)**2.0_dp+q_al(3,LL)**2.0_dp)
            end do
            e123=0.0_dp
            print*,'size q_abs',size(q_abs_2), 3131
            PRINT*,SIZE(e123)
             print*,ngto,size(E123(:,1,1)),size(E123(1,:,1), size(E123(1,1,:)))
            call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,z1t,z2t,&
                    e123,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q_abs,3131)






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


            call tot_integration_aligned(ng,px,py,pz,l,m,n,ddx,ddy,ddz,z1,z2,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q_al,e123,exponent1,exponent2,resultaligned)

            result=0
            open(unit=16, file="result_aligned.dat")
            do LL=1,3131
                write(16,'(5G20.12)') q_al(1,LL), q_al(2,LL), q_al(3,LL), real(resultaligned(LL)),imag(resultaligned(LL))
            end do
            close(16)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ELASTIC ALIGNED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            q_abs=0.0_dp
         CASE(4)


             call maxcoincidence(confs,ep3,ndiff2)
            call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)

             wl=4.0_dp*dacos(-1.0_dp)/maxval(q)
             wl=3.00
             call linspace(0.0_dp,dacos(-1.0_dp)/2,nq,ss);
            ss=(ss-dacos(-1.0_dp))/2.0_dp

            do mm=1,nq
                k0 = wl*1.889726125*0.529
                rr=2*k0*cos(ss(mm))
                q_al(1,mm)=rr*sin(ss(mm))*cos(0.0_dp)
                q_al(2,mm)=rr*sin(ss(mm))*sin(0.0_dp)
                q_al(3,mm)=rr*cos(ss(mm))

                q_abs(mm)=sqrt(q_al(1,mm)**2.0_dp+q_al(2,mm)**2.0_dp+q_al(3,mm)**2.0_dp)
            end do
             co=0

               do i=1,nmomax
                   co=co+onerdm_matrix(i,i)
               end do

             call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,onerdm_matrix,nmomax,q_abs,nq)

            allocate(exponent1(nq,maxval(l)*4+1,maxval(l)*4+1,maxval(l)*4+1))
             do LL=0,maxval(l)*4
                do MM=0,maxval(m)*4
                    do NN=0,maxval(n)*4


                    exponent1(:,LL+1,MM+1,NN+1)=((0,1.00)*q_al(1,:))**(LL) &
                            *((0,1.00)*q_al(2,:))**(MM)*((0,1.00)*q_al(3,:))**(NN)
                     enddo
                enddo
            enddo

             call elastic_integration_alig(ng,px,py,pz,l,m,n,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q_al,e12,exponent1,resultaligned)
                 result=abs(resultaligned)**2




            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TOTAL ELECTRON SCATTERING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            CASE(5)
                P0matrix = 0
                CALL set_P0(P0matrix, 4*maxval(l), q)
                call maxcoincidence(confs,ep3,ndiff2)
                call createtwordm(confs,civecs,ndiff2,ep3,mat,total,state1,state2)
                allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
                m1 = mat(:,1)
                m2 = mat(:,2)
                m3 = mat(:,3)
                m4 = mat(:,4)
                nmat=size(m1)
                call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,z1t,z2t,&
                        e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
                    mmod,m1,m2,m3,m4,nmat, total,q,nq)

                call tot_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z1,z2,z1t,z2t,&
                        group_start,group_count,group, &
                    cutoffz,cutoffmd, cutoffcentre,q,e12,Iee)

                call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)

                call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
                    mmod,onerdm_matrix,nmomax,q,nq)


                call nuclei_electron_integration(Zn,geom,ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,Ine)


                Inn=0.0_dp
                do i=1,size(Zn)
                    do j=1,size(Zn)
                        do k=1,nq
                         Rij=sqrt(sum((geom(i,:)-geom(j,:))**2) )
                         Inn(k)=Inn(k)+Zn(i)*Zn(j)* sinc(q(k)*abs(Rij))

                         end do
                    enddo
                end do

                  result=Iee-2*Ine+Inn
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ELASTIC ELECTRON SCATTERING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            CASE(6)

                P0matrix = 0
                CALL set_P0(P0matrix, 4*maxval(l), q)
                call maxcoincidence(confs,ep3,ndiff2)
                call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)

                call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,onerdm_matrix,nmomax,q,nq)



                call elastic_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,Iee)





                call nuclei_electron_integration(Zn,geom,ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,Ine)


                Inn=0.0_dp
                do i=1,size(Zn)
                    do j=1,size(Zn)
                        do k=1,nq
                         Rij=sqrt(sum((geom(i,:)-geom(j,:))**2) )
                         Inn(k)=Inn(k)+Zn(i)*Zn(j)* sinc(q(k)*abs(Rij))

                         end do
                    enddo
                end do




                result=Iee-2*Ine+Inn







          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TOTAL J2 CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            CASE(7)

               P0matrix = 0.0_dp
                CALL set_P0_j2(P0matrix, 4*maxval(l), q)
           ! CALL set_P0(P0matrix, 4*maxval(l), q)
               call maxcoincidence(confs,ep3,ndiff2)
               call createtwordm(confs,civecs,ndiff2,ep3,mat,total,state1,state2)

               allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
               m1 = mat(:,1)
               m2 = mat(:,2)
               m3 = mat(:,3)
               m4 = mat(:,4)
               nmat=size(m1)
               call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,z1t,z2t, &
               e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq)
               call tot_integration_j2(ng,nq,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z1,z2,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result2)


               result=result2
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ELASTIC J2 CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CASE(8)
                P0matrix = 0

                CALL set_P0_j2(P0matrix, 4*maxval(l), q)
                call maxcoincidence(confs,ep3,ndiff2)
                call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)

                call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
                     mmod,onerdm_matrix,nmomax,q,nq)
                call elastic_integration_j2(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                      cutoffz,cutoffmd, cutoffcentre,q,e12,result2)
            result=result2
        CASE(9)

               P0matrix = 0.0_dp
               print*,'maximum value l', 4*maxval(l)
                CALL set_P0_j1(P0matrix, 4*maxval(l), q)
           ! CALL set_P0(P0matrix, 4*maxval(l), q)
               call maxcoincidence(confs,ep3,ndiff2)
               call createtwordm(confs,civecs,ndiff2,ep3,mat,total,state1,state2)

               allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
               m1 = mat(:,1)
               m2 = mat(:,2)
               m3 = mat(:,3)
               m4 = mat(:,4)
               nmat=size(m1)
               call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,z1t,z2t,&
                       e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, total,q,nq)

               print*,'Going to integration'
               call tot_integration_j1(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z1,z2,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result2)
            result=result2

        end select
        end subroutine total_scattering_calculation


        subroutine total_scattering_calculation_2(type,Zn,geom,state1,state2,maxl,ngto,ng,ga,l,m,n,xx,yy,zz, &
        mmod,q,nq, group,&
        cutoffz,cutoffmd,cutoffcentre,fileJ,numberlines,newdat,start,end,start_2, end_2, fci,irep1,irep2,ordering1 &
                ,ordering2, read2rdm,mcci,q_abs,result)
              use mod_types




    use onerdm

    use twordms

    use variables

    use integrals


    implicit none

        
        
        logical, intent(in) :: fci,read2rdm,mcci
        logical :: read2rdm2
        integer(kind=ikind), intent(in):: ngto, ng,  nq, maxl,type, state1,state2,numberlines
        integer(kind=ikind), intent(in),dimension(:),allocatable :: l, m, n,group
        integer(kind=ikind), intent(in), dimension(:,:):: newdat
        integer(kind=ikind),intent(in), dimension(:):: irep1,start,end, irep2, start_2, end_2,ordering1,ordering2


        real(kind=dp), intent(in),dimension(:), allocatable :: ga, xx, yy, zz, q, Zn
        real(kind=dp), intent(in),dimension(:,:),allocatable :: mmod, geom
        real(kind=dp),dimension(:,:,:),allocatable :: z1,z2,z1t,z2t
        real(kind=dp),dimension(ngto,ngto):: z
        character (len=100), intent(in) :: fileJ
        character (len=60):: fileout_8, fileJ2
        real(kind=dp), intent(in) :: cutoffmd, cutoffz,cutoffcentre

        real(kind=dp), intent(out), dimension(:), allocatable:: result,q_abs
        REAL(kind=dp), DIMENSION(size(q),4*maxval(l)+1,4*maxval(l)+1,4*maxval(l)+1) :: P0matrix
         real(kind=dp),  dimension(maxl*2+1,maxl+1,maxl+1,ng,ng) :: ddx,ddy,ddz
        real(kind=dp), dimension(ng,ng) :: px,py,pz
        real(kind=dp),  dimension(:,:,:,:), allocatable :: zcontr
        real(kind=dp),  dimension(:,:), allocatable :: onerdm_matrix,onerdm_matrix_2
        real(kind=dp),  dimension(:), allocatable :: total,newtotal, totaldiv
        real(kind=dp),  dimension(nq,ng,ng) :: e12
        INTEGER(kind=ikind), DIMENSION(maxval(group))   :: group_start, group_count
        integer(kind=ikind), dimension(:), allocatable :: m1, m2, m3, m4
        integer(kind=ikind), dimension(:,:), allocatable :: mat,ep3,ndiff2
        integer(kind=ikind):: nmat,i,j,nmomax,LL,MM,NN,k,c1,sizenmat,spincont,counter,reads,fin,ini,npoints,ndivs
        real(kind=dp) ::time1,time2,time3,time4,co,wl,rr,k0,rij
        complex(kind=dp), dimension(:,:,:,:), allocatable :: exponent1, exponent2
        complex(kind=dp), dimension(nq):: resultaligned
        real(kind=dp),dimension(3,nq) :: q_al
        real(KIND=dp), dimension(nq) :: ss,ss2,Iee,Inn,Ine
        real(kind=dp), dimension(nq):: result2

        allocate(q_abs(nq), result(nq))
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
                reads=numberlines
                print*,sizenmat,fileJ
!                open(file=fileJ,unit=15)
!                   do i=1,reads
!                    read(15,*)spincont
!                     if (spincont==1) then
!                        sizenmat=sizenmat+1
!                    end if
!
!                   end do
!                close(15)
                open(file=fileJ,unit=15)
                allocate(mat(sizenmat,4), total(sizenmat))
                counter=1
                do i=1,reads
                    read(15,*)mat(counter,1), mat(counter,2), mat(counter,3), mat(counter,4), total(counter)
!                    if (spincont==1) then
!
!                        mat(counter,1)=mat(counter,1)+18
!                        mat(counter,2)=mat(counter,2)+18
!                        total(counter+1)=total(counter)
!                        mat(counter+1,1)=mat(counter,1)
!                        mat(counter+1,2)=mat(counter,2)
!                        mat(counter+1,3)=mat(counter,3)+18
!                        mat(counter+1,4)=mat(counter,4)+18
!                        counter=counter+2
!                    elseif (spincont==2) then
!                        mat(counter,3)=mat(counter,3)+18
!                        mat(counter,4)=mat(counter,4)+18
!                        mat(counter,2)=mat(counter,2)+18
!                        mat(counter,1)=mat(counter,1)+18
!                        counter=counter+1
!                    else
                        counter=counter+1
                    !end if
                end do
                print*,counter,sizenmat
              !  call mcci_to_bit(fileJ,fileout,numberlines)
              !  print*,'Created new file'
                !call one_rdm_two_rdm(mat,total,onerdm_matrix_2)
                close(15)
                else
                    if (fci) then
                    call createtwordm_bit_fci(fileJ,numberlines,newdat,irep1,start,end,mat,total)
                elseif (mcci) then
                       call cpu_time(time2)
                 fileout_8='newdat.dat'
                call mcci_to_bit(fileJ,fileout_8,numberlines)

                write(*,*) 'entering createtwordm_bit'

                 call createtwordm_bit(fileout_8,numberlines,&
                        mat,total)
                 call cpu_time(time3)
                 print*,'Time 2rdm', time3-time2
                else
                  write(*,*) 'going to create es.dat'
                  call cpu_time(time2)
                    fileout_8='es.dat'


                    call mcci_to_bit(fileJ,fileout_8,numberlines)
                    write(*,*) 'done mcci_to_bit'

                  call createtwordm_bit(fileout_8,numberlines,&
                        mat,total)
                    write(*,*) 'created twordm'
               call cpu_time(time3)
            print*,'Time 2rdm', time3-time2,numberlines
            end if
                end if


            ndivs=1
            npoints=numberlines/ndivs
            result=0.d0
            do i=1,ndivs
                print*,i, 'step of', ndivs

                ini=1+(numberlines/ndivs*(i-1))
                fin=1+numberlines/ndivs*i
                if (i==ndivs) then
                    fin=numberlines
                end if
            write(*,*) 'allocating...'
                npoints=fin-ini
                write(*,*) ini,fin,npoints, size(mat)
                allocate(m1(npoints), m2(npoints), m3(npoints), m4(npoints), totaldiv(npoints))
              write(*,*) 'allocated m{1..4}'
                m1 = mat(:,1)
              write(*,*) 'assigned m1'
                m2 = mat(:,2)
              write(*,*) 'assigned m2'
                m3 = mat(:,3)
              write(*,*) 'assigned m3'
                m4 = mat(:,4)
              write(*,*) 'assigned m4'

                nmat=size(m1)

                nmomax=maxval(m1)
                totaldiv=total(:)
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


            write(*,*) 'calling variables_total'
           call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,z1t,z2t, &
                   e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,m1,m2,m3,m4,nmat, totaldiv,q,nq)

            call cpu_time(time3)
            print*,'Time variables', time3-time2
            print*,'size of wavefunction', size(l)

            call tot_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z1,z2,z1t,z2t,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result2)
            result=result+result2
            print*,result(1)
            deallocate(m1,m2,m3,m4, totaldiv)

        enddo

        else if (type==2) then
            call cpu_time(time2)

               P0matrix = 0

               CALL set_P0(P0matrix, 4*maxval(l), q)
           ! call maxcoincidence(confs,ep3,ndiff2)
            !call onerdm_creat(confs,civecs,onerdm_matrix,nmomax,state1,state2)
               nmomax=25
             !  call one_rdm_slow(fileJ,onerdm_matrix,nmomax)

               call cpu_time(time3)
                print*,'time slow one_rdm',time3-time2

                nmomax=25

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
                print*,'calling onerdm_bit'
               call one_rdm_bit(fileJ,onerdm_matrix_2,nmomax,numberlines,newdat,irep1)

            end if

                 call cpu_time(time2)
                print*,'time bit one_rdm', time2-time3

               co=0
               do i=1,nmomax
                   co=co+onerdm_matrix_2(i,i)
               end do

            call variables_elastic(px,py,pz,ddx,ddy,ddz,z,e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
        mmod,onerdm_matrix_2,nmomax,q,nq)
            call elastic_integration(ng,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z,group_start,group_count,group, &
                cutoffz,cutoffmd, cutoffcentre,q,e12,result)



           else if (type==7) then

                P0matrix = 0.0_dp
                CALL set_P0_j2(P0matrix, 4*maxval(l), q)
                ! CALL set_P0(P0matrix, 4*maxval(l), q)
                nmomax=19
                read2rdm2=.False.
                !fileJ2='twordm_fortran_bit_2.dat'
                if (read2rdm2) then
                    sizenmat=10950
                    open(file=fileJ2,unit=15)
                    allocate(mat(sizenmat,4), total(sizenmat))
                    do i=1,sizenmat
                        read(15,*)mat(i,1), mat(i,2), mat(i,3), mat(i,4), total(i)
                    end do
                    print*,'supergood'
                    !  call mcci_to_bit(fileJ,fileout,numberlines)
                    !  print*,'Created new file'
                    !call one_rdm_two_rdm(mat,total,onerdm_matrix_2)
                    close(15)
                else
                    if (fci) then
                        call createtwordm_bit_fci(fileJ,numberlines,newdat,irep1,start,end,mat,total)
                    elseif (mcci) then
                        call cpu_time(time2)
                        fileout_8='newdat.dat'
                        call mcci_to_bit(fileJ,fileout_8,numberlines)

                        call createtwordm_bit(fileout_8,numberlines,&
                                mat,total)
                        call cpu_time(time3)
                        print*,'Time 2rdm', time3-time2
                    else
                        call cpu_time(time2)
                        fileout_8='es.dat'


                         call mcci_to_bit(fileJ,fileout_8,numberlines)

                        call createtwordm_bit(fileJ,numberlines,&
                                mat,total)
                        call cpu_time(time3)
                        print*,'Time 2rdm', time3-time2,numberlines
                    end if
                end if
                 allocate(m1(size(mat(:,1))), m2(size(mat(:,1))), m3(size(mat(:,1))), m4(size(mat(:,1))))
                m1 = mat(:,1)
                m2 = mat(:,2)
                m3 = mat(:,3)
                m4 = mat(:,4)

                nmat=size(m1)

                nmomax=maxval(mat)
                print*,'up to variables'
                call variables_total(px,py,pz,ddx,ddy,ddz,z1,z2,z1t,z2t, &
                e12,maxl, ngto,ng,group_start,group_count,group,ga,l,m,n,xx,yy,zz, &
                        mmod,m1,m2,m3,m4,nmat, total,q,nq)

                call tot_integration_j2(ng,nq,px,py,pz,l,m,n,p0matrix,ddx,ddy,ddz,z1,z2,group_start,group_count,group, &
                        cutoffz,cutoffmd, cutoffcentre,q,e12,result)
            end if
        end subroutine total_scattering_calculation_2



     function sinc (a)
              use mod_types
             
            
             real(kind=dp):: sinc, a
             if (abs(a) < 1.0d-10) then
                sinc = 1
             else
                sinc = sin(a) / (a)
             end if
    end function

    subroutine QCOORD3(k0,qx,qy,qz)
              use mod_types
        implicit none
         
         
         real(kind=dp), intent(in):: k0
         real(kind=dp),  dimension(3131), intent(out):: qx, qy,qz
         real(kind=dp):: ldk0,rc,qqx,qqy,ang,pi
         integer(kind=ikind), parameter::n=30
         integer(kind=ikind):: np, len,i,j,pos
         integer(kind=ikind),dimension(n):: npt



        pi= dacos(-1.00_dp)
        len=7
        npt=0
        ldk0 = dsqrt(3.0_dp)/2.0_dp * k0/(n+1.0_dp)
        !ldk0=k0/(n+1)
        do i = 2,(n+1)

             np = ceiling(pi/asin(1.d0/(2.d0*i)))

             NPT(i-1) = np
             len = len + np
        enddo


        qx=1.0D-10
        qy=1.0D-10
        qz=1.0D-10
        QY(2) = ldk0
        QZ(2) = k0*(1-dsqrt(1.d0-1.d0/(n+1.d0)**2))
        pos = 2

        do i=1,5
            pos=pos+1
            qqx = dsin(i*pi/3.0_dp)*ldk0
            qqy = dcos(i*pi/3.0_dp)*ldk0
            rc=i*ldk0
            qx(pos) = qqx
            qy(pos) = qqy
            qz(pos) = k0*(1.0_dp-sqrt(1.d0-(qqx**2 + qqy**2)/k0**2))
        end do

        do i = 2,(n+1)

            pos = pos + 1
            rc = i*ldk0
            ang = 2.0_dp*pi/NPT(i-1)
            qy(pos) = rc
            qz(pos) = k0*(1.d0 - dsqrt(abs(1.0_dp - i**2.0_dp/(n+1.0_dp)**2)))

            do j = 1,(NPT(i-1)-1)

                pos = pos + 1
                qqx = dsin(j*ang)*rc
                qqy = dcos(j*ang)*rc
                qx(pos) = qqx
                qy(pos) = qqy
                qz(pos) = k0*(1.0_dp - dsqrt(abs(1.0_dp - (qqx**2 + qqy**2)/k0**2)))
            enddo

        enddo


    end subroutine QCOORD3

        end module main_calculation_mod
