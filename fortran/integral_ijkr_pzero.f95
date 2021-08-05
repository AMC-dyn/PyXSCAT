!-----------------------------------------------------------------------
! integral_ijkr_pzero
! Andres Moreno Carrascosa and Mats Simmermacher, 2021
!-----------------------------------------------------------------------

MODULE integral_ijkr_pzero

IMPLICIT NONE

CONTAINS
     
    SUBROUTINE integral_ijkr_pzero(nq,lmax1,lmax2,lmax3,lmax4,p0mat,dx,dy,dz,i,j,k,r,z1,z2,apos,cutoffz,cutoffmd)
    
        ! definition of input
        INTEGER(kind=ikind), INTENT(IN)                       :: nq, lmax1, lmax2, lmax3, lmax4, i, j, k, r
        REAL(kind=dp), INTENT(IN)                             :: cutoffz, cutoffmd
        REAL(kind=dp), INTENT(IN), DIMENSION(:)               :: p0mat
        REAL(kind=dp), INTENT(IN), DIMENSION(:)               :: dx, dy, dz
        REAL(kind=dp), INTENT(IN), DIMENSION(:)               :: z1, z2
        INTEGER(kind=ikind), INTENT(IN), DIMENSION(:)         :: apos
        ! definition of output
        REAL(kind=dp), INTENT(OUT), DIMENSION(nq)             :: itgr 
        ! definition of loop indices    
        INTEGER(kind=ikind)                                   :: l1, m1, l2, m2, l3, m3, l4, m4
        INTEGER(kind=ikind)                                   :: l, m, n, lp, mp, np
        ! definition of internal variables
        INTEGER(kind=ikind)                                   :: n1, n2, n3, n4, h1i, posj, posk, posr
        REAL(kind=dp)                                         :: ztot
        REAL(kind=dp)                                         :: mdl, mdm, mdn, mdlp, mdmp, prodd
        REAL(kind=dp), DIMENSION(nq)                          :: zij1, zkr1
        REAL(kind=dp), DIMENSION(nq)                          :: zij2, zkr2
        REAL(kind=dp), DIMENSION(nq)                          :: f

        ! loop through all possible ways to get total angular momentum lmax1
        do l1 = 0, lmax1
            do m1 = 0, lmax1-l1
                n1 = lmax1-l1-m1
                posj = apos(j)
                ! loop through all possible ways to get total angular momentum lmax2
                do l2 = 0, lmax2
                    do m2 = 0, lmax2-l2
                        n2 = lmax2-l2-m2
                        zij1 = z1(:,posi,posj)
                        zij2 = z2(:,posi,posj)
                        posk = apos(k)
                        ! loop through all possible ways to get total angular momentum lmax3
                        do l3 = 0, lmax3
                            do m3 = 0, lmax3-l3
                                n3 = lmax3-l3-m3
                                posr = apos(r)
                                ! loop through all possible ways to get total angular momentum lmax4
                                do l4 = 0, lmax4
                                    do m4 = 0, lmax4-l4
                                        n4 = lmax4-l4-m4
                                        zkr1 = z1(:,posk,posr)
                                        zkr2 = z2(:,posk,posr)
                                        ! total prefactor    
                                        ztot = sum(zij1*zkr2 + zij2*zkr1) / 8
                                        ! continue only if larger
                                        if abs(ztot) < cutoffz then
                                            posr = posr+1   
                                            cycle 
                                        end if
                                        ! the 6-dimensional sum over MD coefficents
                                        do l = 0, l1+l2
                                            mdl = dx(i,j,l+1,l1+1,l2+1) * ztot
                                            if mdl == 0 cycle
                                            do m = 0, m1+m2
                                                mdm = dy(i,j,m+1,m1+1,m2+1) * mdl
                                                if mdm == 0 cycle
                                                do n = 0, n1+n2
                                                    h1 = (-1)^(l+m+n)
                                                    mdn = dz(i,j,n+1,n1+1,n2+1) * mdm * h1
                                                    if mdn == 0 cycle
                                                    do lp = 0, l3+l4
                                                        mdlp = dx(k,r,lp+1,l3+1,l4+1) * mdn
                                                        if mdlp == 0 cycle
                                                        do mp = 0, m3+m4
                                                            mdmp = dy(k,r,mp+1,m3+1,m4+1) * mdlp
                                                            if mdmp == 0 cycle
                                                            do np = 0, n3+n4
                                                                prodd = dz(k,r,np+1,n3+1,n4+1) * mdmp
                                                                ! cutoff after md
                                                                if abs(prodd) < cutoffmd cycle
                                                                ! add the contribution to the total
                                                                f = p0mat(:,l+lp+1,m+mp+1,n+np+1)
                                                                itgr = itgr + prodd*f
                                                            end do
                                                        end do
                                                    end do
                                                end do
                                            end do
                                        end do
                                        posr = posr+1
                                    end do
                                end do
                                posk = posk+1
                            end do
                        end do
                        posj = posj+1
                    end do
                end do
                posi = posi+1
            end do
        end do

    END SUBROUTINE integral_ijkr_pzero

END MODULE integral_ijkr_pzero

