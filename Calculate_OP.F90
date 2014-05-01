SUBROUTINE Calculate_OP
    USE NEGF_variables
    IMPLICIT NONE
    REAL :: dumVar

    !Reset values to be calculated each iteration after eigenvector calculation
    DeltaAS(xx,yy,:) = c0; DeltaBS(xx,yy,:) = c0; DeltaABS(xx,yy,:) = c0
    DeltaAPUp(xx,yy,:) = c0; DeltaAPDn(xx,yy,:) = c0;
    DeltaBPUp(xx,yy,:) = c0; DeltaBPDn(xx,yy,:) = c0;
    DeltaABPUp(xx,yy,:) = c0; DeltaABPDn(xx,yy,:) = c0;
    DeltaMat(xx,yy,:,:,:) = c0;

    eigvecs = Htot       !Input eigvecs will be overwritten with complex eigenvectors of Htot
    CALL get_evalsevecs  !eigvals listed in ascending order, eigenvectors stored as column vectors in eigvecs matrix

    !Order parameter for s-wave/p-wave coupling in orbital A/B,
    DO ii = NN/2+1, NN  !sum over only positive energy eigenvalues, see HHH's PRB 87 035401 Eq 10
        DO zz = 1, Npz
            index = (zz-1)*Norb
            DO jj = 1, Norb/2
                DO kk = 1, Norb/2
                    DeltaMat(xx,yy,zz,jj,kk) = DeltaMat(xx,yy,zz,jj,kk) &
                        - conjg(eigvecs(index+Norb/2+kk, ii)) * eigvecs(index+jj, ii) * tanh(eigvals(ii)/(kB*T0*2D0))
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    DO ii = 1,Norb/2
        DO jj = 1, Norb/2
            DeltaAS(xx,yy,:) = DeltaAS(xx,yy,:) + DeltaMat(xx,yy,:,ii,jj)*gDAS(ii,jj)/2D0
            DeltaBS(xx,yy,:) = DeltaBS(xx,yy,:) + DeltaMat(xx,yy,:,ii,jj)*gDBS(ii,jj)/2D0
            DeltaABS(xx,yy,:)= DeltaABS(xx,yy,:) + DeltaMat(xx,yy,:,ii,jj)*gDABS(ii,jj)/4D0

            DeltaAPUp(xx,yy,:) = DeltaAPUp(xx,yy,:) + DeltaMat(xx,yy,:,ii,jj)*gDAPUp(ii,jj)
            DeltaAPDn(xx,yy,:) = DeltaAPUp(xx,yy,:) + DeltaMat(xx,yy,:,ii,jj)*gDAPDn(ii,jj)

            DeltaBPUp(xx,yy,:) = DeltaBPUp(xx,yy,:) + DeltaMat(xx,yy,:,ii,jj)*gDBPUp(ii,jj)
            DeltaBPDn(xx,yy,:) = DeltaBPDn(xx,yy,:) + DeltaMat(xx,yy,:,ii,jj)*gDBPDn(ii,jj)

            DeltaABPUp(xx,yy,:) = DeltaABPUp(xx,yy,:) + DeltaMat(xx,yy,:,ii,jj)*gDABPUp(ii,jj)/2D0
            DeltaABPDn(xx,yy,:) = DeltaABPDn(xx,yy,:) + DeltaMat(xx,yy,:,ii,jj)*gDABPDn(ii,jj) /2D0
        ENDDO
    ENDDO

    !Calculate change from previous iteration to determine self-consistency
    IF (ite >= 10) THEN
        AmtChange = 0D0; RMSChange = 0D0

        DO zz = 1, Npz
            DO ii = 1, Norb/2
                DO jj = 1, Norb/2
                    dumVar = abs(DeltaMat(xx,yy,zz,ii,jj)-DeltaMatOld(zz,ii,jj)) !/maxval(abs(DeltaMat(xx,yy,zz,:,:)))
                    RMSChange = RMSChange+abs(dumVar)**2
                    IF (dumVar > AmtChange) THEN
                        AmtChange = dumVar
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

        RMSChange = sqrt(RMSChange/(Npz*Norb/2*Norb/2)) !9D0 because 9 pairing schemes ({A,B,AB}*{s,pup,pdn})
    ENDIF

    DO zz = 1, Npz
        DO ii = 1, Norb/2
            DO jj = 1, Norb/2
                 DeltaMatOld(zz,ii,jj) = DeltaMat(xx,yy,zz,ii,jj);
            ENDDO
        ENDDO
    ENDDO
END SUBROUTINE Calculate_OP


SUBROUTINE CalculateFinal
    USE NEGF_variables
    IMPLICIT NONE
    REAL :: dumVar
    COMPLEX, DIMENSION(:), ALLOCATABLE :: PsiKet

    !!!!!!!!!!Spin-resolved and orbital-resolved LDOS and dispersion relations
    ALLOCATE(PsiKet(Norb/2));  PsiKet=c0
    DO zz = 1, Npz
        index = (zz-1)*Norb
        DO ii = 1, NN
            PsiKet = eigvecs(index+1:index+Norb/2, ii)

            DO ee = 1, NpE
                DO kk = 1, Norb/2
                    DO ll = 1, Norb/2
                        LDOS(xx,yy,zz,ee)   = LDOS(xx,yy,zz,ee)   + exp(-1D3*(Eaxis(ee)-eigvals(ii))**2) * &
                            real( conjg(PsiKet(kk))*PsiKet(ll)* id(kk,ll) )

                        LDOSMx(xx,yy,zz,ee) = LDOSMx(xx,yy,zz,ee) + exp(-1D3*(Eaxis(ee)-eigvals(ii))**2) * &
                            real( conjg(PsiKet(kk))*PsiKet(ll)*gMx(kk,ll) )

                        LDOSMy(xx,yy,zz,ee) = LDOSMy(xx,yy,zz,ee) + exp(-1D3*(Eaxis(ee)-eigvals(ii))**2) * &
                            real( conjg(PsiKet(kk))*PsiKet(ll)*gMy(kk,ll) )

                        LDOSMz(xx,yy,zz,ee) = LDOSMz(xx,yy,zz,ee) + exp(-1D3*(Eaxis(ee)-eigvals(ii))**2) * &
                            real( conjg(PsiKet(kk))*PsiKet(ll)*gMz(kk,ll) )
                            
                            
                        LDOSPx(xx,yy,zz,ee) = LDOSPx(xx,yy,zz,ee) + exp(-1D3*(Eaxis(ee)-eigvals(ii))**2) * &
                            real( conjg(PsiKet(kk))*PsiKet(ll)*gPx(kk,ll) )

                        LDOSPy(xx,yy,zz,ee) = LDOSPy(xx,yy,zz,ee) + exp(-1D3*(Eaxis(ee)-eigvals(ii))**2) * &
                            real( conjg(PsiKet(kk))*PsiKet(ll)*gPy(kk,ll) )

                        LDOSPz(xx,yy,zz,ee) = LDOSPz(xx,yy,zz,ee) + exp(-1D3*(Eaxis(ee)-eigvals(ii))**2) * &
                            real( conjg(PsiKet(kk))*PsiKet(ll)*gPz(kk,ll) )
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    DEALLOCATE(PsiKet);
    !!!!!!!!!!Spin-resolved LDOS and dispersion relations

    BStruc(xx, yy, 1:5) = eigvals(NN/2+1:NN/2+5) !Only lowest-energy eigenvalues > 0

    !Ground state energy at each Kx, Ky, z-point
    DO ii = 1, NN/2
        EGnd(xx,yy) = EGnd(xx,yy) + eigvals(ii)
    ENDDO
    DO zz = 1, Npz
        IF (zz <= NpzSC) THEN
            EGnd(xx,yy) = EGnd(xx,yy) - 4D0*muSC
        ELSE
            EGnd(xx,yy) = EGnd(xx,yy) - 4D0*muTI
        ENDIF
    ENDDO
    EGnd0(xx,yy) = EGnd(xx,yy)/Npz  !EGnd to zeroth order, before MFT subtraction

    !Subtract MFT component to get ground state energy
    DO zz = 1, Npz
        IF (zz <= NpzSC) THEN
            EGnd(xx,yy) = EGnd(xx,yy) - abs(UIntSC)* &
                   ( 2D0*abs(DeltaAS(xx,yy,zz))**2 + 2D0*abs(DeltaBS(xx,yy,zz))**2 )
        ELSE
            EGnd(xx,yy) = EGnd(xx,yy) - abs(UIntTI)* &
                   ( 2D0*abs(DeltaAS(xx,yy,zz))**2    + 2D0*abs(DeltaBS(xx,yy,zz))**2  &
                       + abs(DeltaAPUp(xx,yy,zz))**2  + abs(DeltaAPDn(xx,yy,zz))**2 &
                       + abs(DeltaBPUp(xx,yy,zz))**2  + abs(DeltaBPDn(xx,yy,zz))**2 )
        ENDIF
    ENDDO
    EGnd(xx,yy) = EGnd(xx,yy)/Npz
END SUBROUTINE CalculateFinal
