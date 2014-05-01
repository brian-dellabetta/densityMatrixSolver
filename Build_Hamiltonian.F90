SUBROUTINE Build_Hamiltonian
    USE NEGF_variables
    IMPLICIT NONE

    IF (MatType == 0) THEN
        MM = -265D-2
        CALL BuildTICells
    ELSEIF (MatType == 1) THEN
        CALL BuildLiangTCICells
    ELSEIF (MatType == 2) THEN
        CALL BuildChenTCICells
    ELSE
        CALL BuildBiSeCells
    ENDIF
    CALL BuildSCCells

    !!!!!!!!!On-diagonal (kx/ky-direction) Interaction (1D)
    DO zz = 1, Npz-1
        index = (zz-1)*Norb
        IF (zz < NpzSC) THEN
            Htot(index+1:index+Norb,index+1:index+Norb) = AASC
            Htot(index+1:index+Norb,index+Norb+1:index+Norb*2) = HopZSC
            Htot(index+Norb+1:index+Norb*2,index+1:index+Norb) = conjg(transpose(HopZSC))
        ELSEIF (zz > NpzSC) THEN
            Htot(index+1:index+Norb,index+1:index+Norb) = AATI
            Htot(index+1:index+Norb,index+Norb+1:index+Norb*2) = HopZTI
            Htot(index+Norb+1:index+Norb*2,index+1:index+Norb) = conjg(transpose(HopZTI))
        ELSE
            Htot(index+1:index+Norb,index+1:index+Norb) = AASC
            Htot(index+1:index+Norb,index+Norb+1:index+Norb*2) = HopZSC !(HopZSC+HopZTI)/2D0
            Htot(index+Norb+1:index+Norb*2,index+1:index+Norb) = conjg(transpose(HopZSC)) !conjg(transpose((HopZSC+HopZTI)/2D0))
        ENDIF
    ENDDO
    DO zz = NpzSC+1, Npz-2
        index = (zz-1)*Norb
        Htot(index+1:index+Norb,index+2*Norb+1:index+Norb*3) = HopZZTI
        Htot(index+2*Norb+1:index+Norb*3,index+1:index+Norb) = conjg(transpose(HopZZTI))
    ENDDO

    IF (NpzSC == Npz) THEN
        Htot(NN-Norb+1:NN,NN-Norb+1:NN) = AASC
    ELSE
        Htot(NN-Norb+1:NN,NN-Norb+1:NN) = AATI
    ENDIF

    !!!!!!!!!!!!Add attractive Hubbard-type interaction to BdG off-diagonal blocks
    DO zz = 1, Npz
        index = (zz-1)*Norb
        IF (zz <= NpzSC) THEN
            Htot(index+1:index+Norb/2, index+1+Norb/2:index+Norb) = &
                     -abs(UIntSC) * ( DeltaAS(xx,yy,zz) * gDAS + DeltaBS(xx,yy,zz) * gDBS )
        ELSE
            Htot(index+1:index+Norb/2, index+1+Norb/2:index+Norb) = &
                     -abs(UIntTI) * ( DeltaMat(xx,yy,zz,:,:) )
        ENDIF

        Htot(index+1+Norb/2:index+Norb, index+1:index+Norb/2) = conjg(transpose(Htot(index+1:index+Norb/2, index+1+Norb/2:index+Norb)))
    ENDDO
    !!!!!!!!!!!!Add BdG components to AATI off-diagonal blocks
END SUBROUTINE Build_Hamiltonian


!Topological Crystalline Insulator Hamiltonian, from arXiv 1308.2424 Eq 12
SUBROUTINE BuildLiangTCICells
    USE NEGF_variables
    IMPLICIT NONE
    !Reset matrices
    HOnSite=c0;HSOC=c0;HZeeman=c0;
    AATI=c0;HopZTI=c0;HopZZTI=c0

    !!!!!!!!!!On-site mass difference between Sn, Te atoms
    HOnSite(1:6,1:6)   =  LFm*id(1:6,1:6)
    HOnSite(7:12,7:12) = -LFm*id(1:6,1:6)
    !!!!!!!!!!On-site mass difference between Sn, Te atoms

    !!!!!!!!!!Spin orbit coupling, L \cdot s
    HSOC(1:6,1:6)   = LFlambda1 * HLdotS
    HSOC(7:12,7:12) = LFlambda2 * HLdotS
    !!!!!!!!!!Spin orbit coupling, L \cdot s

    !!!!!!!!!!Zeeman splitting
    HZeeman = Mx*gMx + My*gMy + Mz*gMz
    !!!!!!!!!!Zeeman splitting

    !!!Electron component to BdG Hamiltonian (see HHH note Eq. 29)
    CALL BuildTCIHopping(Kx, Ky)  !Get HHop(k)

    AATI(1:12,1:12)    = HOnSite + HHop + HSOC + HZeeman - muTI*id
    HopZTI(1:12,1:12)  = HHopz          !NN, NNN hopping in z-direction
    HopZZTI(1:12,1:12) = c0             !3NN hopping neglected


    !!!Hole component to BdG Hamiltonian becomes H(k) => -H*(-k) in k-space
    CALL BuildTCIHopping(-Kx, -Ky)  !Get HHop(-k)

    AATI(13:24,13:24)    = -conjg( HOnSite + HHop + HSOC + HZeeman - muTI*id )
    HopZTI(13:24,13:24)  = -conjg( HHopz )
    HopZZTI(13:24,13:24) = -conjg( c0 )
END SUBROUTINE BuildLiangTCICells


!!!!!!!!!!Hopping interaction for TCI
SUBROUTINE BuildTCIHopping(kxval, kyval)
    USE NEGF_variables
    IMPLICIT NONE
    REAL  :: kxval, kyval, ckx, cky, skx, sky

    HNN=c0; HNNN=c0; HHop=c0;
    HNNz=c0;HNNNz=c0;HHopz=c0;

    ckx = cos(a0*kxval)     !kxval   !(kxval+kyval)/dsqrt(2D0)
    cky = cos(a0*kyval)     !kyval   !(kxval-kyval)/dsqrt(2D0)
    skx = sin(a0*kxval)
    sky = sin(a0*kyval)

    !-2sin(a0*Kz) -> <z+1|H|z> = i = -<z-1|H|z>
    ! 2cos(a0*Kz) -> <z+1|H|z> = 1 =  <z-1|H|z>

    !Nearest neighbor hopping
    HNN(1,1) = ckx/4D0
    HNN(2,2) = cky/4D0

    !2cos(kz)
    HNNz(3,3) = c1/8D0

    !Next-Nearest neighbor hopping
    HNNN(1,1) =  ckx*cky/2D0
    HNNN(1,2) = -skx*sky/2D0
    HNNN(2,1) = -sky*skx/2D0
    HNNN(2,2) =  cky*ckx/2D0

    !cos(kz)*(cos(kx)+cos(ky))
    HNNNz(1,1) = ckx/4D0
    HNNNz(2,2) = cky/4D0
    HNNNz(3,3) = (ckx+cky)/4D0
    !-sin(kx)*sin(kz)
    HNNNz(1,3) = ci*skx/4D0
    HNNNz(3,1) = ci*skx/4D0
    !-sin(ky)*sin(kz)
    HNNNz(2,3) = ci*sky/4D0
    HNNNz(3,2) = ci*sky/4D0

    HNN(4:6,4:6) = HNN(1:3,1:3)
    HNNz(4:6,4:6) = HNNz(1:3,1:3)
    HNNN(4:6,4:6) = HNNN(1:3,1:3)
    HNNNz(4:6,4:6) = HNNNz(1:3,1:3)


    HHop(1:6,1:6)   = LFt11*HNNN
    HHop(7:12,7:12) = LFt22*HNNN
    HHop(1:6,7:12)  = LFt12*HNN
    HHop(7:12,1:6)  = LFt12*HNN

    HHopz(1:6,1:6)   = LFt11*HNNNz
    HHopz(7:12,7:12) = LFt22*HNNNz
    HHopz(1:6,7:12)  = LFt12*HNNz
    HHopz(7:12,1:6)  = LFt12*HNNz
END SUBROUTINE BuildTCIHopping
!!!!!!!!!!Hopping interaction

!!!!!!!!!!Chen's TCI Hamiltonian, assumes Norb=8
SUBROUTINE BuildChenTCICells
    USE NEGF_variables
    IMPLICIT NONE
    
    !!!!!!!!!!Zeeman splitting
    HZeeman = Mx*gMx + My*gMy + Mz*gMz
    !!!!!!!!!!Zeeman splitting

    !Topological Insulator on-site term
    AATI(1:4,1:4) = (CFm - CFt1*(cos(2D0*Kx)+cos(2D0*Ky))) * SigZ0 &
                  + (CFt2*sin(Kx)*cos(Ky)) * SigXX &
                  + (CFt2*sin(Ky)*cos(Kx)) * SigXy &
                  + HZeeman  &
                  - muTI*id

    !Hole component to BdG Hamiltonian becomes H(k) => -H*(-k) in k-space  (see HHH note Eq. 29)
    AATI(5:8,5:8) = -conjg( &
                    (CFm - CFt1*(cos(-2D0*Kx)+cos(-2D0*Ky))) * SigZ0 &
                  + (CFt2*sin(-Kx)*cos(-Ky)) * SigXX &
                  + (CFt2*sin(-Ky)*cos(-Kx)) * SigXy &
                  + HZeeman  &
                  - muTI*id )

    !Topological Insulator z-direction hopping
    HopZTI(1:4,1:4) = CFt2*sin(Kx) * SigXX &
                    + CFt2*sin(Ky) * SigXY &
                  -ci*CFt2*(cos(Kx)+cos(Ky)) * SigXZ
                  
    HopZTI(5:8,5:8) = -conjg( &
                      CFt2*sin(Kx) * SigXX &
                    + CFt2*sin(Ky) * SigXY &
                  -ci*CFt2*(cos(Kx)+cos(Ky)) * SigXZ )
                  
    HopZZTI(1:4,1:4) = -CFt1/2D0/a0 * SigZ0

    HopZZTI(5:8,5:8) = -conjg( -CFt1/2D0/a0 * SigZ0 )
    
END SUBROUTINE BuildChenTCICells
!!!!!!!!!!Chen's TCI Hamiltonian, assumes Norb=8


!Topological Insulator Hamiltonian Unit Cell, assumes Norb=8
SUBROUTINE BuildTICells
    USE NEGF_variables
    IMPLICIT NONE

    !!!!!!!!!!Zeeman splitting
    HZeeman = Mx*gMx + My*gMy + Mz*gMz
    !!!!!!!!!!Zeeman splitting

    !Topological Insulator on-site term
    AATI(1:4,1:4) = (MM + (cos(Kx)+cos(Ky)))*g0 &
                  + sin(Kx)*g1 &
                  + sin(Ky)*g2 &
                  + HZeeman  &
                  - muTI*id

    !Hole component to BdG Hamiltonian becomes H(k) => -H*(-k) in k-space  (see HHH note Eq. 29)
    AATI(5:8,5:8) = -conjg( &
                    (MM + (cos(-Kx)+cos(-Ky)))*g0 &
                  + sin(-Kx)*g1 &
                  + sin(-Ky)*g2 &
                  + HZeeman &
                  - muTI*id )

    !Topological Insulator z-direction hopping
    HopZTI(1:4,1:4) = 5D-1*(g0 + ci*g3)
    HopZTI(5:8,5:8) = -conjg( 5D-1*(g0 + ci*g3) )

    HopZZTI = c0
    
END SUBROUTINE BuildTICells


!Bi2Se3 Hamiltonian Unit Cell, assumes Norb=8
SUBROUTINE BuildBiSeCells
    USE NEGF_variables
    !Material Parameters Specific to Bi2Se3
    REAL(8), PARAMETER          :: mpA1 = 2.26                  !eV*Angstroms
    REAL(8), PARAMETER          :: mpA2 = 3.33                  !eV*Angstroms
    REAL(8), PARAMETER          :: mpC = -0.0083                !eV
    REAL(8), PARAMETER          :: mpD1 = 5.74                  !eV*Angstroms^2
    REAL(8), PARAMETER          :: mpD2 = 30.4                  !ev*Angstroms^2
    REAL(8), PARAMETER          :: mpMTI = 0.28                   !eV
    REAL(8), PARAMETER          :: mpB1 = 6.86                  !eV*Angstroms^2
    REAL(8), PARAMETER          :: mpB2 = 44.5                  !eV*Angstroms^2

    !Topological Insulator on-site term
    AATI(1:4,1:4) = ( mpC + (mpD1+2.0*mpD2)/(a0*a0) - mpD2/(a0*a0)*(cos(Kx*a0)+cos(Ky*a0)) )*id &
                  + ( mpMTI - (mpB1+2.0*mpB2)/(a0*a0) + mpB2/(a0*a0)*(cos(Kx*a0)+cos(Ky*a0)) )*g0 &
                  + mpA2/a0*sin(Kx*a0)*g1 &
                  + mpA2/a0*sin(Ky*a0)*g2 &
                  + Mx*gMx + My*gMy + Mz*gMz  &
                  - muTI*id

    !Hole component to BdG Hamiltonian becomes H(k) => -H*(-k) in k-space  (see HHH note Eq. 29)
    AATI(5:8,5:8) = -conjg( &
                    ( mpC + (mpD1+2.0*mpD2)/(a0*a0) - mpD2/(a0*a0)*(cos(-Kx*a0)+cos(-Ky*a0)) )*id &
                  + ( mpMTI - (mpB1+2.0*mpB2)/(a0*a0) + mpB2/(a0*a0)*(cos(-Kx*a0)+cos(-Ky*a0)) )*g0 &
                  + mpA2/a0*sin(-Kx*a0)*g1 &
                  + mpA2/a0*sin(-Ky*a0)*g2 &
                  + Mx*gMx + My*gMy + Mz*gMz ) &
                  + muTI*id

    !Topological Insulator z-direction hopping
    HopZTI(1:4,1:4) = 5D-1*(mpB1/(a0*a0)*g0 - mpD1/(a0*a0)*id - ci*mpA1/a0*g3)
    HopZTI(5:8,5:8) = -conjg(HopZTI(1:4,1:4))

    HopZZTI = c0
END SUBROUTINE BuildBiSeCells


!Hamiltonain for model s-wave superconductor, assumes Norb=8
SUBROUTINE BuildSCCells
    USE NEGF_variables
    IMPLICIT NONE

    !S-Wave Supercondcutor on-site term !k^2 approx 1-cos(k)
    AASC(1:Norb/2,1:Norb/2) = (c1 - cos(Kx) &
                   + c1 - cos(Ky)) * id &
                   - muSC*id

    !Hole component to BdG Hamiltonian becomes H(k) => -H*(-k) in k-space  (see HHH note Eq. 29)
    AASC(1+Norb/2:Norb,1+Norb/2:Norb) = -conjg( &
                     c1 - cos(-Kx) &
                   + c1 - cos(-Ky) ) * id &
                   + muSC*id

    IF (fullZeeman == 1) THEN
        AASC(1:Norb/2,1:Norb/2) = AASC(1:Norb/2,1:Norb/2) + Mx*gMx + My*gMy + Mz*gMz
        AASC(1+Norb/2:Norb,1+Norb/2:Norb) = AASC(1+Norb/2:Norb,1+Norb/2:Norb) - conjg(Mx*gMx + My*gMy + Mz*gMz)
    ENDIF

    !Topological Insulator z-direction hopping
    HopZSC(1:Norb/2,1:Norb/2) = tSC*id
    HopZSC(1+Norb/2:Norb,1+Norb/2:Norb) = -conjg(tSC*id)

    !Hopping between SC and TI slices in z, connect all orbitals of same spin together?
    HopZSCTI = HopZSC !+ HopZTI
END SUBROUTINE BuildSCCells
