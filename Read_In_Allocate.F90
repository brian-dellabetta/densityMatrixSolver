SUBROUTINE Read_In_Allocate
    USE NEGF_variables
    IMPLICIT NONE
    REAL :: Kxmin, Kxmax, Kymin, Kymax, Emin, Emax

    !Read in parameters
    nargs = iargc()
    IF (nargs == 14) THEN
        CALL getarg(1, arg); READ (arg, *) Nprocs
        CALL InitializeMPI
        CALL getarg(2, arg); READ (arg, *) NpKx
        CALL getarg(3, arg); READ (arg, *) NpKy
        CALL getarg(4, arg); READ (arg, *) Npz
        CALL getarg(5, arg); READ (arg, *) NpzSC
        CALL getarg(6, arg); READ (arg, *) UIntTI
        CALL getarg(7, arg); READ (arg, *) UIntSC
        CALL getarg(8, arg); READ (arg, *) muTI
        CALL getarg(9, arg); READ (arg, *) muSC
        CALL getarg(10, arg); READ (arg, *) MatType
        CALL getarg(11, arg); READ (arg, *) Mx
        CALL getarg(12, arg); READ (arg, *) My
        CALL getarg(13, arg); READ (arg, *) Mz
        CALL getarg(14, outputfiledir)
    ELSE                    !default Values
        WRITE(*,*) "Insufficient command line arguments"
        WRITE(*,*) "nargs = ", nargs, "Using default values..."
        Nprocs = 1
        CALL InitializeMPI

        NpKx = 10
        NpKy = 1
        Npz = 10
        NpzSC = 0
        muTI = 0D0
        muSC = 0D0 !-5D-1
        UIntTI = 0D0
        UIntSC = 0D0
        MatType = 0
        Mx = 0D0
        My = 0D0
        Mz = 0D0
        outputfiledir = 'Debug/'
    ENDIF
    UIntTI = abs(UIntTI)
    UIntSC = abs(UIntSC)
    IF (MatType == 1) THEN
        Norb = 24
    ELSE
        Norb = 8
    ENDIF
    NN  = Npz*Norb

    !!!!!Variables to Set
    a0 = 1D0
    tSC = c1
    DeltaSC = 1D-3
    iteTot = 1000
    Accuracy = 1D-6
    T0 = 1D-6    !degrees Kelvin
    NpE = 1 !41
    fullDeltaMat = 0
    fullZeeman = 0 !if 1, include Zeeman in SC as well as TI

    IF (NpzSC > Npz) THEN
        NpzSC = Npz/2
    ENDIF

    !Allocate arrays...this is done to prevent stack overflow (automatic arrays are placed on the stack)
    ALLOCATE(Htot(NN,NN))            ;Htot = c0
    ALLOCATE(HopZTI(Norb,Norb))      ;HopZTI = c0
    ALLOCATE(HopZZTI(Norb,Norb))     ;HopZZTI = c0
    ALLOCATE(AATI(Norb,Norb))        ;AATI=c0

    ALLOCATE(AASC(Norb,Norb))        ;AASC=c0
    ALLOCATE(HopZSC(Norb,Norb))      ;HopZSC = c0
    ALLOCATE(HopZSCTI(Norb,Norb))    ;HopZSCTI = c0

    ALLOCATE(eigvecs(NN,NN))         ;eigvecs = c0
    ALLOCATE(eigvals(NN))            ;eigvals = 0D0
    ALLOCATE(BStruc(NpKx,NpKy,5))    ;BStruc = 0D0

    ALLOCATE(LDOS(NpKx,NpKy,Npz,NpE))   ;LDOS = 0D0
    ALLOCATE(LDOSMx(NpKx,NpKy,Npz,NpE)) ;LDOSMx = 0D0
    ALLOCATE(LDOSMy(NpKx,NpKy,Npz,NpE)) ;LDOSMy = 0D0
    ALLOCATE(LDOSMz(NpKx,NpKy,Npz,NpE)) ;LDOSMz = 0D0
    ALLOCATE(LDOSPx(NpKx,NpKy,Npz,NpE)) ;LDOSPx = 0D0
    ALLOCATE(LDOSPy(NpKx,NpKy,Npz,NpE)) ;LDOSPy = 0D0
    ALLOCATE(LDOSPz(NpKx,NpKy,Npz,NpE)) ;LDOSPz = 0D0
    
    ALLOCATE(KxAxis(NpKx))           ;KxAxis = 0D0
    ALLOCATE(KyAxis(NpKy))           ;KyAxis = 0D0
    ALLOCATE(Eaxis(NpE))             ;Eaxis = 0D0

    ALLOCATE(EGnd(NpKx,NpKy))        ;EGnd = 0D0
    ALLOCATE(EGnd0(NpKx,NpKy))       ;EGnd0 = 0D0

    !Pauli matrices
    ALLOCATE(Sig0(2, 2))   ;Sig0(1,:) = (/c1, c0/); Sig0(2,:) = (/c0, c1/)
    ALLOCATE(Sig1(2, 2))   ;Sig1(1,:) = (/c0, c1/); Sig1(2,:) = (/c1, c0/)
    ALLOCATE(Sig2(2, 2))   ;Sig2(1,:) = (/c0,-ci/); Sig2(2,:) = (/ci, c0/)
    ALLOCATE(Sig3(2, 2))   ;Sig3(1,:) = (/c1, c0/); Sig3(2,:) = (/c0,-c1/)

    ALLOCATE(HZeeman(Norb/2,Norb/2))  ;HZeeman=c0
    ALLOCATE(gMx(Norb/2,Norb/2))      ;gMx=c0
    ALLOCATE(gMy(Norb/2,Norb/2))      ;gMy=c0
    ALLOCATE(gMz(Norb/2,Norb/2))      ;gMz=c0
    ALLOCATE(gPx(Norb/2,Norb/2))      ;gPx=c0
    ALLOCATE(gPy(Norb/2,Norb/2))      ;gPy=c0
    ALLOCATE(gPz(Norb/2,Norb/2))      ;gPz=c0
    ALLOCATE(gDAS(Norb/2, Norb/2))    ;gDAS=c0
    ALLOCATE(gDBS(Norb/2, Norb/2))    ;gDBS=c0
    ALLOCATE(gDABS(Norb/2, Norb/2))   ;gDABS=c0
    ALLOCATE(gDAPUp(Norb/2, Norb/2))  ;gDAPUp=c0
    ALLOCATE(gDAPDn(Norb/2, Norb/2))  ;gDAPDn=c0
    ALLOCATE(gDBPUp(Norb/2, Norb/2))  ;gDBPUp=c0
    ALLOCATE(gDBPDn(Norb/2, Norb/2))  ;gDBPDn=c0
    ALLOCATE(gDABPUp(Norb/2, Norb/2)) ;gDBPUp=c0
    ALLOCATE(gDABPDn(Norb/2, Norb/2)) ;gDBPDn=c0

    !Identity Matrix
    ALLOCATE(id(Norb/2,Norb/2))       ;id=c0
    DO ii = 1, Norb/2
        id(ii,ii) = c1
    ENDDO

    IF (Norb == 8) THEN
        CALL Allocate8OrbCells
    ELSEIF (Norb == 24) THEN
        CALL Allocate24OrbCells
    ENDIF

    ALLOCATE(DeltaAS(NpKx,NpKy,Npz));    DeltaAS=c0
    ALLOCATE(DeltaBS(NpKx,NpKy,Npz));    DeltaBS=c0
    ALLOCATE(DeltaABS(NpKx,NpKy,Npz));   DeltaABS=c0
    ALLOCATE(DeltaAPUp(NpKx,NpKy,Npz));  DeltaAPUp=c0
    ALLOCATE(DeltaBPUp(NpKx,NpKy,Npz));  DeltaBPUp=c0
    ALLOCATE(DeltaAPDn(NpKx,NpKy,Npz));  DeltaAPDn=c0
    ALLOCATE(DeltaBPDn(NpKx,NpKy,Npz));  DeltaBPDn=c0
    ALLOCATE(DeltaABPUp(NpKx,NpKy,Npz)); DeltaABPUp=c0
    ALLOCATE(DeltaABPDn(NpKx,NpKy,Npz)); DeltaABPDn=c0

    ALLOCATE(DeltaMat(NpKx,NpKy,Npz,Norb/2,Norb/2)); DeltaMat=c0
    ALLOCATE(DeltaMatOld(Npz,Norb/2,Norb/2)); DeltaMatOld=c0

    !Set up momentum range.  Full BZ is (-pi/a0 < kx,ky < pi/a0)
    IF (NpKx == 1) THEN
        KxAxis(1) = 0D0
    ELSE
        IF (MatType == 2) THEN !TCI
            Kxmin = 0     !-pi/(2D0*a0) !
            Kxmax = pi      !pi/(2D0*a0) !
        ELSE
            Kxmin = -pi/(2D0*a0)
            Kxmax = pi/(2D0*a0)
        ENDIF
        deltaKx = (Kxmax-Kxmin)/(NpKx-1);
        DO xx = 1, NpKx
            KxAxis(xx) = Kxmin + deltaKx*(xx-1)
        ENDDO
    ENDIF

    IF (NpKy == 1) THEN
        KyAxis(1) = 0D0
        KyMid = 1
    ELSE
        IF (MatType == 2) THEN !TCI
            Kymin = 0     !-pi/(2D0*a0) !
            Kymax = pi      !pi/(2D0*a0) !
        ELSE
            Kymin = -pi/(2D0*a0)
            Kymax = pi/(2D0*a0)
        ENDIF
        deltaKy = (Kymax-Kymin)/(NpKy-1);
        DO yy = 1, NpKy
            KyAxis(yy) = Kymin + deltaKy*(yy-1)
        ENDDO

        KyMid = int(NpKy/2)
    ENDIF

    !Set up energy range for LDOS(kx,ky,z,E)
    IF (NpE == 1) THEN
        Eaxis(1) = 0D0
    ELSE
        Emin = -1D0
        Emax = 1D0
        deltaE = (Emax-Emin)/(NpE-1);
        DO ee = 1, NpE
            Eaxis(ee) = Emin + deltaE*(ee-1)
        ENDDO
    ENDIF

END SUBROUTINE Read_In_Allocate

SUBROUTINE Allocate8OrbCells
    USE NEGF_variables
    IMPLICIT NONE

    !Gamma matrices for TI
    ALLOCATE(g1(Norb/2, Norb/2))      ;g1=c0
    ALLOCATE(g2(Norb/2, Norb/2))      ;g2=c0
    ALLOCATE(g3(Norb/2, Norb/2))      ;g3=c0
    ALLOCATE(g0(Norb/2, Norb/2))      ;g0=c0

    ALLOCATE(SigZ0(Norb/2, Norb/2))   ;SigZ0=c0
    ALLOCATE(SigXX(Norb/2, Norb/2))   ;SigXX=c0
    ALLOCATE(SigXY(Norb/2, Norb/2))   ;SigXY=c0
    ALLOCATE(SigXZ(Norb/2, Norb/2))   ;SigXZ=c0
    
    !!!!!!!!!!!!!!!Chen's TI Model
    g1(1,:) = (/ c0, c0, c0,-ci/)
    g1(2,:) = (/ c0, c0, ci, c0/)
    g1(3,:) = (/ c0,-ci, c0, c0/)
    g1(4,:) = (/ ci, c0, c0, c0/)

    g2(1,:) = (/ c0, c0, c0,-c1/)
    g2(2,:) = (/ c0, c0,-c1, c0/)
    g2(3,:) = (/ c0,-c1, c0, c0/)
    g2(4,:) = (/-c1, c0, c0, c0/)

    g3(1,:) = (/c0, c0,-ci, c0/)
    g3(2,:) = (/c0, c0, c0,-ci/)
    g3(3,:) = (/ci, c0, c0, c0/)
    g3(4,:) = (/c0, ci, c0, c0/)

    g0(1,:) = (/c1, c0, c0, c0/)
    g0(2,:) = (/c0, c1, c0, c0/)
    g0(3,:) = (/c0, c0,-c1, c0/)
    g0(4,:) = (/c0, c0, c0,-c1/)
    !!!!!!!!!!!!!!!Chen's TI Model

    !!!!!!!!!!!!!!!Chen's TCI Model
    SigZ0(1,:) = (/ c1,  c0,  c0,  c0/)
    SigZ0(2,:) = (/ c0,  c1,  c0,  c0/)
    SigZ0(3,:) = (/ c0,  c0, -c1,  c0/)
    SigZ0(4,:) = (/ c0,  c0,  c0, -c1/)

    SigXX(1,:) = (/ c0,  c0,  c0,  c1/)
    SigXX(2,:) = (/ c0,  c0,  c1,  c0/)
    SigXX(3,:) = (/ c0,  c1,  c0,  c0/)
    SigXX(4,:) = (/ c1,  c0,  c0,  c0/)

    SigXY(1,:) = (/c0,  c0,  c0, -ci/)
    SigXY(2,:) = (/c0,  c0,  ci,  c0/)
    SigXY(3,:) = (/c0, -ci,  c0,  c0/)
    SigXY(4,:) = (/ci,  c0,  c0,  c0/)

    SigXZ(1,:) = (/c0,  c0,  c1,  c0/)
    SigXZ(2,:) = (/c0,  c0,  c0, -c1/)
    SigXZ(3,:) = (/c1,  c0,  c0,  c0/)
    SigXZ(4,:) = (/c0, -c1,  c0,  c0/)
    !!!!!!!!!!!!!!!Chen's TCI Model

    !!!!!!!!!!!!!!!Zeeman splitting fields
    gMx(1,:) = (/c0, c1, c0, c0/)
    gMx(2,:) = (/c1, c0, c0, c0/)
    gMx(3,:) = (/c0, c0, c0, c1/)
    gMx(4,:) = (/c0, c0, c1, c0/)

    gMy(1,:) = (/c0,-ci, c0, c0/)
    gMy(2,:) = (/ci, c0, c0, c0/)
    gMy(3,:) = (/c0, c0, c0,-ci/)
    gMy(4,:) = (/c0, c0, ci, c0/)

    gMz(1,:) = (/c1, c0, c0, c0/)
    gMz(2,:) = (/c0,-c1, c0, c0/)
    gMz(3,:) = (/c0, c0, c1, c0/)
    gMz(4,:) = (/c0, c0, c0,-c1/)
    !!!!!!!!!!!!!!!Zeeman splitting fields
    

    !!!!!!!!S-wave pairs k,up with (-k,down)
    gDAS(1,:) = (/ c0, c1, c0, c0/)
    gDAS(2,:) = (/-c1, c0, c0, c0/)
    gDAS(3,:) = (/ c0, c0, c0, c0/)
    gDAS(4,:) = (/ c0, c0, c0, c0/)

    gDBS(1,:) = (/ c0, c0, c0, c0/)
    gDBS(2,:) = (/ c0, c0, c0, c0/)
    gDBS(3,:) = (/ c0, c0, c0, c1/)
    gDBS(4,:) = (/ c0, c0,-c1, c0/)

    gDABS(1,:) = (/ c0, c0, c0, c1/)
    gDABS(2,:) = (/ c0, c0,-c1, c0/)
    gDABS(3,:) = (/ c0, c1, c0, c0/)
    gDABS(4,:) = (/-c1, c0, c0, c0/)
    !!!!!!!!S-wave pairs k,up with (-k,down)

    !!!!!!!!P-wave pairs k,up with (-k,up)
    gDAPUp(1,:) = (/ c1, c0, c0, c0/)
    gDAPUp(2,:) = (/ c0, c0, c0, c0/)
    gDAPUp(3,:) = (/ c0, c0, c0, c0/)
    gDAPUp(4,:) = (/ c0, c0, c0, c0/)

    gDAPDn(1,:) = (/ c0, c0, c0, c0/)
    gDAPDn(2,:) = (/ c0, c1, c0, c0/)
    gDAPDn(3,:) = (/ c0, c0, c0, c0/)
    gDAPDn(4,:) = (/ c0, c0, c0, c0/)

    gDBPUp(1,:) = (/ c0, c0, c0, c0/)
    gDBPUp(2,:) = (/ c0, c0, c0, c0/)
    gDBPUp(3,:) = (/ c0, c0, c1, c0/)
    gDBPUp(4,:) = (/ c0, c0, c0, c0/)

    gDBPDn(1,:) = (/ c0, c0, c0, c0/)
    gDBPDn(2,:) = (/ c0, c0, c0, c0/)
    gDBPDn(3,:) = (/ c0, c0, c0, c0/)
    gDBPDn(4,:) = (/ c0, c0, c0, c1/)

    gDABPUp(1,:) = (/ c0, c0, c1, c0/)
    gDABPUp(2,:) = (/ c0, c0, c0, c0/)
    gDABPUp(3,:) = (/-c1, c0, c0, c0/)
    gDABPUp(4,:) = (/ c0, c0, c0, c0/)

    gDABPDn(1,:) = (/ c0, c0, c0, c0/)
    gDABPDn(2,:) = (/ c0, c0, c0, c1/)
    gDABPDn(3,:) = (/ c0, c0, c0, c0/)
    gDABPDn(4,:) = (/ c0,-c1, c0, c0/)
    !!!!!!!!P-wave pairs k,up with (-k,up)
END SUBROUTINE Allocate8OrbCells

SUBROUTINE Allocate24OrbCells
    USE NEGF_variables
    IMPLICIT NONE

    !!!!!!!!SOC Hamiltonian
    ALLOCATE(HLdotS(6,6))
    HLdotS(1,:) = (/ c0, ci, c0, c0, c0,-c1/)
    HLdotS(2,:) = (/-ci, c0, c0, c0, c0, ci/)
    HLdotS(3,:) = (/ c0, c0, c0, c1,-ci, c0/)
    HLdotS(4,:) = (/ c0, c0, c1, c0,-ci, c0/)
    HLdotS(5,:) = (/ c0, c0, ci, ci, c0, c0/)
    HLdotS(6,:) = (/-c1,-ci, c0, c0, c0, c0/)

    ALLOCATE(HSOC(12,12)); HSOC=c0
    ALLOCATE(HOnSite(12,12)); HOnSite=c0

    ALLOCATE(HNN(6,6)); HNN=c0
    ALLOCATE(HNNz(6,6)); HNNz=c0
    ALLOCATE(HNNN(6,6)); HNNN=c0
    ALLOCATE(HNNNz(6,6)); HNNNz=c0

    ALLOCATE(HHop(12,12)); HHop=c0
    ALLOCATE(HHopz(12,12)); HHopz=c0


    !!!!!!!!!!!!!!!Zeeman splitting fields
    gMx(1,1:6) = (/c0, c0, c0, c1, c0, c0/)
    gMx(2,1:6) = (/c0, c0, c0, c0, c1, c0/)
    gMx(3,1:6) = (/c0, c0, c0, c0, c0, c1/)
    gMx(4,1:6) = (/c1, c0, c0, c0, c0, c0/)
    gMx(5,1:6) = (/c0, c1, c0, c0, c0, c0/)
    gMx(6,1:6) = (/c0, c0, c1, c0, c0, c0/)
    gMx(7:12,7:12) = gMx(1:6,1:6)

    gMy(1,1:6) = (/c0, c0, c0,-ci, c0, c0/)
    gMy(2,1:6) = (/c0, c0, c0, c0,-ci, c0/)
    gMy(3,1:6) = (/c0, c0, c0, c0, c0,-ci/)
    gMy(4,1:6) = (/ci, c0, c0, c0, c0, c0/)
    gMy(5,1:6) = (/c0, ci, c0, c0, c0, c0/)
    gMy(6,1:6) = (/c0, c0, ci, c0, c0, c0/)
    gMy(7:12,7:12) = gMy(1:6,1:6)

    gMz(1,1:6) = (/c1, c0, c0, c0, c0, c0/)
    gMz(2,1:6) = (/c0, c1, c0, c0, c0, c0/)
    gMz(3,1:6) = (/c0, c0, c1, c0, c0, c0/)
    gMz(4,1:6) = (/c0, c0, c0,-c1, c0, c0/)
    gMz(5,1:6) = (/c0, c0, c0, c0,-c1, c0/)
    gMz(6,1:6) = (/c0, c0, c0, c0, c0,-c1/)
    gMz(7:12,7:12) = gMz(1:6,1:6)
    !!!!!!!!!!!!!!!Zeeman splitting fields

    !!!!!!!!!!!!!!!Orbital-resolved weights
    gPx(1,1) = c1; gPx(4,4) = c1; gPx(7,7) = c1; gPx(10,10) = c1
    gPy(2,2) = c1; gPy(5,5) = c1; gPy(8,8) = c1; gPy(11,11) = c1
    gPz(3,3) = c1; gPz(6,6) = c1; gPz(9,9) = c1; gPz(12,12) = c1
    !!!!!!!!!!!!!!!Orbital-resolved weights

    !!!!!!!!S-wave pairs k,up with (-k,down)
    gDAS(1,1:6) = (/ c0, c0, c0, c1, c0, c0/)
    gDAS(2,1:6) = (/ c0, c0, c0, c0, c1, c0/)
    gDAS(3,1:6) = (/ c0, c0, c0, c0, c0, c1/)
    gDAS(4,1:6) = (/-c1, c0, c0, c0, c0, c0/)
    gDAS(5,1:6) = (/ c0,-c1, c0, c0, c0, c0/)
    gDAS(6,1:6) = (/ c0, c0,-c1, c0, c0, c0/)

    gDBS(7:12,7:12) =  gDAS(1:6,1:6)

    gDABS(1:6,7:12) =  gDAS(1:6,1:6)
    gDABS(7:12,1:6) = -gDAS(1:6,1:6)
    !!!!!!!!S-wave pairs k,up with (-k,down)

    !!!!!!!!P-wave pairs k,up with (-k,up)
    gDAPUp(1,1:3) = (/c1, c0, c0/)
    gDAPUp(2,1:3) = (/c0, c1, c0/)
    gDAPUp(3,1:3) = (/c0, c0, c1/)

    gDAPDn(4:6,4:6) = gDAPUp(1:3,1:3)

    gDBPUp(7:9,7:9) = gDAPUp(1:3,1:3)

    gDBPDn(10:12,10:12) = gDAPUp(1:3,1:3)

    gDABPUp(1:3,7:9) = gDAPUp(1:3,1:3)
    gDABPUp(7:9,1:3) = -gDAPUp(1:3,1:3)

    gDABPDn(4:6,10:12) = gDAPUp(1:3,1:3)
    gDABPDn(10:12,4:6) = -gDAPUp(1:3,1:3)
    !!!!!!!!P-wave pairs k,up with (-k,up)
END SUBROUTINE Allocate24OrbCells










    !!!!!!!!!!!!!!!!Original Model Gamma matrices
!g1(1,:) = (/c0,  c0,  c0,  c1/)
!g1(2,:) = (/c0,  c0,  c1,  c0/)
!g1(3,:) = (/c0,  c1,  c0,  c0/)
!g1(4,:) = (/c1,  c0,  c0,  c0/)
!
!g2(1,:) = (/c0,  c0,  c0, -ci/)
!g2(2,:) = (/c0,  c0,  ci,  c0/)
!g2(3,:) = (/c0, -ci,  c0,  c0/)
!g2(4,:) = (/ci,  c0,  c0,  c0/)
!
!g3(1,:) = (/c0,  c0,  c1,  c0/)
!g3(2,:) = (/c0,  c0,  c0, -c1/)
!g3(3,:) = (/c1,  c0,  c0,  c0/)
!g3(4,:) = (/c0, -c1,  c0,  c0/)
!
!g0(1,:) = (/c1,  c0,  c0,  c0/)
!g0(2,:) = (/c0,  c1,  c0,  c0/)
!g0(3,:) = (/c0,  c0, -c1,  c0/)
!g0(4,:) = (/c0,  c0,  c0, -c1/)
!
!
!gMx(1,:) = (/c0,  c1,  c0,  c0/)
!gMx(2,:) = (/c1,  c0,  c0,  c0/)
!gMx(3,:) = (/c0,  c0,  c0, -c1/)
!gMx(4,:) = (/c0,  c0, -c1,  c0/)
!
!gMy(1,:) = (/c0, -ci,  c0,  c0/)
!gMy(2,:) = (/ci,  c0,  c0,  c0/)
!gMy(3,:) = (/c0,  c0,  c0,  ci/)
!gMy(4,:) = (/c0,  c0, -ci,  c0/)
!
!gMz(1,:) = (/c1,  c0,  c0,  c0/)
!gMz(2,:) = (/c0, -c1,  c0,  c0/)
!gMz(3,:) = (/c0,  c0, -c1,  c0/)
!gMz(4,:) = (/c0,  c0,  c0,  c1/)
    !!!!!!!!!!!!!!!Original Model Gamma matrices


    !!!!!!!!!!!!!!!BHZ Model
!g1(1,:) = (/c0,  c0,  c1,  c0/)
!g1(2,:) = (/c0,  c0,  c0, -c1/)
!g1(3,:) = (/c1,  c0,  c0,  c0/)
!g1(4,:) = (/c0, -c1,  c0,  c0/)
!
!g2(1,:) = (/c0,  c0,  c0, -ci/)
!g2(2,:) = (/c0,  c0,  ci,  c0/)
!g2(3,:) = (/c0, -ci,  c0,  c0/)
!g2(4,:) = (/ci,  c0,  c0,  c0/)
!
!g3(1,:) = (/c0,  c0, -ci,  c0/)
!g3(2,:) = (/c0,  c0,  c0, -ci/)
!g3(3,:) = (/ci,  c0,  c0,  c0/)
!g3(4,:) = (/c0,  ci,  c0,  c0/)
!
!g0(1,:) = (/c1,  c0,  c0,  c0/)
!g0(2,:) = (/c0,  c1,  c0,  c0/)
!g0(3,:) = (/c0,  c0, -c1,  c0/)
!g0(4,:) = (/c0,  c0,  c0, -c1/)
!
!
!gMx(1,:) = (/c0,  c1,  c0,  c0/)
!gMx(2,:) = (/c1,  c0,  c0,  c0/)
!gMx(3,:) = (/c0,  c0,  c0,  c1/)
!gMx(4,:) = (/c0,  c0,  c1,  c0/)
!
!gMy(1,:) = (/c0, -ci,  c0,  c0/)
!gMy(2,:) = (/ci,  c0,  c0,  c0/)
!gMy(3,:) = (/c0,  c0,  c0, -ci/)
!gMy(4,:) = (/c0,  c0,  ci,  c0/)
!
!gMz(1,:) = (/c1,  c0,  c0,  c0/)
!gMz(2,:) = (/c0, -c1,  c0,  c0/)
!gMz(3,:) = (/c0,  c0,  c1,  c0/)
!gMz(4,:) = (/c0,  c0,  c0, -c1/)
    !!!!!!!!!!!!!!!BHZ Model
