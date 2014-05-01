MODULE  NEGF_variables
    IMPLICIT NONE

    !constants
    COMPLEX, PARAMETER :: ci = cmplx(0D0,1D0), c0 = cmplx(0D0,0D0), c1 = cmplx(1D0,0D0)
    REAL, PARAMETER    :: q = 1.602D-19, m0 = 9.1019D-31, pi = 4D0*atan(1D0)
    REAL, PARAMETER    :: eps0 = 8.85418782D-12, kB = 1.38066D-23
    REAL, PARAMETER    :: h = 6.626D-34, hbar = h/(2D0*pi)

    !Device parameters for TI
    REAL                     :: MM                           !Dirac Mass, MM = -1.5 (-4.5) sets the topologically nontrivial (trivial) phase
    INTEGER                  :: Norb                         !4 orbitals for 3D dirac Hamiltonian
    REAL                     :: muTI, muSC                   !Bulk chemical potential for TI and SC
    REAL                     :: UIntTI, UIntSC               !Attractive Mean-Field interaction for TI half and SC half
    COMPLEX                  :: tSC                          !Hopping for S-wave superconductor
    REAL                     :: T0                           !Temperature must be nonzero to avoid intrinsic TI superconductivity
    REAL                     :: Mx, My, Mz                   !Zeeman field added to TI Hamiltonian (uniform, homogeneous)
    REAL                     :: Myz                          !yz Mirror-plane breaking perturbation
    INTEGER                  :: MatType                      !0=Model TI, 1=Liang Fu TCI model, 2=Bi2Se3
    REAL                     :: a0                           !Lattice Constant (Angstroms)

    !Tight-binding parameters for LF model TCI SnTe, from Nat Comm 3, 982(2013)
    REAL, PARAMETER :: LFm   =  0.8;
    REAL, PARAMETER :: LFt11 =  0.5;
    REAL, PARAMETER :: LFt12 =  1.8;
    REAL, PARAMETER :: LFt22 = -0.5;
    REAL, PARAMETER :: LFlambda1 = 0.5;
    REAL, PARAMETER :: LFlambda2 = 0.5;
    
    !Tight-binding parameters for CF model TCI, from arXiv 1308.2424
    REAL, PARAMETER :: CFm  =  2.5;
    REAL, PARAMETER :: CFt1 = -1.0;
    REAL, PARAMETER :: CFt2 =  0.5;

    !input parameters from command line
    INTEGER             :: nargs
    CHARACTER           :: arg*200, outputfiledir*200
    INTEGER             :: NpKx, NpKy, Npz, NpzSC, NN            !number of points per layer

    !iteration/convergence variables
    INTEGER             :: ite, iteTot
    REAL                :: AmtChange  !Magnitude in change from eigvecs(ite-1) to eigvecs(ite)
    REAL                :: RMSChange  !Root Means Square change
    REAL                :: Accuracy   !Maximum allowable difference between iterations before convergence determined

    !gamma matrices
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: g1,g2,g3,g0,id

    !TCI matrices
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: Sig0,Sig1,Sig2,Sig3,SigZ0,SigXX,SigXY,SigXZ
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: HLdotS, HSOC, HOnSite, HZeeman, HNN, HNNN, HHop
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: HNNz, HNNNz, HHopz

    !Superconducting Parameters
    COMPLEX, DIMENSION(:,:), ALLOCATABLE   :: gDAS, gDBS, gDAPUp, gDBPUp, gDAPDn, gDBPDn, gDABS, gDABPUp, gDABPDn  !SC Pairing Gamma matrices (4x4)
    COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: DeltaAS, DeltaBS, DeltaAPUp, DeltaAPDn, DeltaBPUp, DeltaBPDn, DeltaABS, DeltaABPUp, DeltaABPDn    !\Delta_A(\vec{r})
    COMPLEX, ALLOCATABLE                   :: DeltaMat(:,:,:,:,:), DeltaMatOld(:,:,:)       !Entire pairing matrix for each (Kx,Ky,z)
    REAL                                   :: DeltaSC    !Initial condition for S-wave order parameter inside superconductor
    REAL, DIMENSION(:,:), ALLOCATABLE      :: EGnd, EGnd0 !Ground State Energy
    INTEGER                                :: fullDeltaMat, fullZeeman

    !Hamiltonian matrices
    COMPLEX, DIMENSION(:,:), ALLOCATABLE   :: Htot
    COMPLEX, DIMENSION(:,:), ALLOCATABLE   :: AATI, HopZTI, HopZZTI, AASC, HopZSC, HopZSCTI   !hopping matrices along x,y,z direction for TI and SC components
    COMPLEX, DIMENSION(:,:), ALLOCATABLE   :: eigvecs
    COMPLEX, DIMENSION(:,:), ALLOCATABLE   :: gMx, gMy, gMz  !Gamma matrices to apply Zeeman Field
    COMPLEX, DIMENSION(:,:), ALLOCATABLE   :: gPx, gPy, gPz  !To calculate Px,Py,Pz-orbital weights 
    REAL, DIMENSION(:), ALLOCATABLE        :: eigvals
    REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: LDOS, LDOSMx, LDOSMy, LDOSMz, LDOSPx, LDOSPy, LDOSPz
    REAL, DIMENSION(:,:,:), ALLOCATABLE    :: BStruc

    !MPI variables
    INTEGER    :: Nprocs
    INTEGER    :: mpirank, mpisize, mpiroot, mpierror
    INTEGER, ALLOCATABLE :: mpistatus(:)

    !Momentum parameters
    REAL                                 :: Kx, Ky, deltaKx, deltaKy, KxMid, KyMid
    REAL, DIMENSION(:), ALLOCATABLE      :: KxAxis, KyAxis

    !Energy parameters
    INTEGER                              :: NpE
    REAL                                 :: deltaE
    REAL, DIMENSION(:), ALLOCATABLE      :: Eaxis

    !temp/dummy variables
    INTEGER                              :: ii, jj, kk, ll, ee, ctr, xx, yy, zz, index

    CONTAINS
        FUNCTION Kron(A, B)
            IMPLICIT NONE
            COMPLEX, DIMENSION(:,:), INTENT(IN) :: A, B
            COMPLEX, DIMENSION(size(A,1)*size(B,1), size(A,2)*size(B,2)) :: Kron
            INTEGER :: mm, nn
            DO mm = 1, size(B,1)
                DO nn = 1, size(B,2)
                    Kron(1+(mm-1)*size(B,1):mm*size(B,1),1+(nn-1)*size(B,2):nn*size(B,2)) = A(mm,nn)*B
                ENDDO
            ENDDO
        END FUNCTION Kron
END MODULE NEGF_variables



