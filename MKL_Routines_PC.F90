!solve Htot*eigvecs = eigval*eigvecs for Hermitian matrix Htot
SUBROUTINE get_evalsevecs
    USE NEGF_variables
    IMPLICIT NONE
    INTEGER               :: LWORK, LRWORK, info
    COMPLEX, ALLOCATABLE  :: WORK(:)
    REAL, ALLOCATABLE     :: RWORK(:)
    
    LWORK = NN*NN + 2*NN
    LRWORK = 3*NN-2 !4*NN*NN
    
    ALLOCATE(WORK(LWORK));   WORK = c0
    ALLOCATE(RWORK(LRWORK)); RWORK = 0D0
    
    !CALL CHEEVD(JOBZ, UPLO, NN, eigvecs, NN, eigvals, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, info)
    CALL CHEEV('V', 'U', NN, eigvecs, NN, eigvals, WORK, LWORK, RWORK, info)
    
    IF (info /= 0) THEN
        WRITE(*,*) "Error in get_Eigvals computation, error in parameter", info
    ENDIF
    
    DEALLOCATE(WORK); DEALLOCATE(RWORK);
END SUBROUTINE get_evalsevecs

!dummy subroutines needed in cluster version of code
SUBROUTINE InitializeMPI
    USE NEGF_variables
    mpisize = 1
    mpirank = 0
    mpiroot = 0
END SUBROUTINE

SUBROUTINE BroadcastDataMPI
END SUBROUTINE

SUBROUTINE GatherDataMPI
END SUBROUTINE

SUBROUTINE FinalizeMPI
END SUBROUTINE

SUBROUTINE backup_evalsevecs
    USE NEGF_variables
    IMPLICIT NONE
    INTEGER     :: LWORK, LRWORK, LIWORK, info, numeigvals
    COMPLEX, DIMENSION(:), ALLOCATABLE :: WORK
    REAL, DIMENSION(:), ALLOCATABLE     :: RWORK
    INTEGER, DIMENSION(:), ALLOCATABLE     :: IWORK, isuppz
    REAL                                :: abstol
    
    abstol = 1D-6
    LWORK = 4*NN
    LRWORK = 48*NN
    LIWORK = 20*NN
    
    ALLOCATE(WORK(LWORK));   WORK = dcmplx(0D0,0D0)
    ALLOCATE(RWORK(LRWORK)); RWORK = 0D0
    ALLOCATE(IWORK(LIWORK)); IWORK = 0
    ALLOCATE(isuppz(2*NN));  isuppz = 0
    
    eigvecs=Htot; !Input Hamiltonian will be overwritten with eigenvectors
    CALL ZHEEVR('V', 'A', 'L', NN, Htot, NN, 0D0, 0D0, 0, 0, abstol, numeigvals, eigvals, eigvecs, NN, isuppz, work, LWORK, RWORK, LRWORK, IWORK, LIWORK, info)
    IF (info /= 0) THEN
        WRITE(*,*) "Error in get_Eigvals computation, error in parameter", info
    ENDIF
    WRITE(*,*) numeigvals, abstol, info
    WRITE(*,*) eigvals
    DEALLOCATE(WORK); DEALLOCATE(RWORK); DEALLOCATE(IWORK); DEALLOCATE(isuppz)
END SUBROUTINE backup_evalsevecs