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
