SUBROUTINE InitializeMPI
    USE NEGF_variables
    USE MPI
    IMPLICIT NONE

    IF (Nprocs > 1) THEN
        CALL MPI_Init( mpierror ); CALL ErrorCheckMPI(mpierror, 1)
        CALL MPI_Comm_size( MPI_COMM_WORLD, mpisize, mpierror ); CALL ErrorCheckMPI(mpierror, 1)
        CALL MPI_Comm_Rank( MPI_COMM_WORLD, mpirank, mpierror ); CALL ErrorCheckMPI(mpierror, 1)
        ALLOCATE (mpistatus(MPI_STATUS_SIZE))
    ELSE
        mpisize = 1
        mpirank = 0
    ENDIF
    mpiroot = 0
END SUBROUTINE


SUBROUTINE BroadcastDataMPI
    USE NEGF_variables
    USE MPI
    IMPLICIT NONE

    IF (Nprocs > 1) THEN
        CALL MPI_BARRIER(MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 2)

        !CALL MPI_BCAST(Htot, NNz*NNz*Npx, MPI_COMPLEX16, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 2)
    ENDIF
END SUBROUTINE BroadcastDataMPI


SUBROUTINE GatherDataMPI
    USE NEGF_variables
    USE MPI
    IMPLICIT NONE

    REAL, ALLOCATABLE :: rmat1(:), rmat2(:,:), rmat3(:,:,:), rmat4(:,:,:,:)
    COMPLEX, ALLOCATABLE :: cmat1(:), cmat2(:,:), cmat3(:,:,:), cmat4(:,:,:,:), cmat5(:,:,:,:,:)

    IF (Nprocs > 1) THEN
        CALL MPI_BARRIER(MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)

        ALLOCATE(rmat4(NpKx, NpKy, Npz, NpE))
        rmat4 = LDOS
        CALL MPI_REDUCE(rmat4, LDOS, NpKx*NpKy*Npz*NpE, MPI_REAL, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        rmat4 = LDOSMx
        CALL MPI_REDUCE(rmat4, LDOSMx, NpKx*NpKy*Npz*NpE, MPI_REAL, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        rmat4 = LDOSMy
        CALL MPI_REDUCE(rmat4, LDOSMy, NpKx*NpKy*Npz*NpE, MPI_REAL, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        rmat4 = LDOSMz
        CALL MPI_REDUCE(rmat4, LDOSMz, NpKx*NpKy*Npz*NpE, MPI_REAL, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        rmat4 = LDOSPx
        CALL MPI_REDUCE(rmat4, LDOSPx, NpKx*NpKy*Npz*NpE, MPI_REAL, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        rmat4 = LDOSPy
        CALL MPI_REDUCE(rmat4, LDOSPy, NpKx*NpKy*Npz*NpE, MPI_REAL, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        rmat4 = LDOSPz
        CALL MPI_REDUCE(rmat4, LDOSPz, NpKx*NpKy*Npz*NpE, MPI_REAL, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        DEALLOCATE(rmat4)

        ALLOCATE(rmat2(NpKx, NpKy))
        rmat2 = EGnd
        CALL MPI_REDUCE(rmat2, EGnd, NpKx*NpKy, MPI_REAL, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        rmat2 = EGnd0
        CALL MPI_REDUCE(rmat2, EGnd0, NpKx*NpKy, MPI_REAL, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        DEALLOCATE(rmat2)

        ALLOCATE(cmat3(NpKx, NpKy, Npz))
        cmat3 = DeltaAS
        CALL MPI_REDUCE(cmat3, DeltaAS,   NpKx*NpKy*Npz, MPI_COMPLEX, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        cmat3 = DeltaBS
        CALL MPI_REDUCE(cmat3, DeltaBS,   NpKx*NpKy*Npz, MPI_COMPLEX, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        cmat3 = DeltaABS
        CALL MPI_REDUCE(cmat3, DeltaABS,  NpKx*NpKy*Npz, MPI_COMPLEX, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        cmat3 = DeltaAPUp
        CALL MPI_REDUCE(cmat3, DeltaAPUp, NpKx*NpKy*Npz, MPI_COMPLEX, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        cmat3 = DeltaAPDn
        CALL MPI_REDUCE(cmat3, DeltaAPDn, NpKx*NpKy*Npz, MPI_COMPLEX, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        cmat3 = DeltaBPUp
        CALL MPI_REDUCE(cmat3, DeltaBPUp, NpKx*NpKy*Npz, MPI_COMPLEX, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        cmat3 = DeltaBPDn
        CALL MPI_REDUCE(cmat3, DeltaBPDn, NpKx*NpKy*Npz, MPI_COMPLEX, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        cmat3 = DeltaABPUp
        CALL MPI_REDUCE(cmat3, DeltaABPUp,NpKx*NpKy*Npz, MPI_COMPLEX, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        cmat3 = DeltaABPDn
        CALL MPI_REDUCE(cmat3, DeltaABPDn,NpKx*NpKy*Npz, MPI_COMPLEX, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        DEALLOCATE(cmat3)

        ALLOCATE(cmat5(NpKx, NpKy, Npz, Norb/2, Norb/2))
        cmat5 = DeltaMat
        CALL MPI_REDUCE(cmat5, DeltaMat, NpKx*NpKy*Npz*Norb/2*Norb/2, MPI_COMPLEX, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        DEALLOCATE(cmat5)

        ALLOCATE(rmat3(NpKx, NpKy, 5))
        rmat3 = BStruc
        CALL MPI_REDUCE(rmat3, BStruc, NpKx*NpKy*5, MPI_REAL, MPI_SUM, mpiroot, MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 3)
        DEALLOCATE(rmat3)
    ENDIF
END SUBROUTINE GatherDataMPI


SUBROUTINE FinalizeMPI
    USE NEGF_variables
    USE MPI
    IMPLICIT NONE

    IF (Nprocs > 1) THEN
        CALL MPI_BARRIER(MPI_COMM_WORLD, mpierror); CALL ErrorCheckMPI(mpierror, 4)

        CALL MPI_FINALIZE(mpierror); CALL ErrorCheckMPI(mpierror, 4)

        IF (mpirank /= mpiroot) STOP
    ENDIF
END SUBROUTINE FinalizeMPI


SUBROUTINE ErrorCheckMPI(mpierror, condition)
    IMPLICIT NONE
    INTEGER :: mpierror, condition

    IF (mpierror /= 0) THEN
        IF (condition == 1) THEN
            WRITE(*,*) "MPI Error in Initialize MPI Subroutine.  mpierror = ", mpierror
        ELSEIF (condition == 2) THEN
            WRITE(*,*) "MPI Error in Broadcast MPI Subroutine.  mpierror = ", mpierror
        ELSEIF (condition == 3) THEN
            WRITE(*,*) "MPI Error in GatherData MPI Subroutine.  mpierror = ", mpierror
        ELSEIF (condition == 4) THEN
            WRITE(*,*) "MPI Error in Finalize MPI Subroutine.  mpierror = ", mpierror
        ELSE
            WRITE(*,*) "MPI Error in unknown MPI Subroutine.  mpierror = ", mpierror
        ENDIF
    ENDIF  !If error occurred
END SUBROUTINE ErrorCheckMPI
