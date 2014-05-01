PROGRAM NEGF_main
    USE NEGF_variables
    IMPLICIT NONE

    CALL Read_In_Allocate

    CALL BroadcastDataMPI !Broadcast mpiroot information to non-root processes

    DO xx = 1, NpKx
        DO yy = 1, NpKy
            !Split up k-space calculations over each process
            IF ((mpisize == 1) .OR. (mod(yy+NpKy*(xx-1), mpisize) == mpirank)) THEN
                Kx = KxAxis(xx); Ky = KyAxis(yy);

                !Set initial conditions for Order Parameter in Superconductor
                DeltaAS(xx,yy,1:NpzSC) = DeltaSC; DeltaBS(xx,yy,1:NpzSC) = DeltaSC;
                DO zz = 1, Npz
                    DeltaMat(xx,yy,zz,:,:) = DeltaSC * (gDAS + gDBS)
                ENDDO

                ite = 0; AmtChange = 1D0
                DO WHILE ((ite <= iteTot) .AND. (abs(AmtChange)>Accuracy)) !Wait for Self-consistent solution to order parameter
                    ite = ite+1

                    CALL Build_Hamiltonian !Build Htot(Kx, Ky, z)

                    CALL Calculate_OP
                ENDDO

                CALL CalculateFinal

                IF (mpirank == mpiroot) THEN
                    WRITE(*,*) "ite xx yy, AmtChange, EGndChange = ", ite, xx, yy, AmtChange, EGnd(xx,yy)-EGnd0(xx,yy)
                ENDIF
            ENDIF
        ENDDO
    ENDDO

    !Gather Data from each process, store in mpiroot process
    CALL GatherDataMPI

    CALL FinalizeMPI

    IF (mpirank==mpiroot) THEN !only one process needs to do this
        CALL DumpData !dump data to file and stop
    ENDIF
END PROGRAM NEGF_main
