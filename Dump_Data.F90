SUBROUTINE dumpdata
    USE NEGF_variables
    IMPLICIT NONE
    CHARACTER :: Npxstring*10, Npystring*10

    !Open files to store scalar data
    OPEN(unit=13, file=trim(outputfiledir)//"/scalars.txt")
    WRITE(13, *) NpKx
    WRITE(13, *) NpKy
    WRITE(13, *) Npz
    WRITE(13, *) NpzSC
    WRITE(13, *) NpE
    WRITE(13, *) Norb
    WRITE(13, 43) muTI
    WRITE(13, 43) muSC
    WRITE(13, 43) UIntTI
    WRITE(13, 43) UIntSC
    WRITE(13, *) MatType
    WRITE(13, *) MM
    WRITE(13, *) a0
    WRITE(13, *) DeltaSC
    WRITE(13, *) Mx
    WRITE(13, *) My
    WRITE(13, *) Mz
    CLOSE(13)

    OPEN(unit=13, file=trim(outputfiledir)//"/Eigenstates.txt")
    IF (fullDeltaMat==1) THEN
        DO xx = 1, NpKx
            DO yy = 1, NpKy
                DO zz = 1, Npz
                    DO ii = 1, Norb/2
                        DO jj = 1, Norb/2
                            WRITE(13, 43) real(DeltaMat(xx,yy,zz,ii,jj))
                            WRITE(13, 43) imag(DeltaMat(xx,yy,zz,ii,jj))
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ELSE
        DO xx = 1, NpKx
            DO yy = 1, NpKy
                DO zz = 1, Npz
                    WRITE(13, 43) real(DeltaAS(xx,yy,zz))
                    WRITE(13, 43) imag(DeltaAS(xx,yy,zz))

                    WRITE(13, 43) real(DeltaBS(xx,yy,zz))
                    WRITE(13, 43) imag(DeltaBS(xx,yy,zz))

                    WRITE(13, 43) real(DeltaAPUp(xx,yy,zz))
                    WRITE(13, 43) imag(DeltaAPUp(xx,yy,zz))

                    WRITE(13, 43) real(DeltaAPDn(xx,yy,zz))
                    WRITE(13, 43) imag(DeltaAPDn(xx,yy,zz))

                    WRITE(13, 43) real(DeltaBPUp(xx,yy,zz))
                    WRITE(13, 43) imag(DeltaBPUp(xx,yy,zz))

                    WRITE(13, 43) real(DeltaBPDn(xx,yy,zz))
                    WRITE(13, 43) imag(DeltaBPDn(xx,yy,zz))

                    WRITE(13, 43) real(DeltaABPUp(xx,yy,zz))
                    WRITE(13, 43) imag(DeltaABPUp(xx,yy,zz))

                    WRITE(13, 43) real(DeltaABPDn(xx,yy,zz))
                    WRITE(13, 43) imag(DeltaABPDn(xx,yy,zz))
                ENDDO
            ENDDO
        ENDDO
    ENDIF
    CLOSE(13)

    OPEN(unit=13, file=trim(outputfiledir)//"/KAxes.txt")
    DO xx = 1, NpKx
        WRITE(13, *) KxAxis(xx)
    ENDDO
    DO yy = 1, NpKy
        WRITE(13, *) KyAxis(yy)
    ENDDO
    DO ee = 1, NpE
        WRITE(13, *) Eaxis(ee)
    ENDDO
    DO xx = 1, NpKx
        DO yy = 1, NpKy
            WRITE(13, *) EGnd(xx,yy)
            WRITE(13, *) EGnd0(xx,yy)
        ENDDO
    ENDDO
    DO xx = 1, NpKx
        DO yy = 1, NpKy
            DO zz = 1, Npz
                DO ee = 1, NpE
                    WRITE(13, 43) LDOS(xx,yy,zz,ee)
                    WRITE(13, 43) LDOSMx(xx,yy,zz,ee)
                    WRITE(13, 43) LDOSMy(xx,yy,zz,ee)
                    WRITE(13, 43) LDOSMz(xx,yy,zz,ee)
                    IF (MatType == 1) THEN
                        WRITE(13, 43) LDOSPx(xx,yy,zz,ee)
                        WRITE(13, 43) LDOSPy(xx,yy,zz,ee)
                        WRITE(13, 43) LDOSPz(xx,yy,zz,ee)
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO
    CLOSE(13)

    OPEN(unit=13, file=trim(outputfiledir)//"/BStruc.txt")
    DO xx = 1, NpKx
        DO yy = 1, NpKy
            DO ii = 1, 5
                WRITE(13, *) BStruc(xx,yy,ii)
            ENDDO
        ENDDO
    ENDDO
    CLOSE(13)

    WRITE(*,*)"COMPLETE:  Data dumped to directory ", outputfiledir
    STOP "Run Complete!"

43 format(ES11.4)
END SUBROUTINE dumpdata
