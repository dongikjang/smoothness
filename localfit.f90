SUBROUTINE LOCALFIT (LAMBDA, Y, N, D, QV, FITTED, TRA)

    INTEGER :: N
    DOUBLE PRECISION :: LAMBDA(N), Y(N), D(N), QV(N,N)
    DOUBLE PRECISION :: FITTED(N), TRA(N)

    INTEGER :: I, J
    DOUBLE PRECISION :: DTMP(N), TMPA, TMPVAL
    INTEGER :: NTH = 4

    DO I=1,N
        DTMP(:) = 1./(1. + LAMBDA(I)*D(:)) 
        TMPVAL = 0.
        DO J=1, N
            TMPA = SUM(QV(I, :) * DTMP(:) * QV(J, :))
            TMPVAL = TMPVAL + TMPA * Y(J)
            IF(J .EQ. I) TRA(I) = TMPA
        END DO    
        FITTED(I) = TMPVAL
    END DO

END SUBROUTINE