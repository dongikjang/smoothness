SUBROUTINE LOCALRISK1 (LAMBDA, NH, D, QV, N, PILOTF, PILOTS, RISK)
! Estimated or true

    IMPLICIT NONE
    INTEGER :: NH, N
    DOUBLE PRECISION :: LAMBDA(NH), D(N), QV(N,N), PILOTF(N), PILOTS
    DOUBLE PRECISION :: RISK(N, NH)
    INTEGER :: I, J
    DOUBLE PRECISION :: IDL(N), DL(N), TRA2(N), IA(N,N), TQV(N,N)
    INTEGER :: NTH = 4
    

    TQV = TRANSPOSE(QV) 

    DO I = 1, NH
        IDL(:) =  LAMBDA(I) *D(:)/(1. + LAMBDA(I)*D(:)) 
        DL(:) = 1./(1. + LAMBDA(I)*D(:))
        
        DO J = 1, N
            TRA2(J) = SUM((QV(J,:)**2) * (DL(:)**2))
            IA(J,:) = QV(J,:) * IDL(:) 
        END DO
        
        IA = MATMUL(IA, TQV)
        RISK(:,I) = (MATMUL(IA, PILOTF)**2) + (PILOTS**2) * TRA2(:)
    END DO

END SUBROUTINE LOCALRISK1

SUBROUTINE LOCALRISK2 (LAMBDA, NH, D, QV, N, PILOTF, PILOTS, RISK)
! Bias-corrected
! PILOTF is y

    IMPLICIT NONE
    INTEGER :: NH, N
    DOUBLE PRECISION :: LAMBDA(NH), D(N), QV(N,N), PILOTF(N), PILOTS
    DOUBLE PRECISION :: RISK(N, NH)
    INTEGER :: I, J
    DOUBLE PRECISION :: IDL(N), DL(N), TRA2(N), TRIA2(N), IA(N,N), TQV(N,N)
    INTEGER :: NTH = 4
    

    TQV = TRANSPOSE(QV) 

    DO I = 1, NH
        IDL(:) =  LAMBDA(I) *D(:)/(1. + LAMBDA(I)*D(:)) 
        DL(:) = 1./(1. + LAMBDA(I)*D(:))
        
        DO J = 1, N
            TRA2(J) = SUM((QV(J,:)**2) * (DL(:)**2))
            IA(J,:) = QV(J,:) * IDL(:) 
            TRIA2(J) = SUM((QV(J,:)**2) * (IDL(:)**2))
        END DO
        
        IA = MATMUL(IA, TQV)
        RISK(:,I) = (MATMUL(IA, PILOTF)**2) + (PILOTS**2) * (TRA2(:) - TRIA2(:))
    END DO

END SUBROUTINE LOCALRISK2


SUBROUTINE IATRA2 (LAMBDA, NH, D, QV, N, IAALL, TRA2ALL)
! Estimated or true

    IMPLICIT NONE
    INTEGER :: NH, N
    DOUBLE PRECISION :: LAMBDA(NH), D(N), QV(N,N), PILOTF(N), PILOTS
    DOUBLE PRECISION :: RISK(N, NH)
    INTEGER :: I, J
    DOUBLE PRECISION :: IDL(N), DL(N), TRA2(N), IA(N,N), TQV(N,N)
    INTEGER :: NTH = 4
    DOUBLE PRECISION :: IAALL(N*N, NH), TRA2ALL(N, NH)

    TQV = TRANSPOSE(QV) 

    DO I = 1, NH
        IDL(:) =  LAMBDA(I) *D(:)/(1. + LAMBDA(I)*D(:)) 
        DL(:) = 1./(1. + LAMBDA(I)*D(:))
        
        DO J = 1, N
            TRA2(J) = SUM((QV(J,:)**2) * (DL(:)**2))
            IA(J,:) = QV(J,:) * IDL(:) 
        END DO
        
        IA = MATMUL(IA, TQV)
        DO J = 1, N
            IAALL( (N*(J-1)+1):(N*J), I) = IA(:, J )
        END DO
        !IAALL(:, I) = RESHAPE(IA, N*N )
        TRA2ALL(:, I) = TRA2(:)
        
        ! RISK(:,I) = (MATMUL(IA, PILOTF)**2) + (PILOTS**2) * TRA2(:)
    END DO

END SUBROUTINE IATRA2


SUBROUTINE LOCALRISKALL (N, NH, PILOTF, PILOTS, IAALL, TRA2ALL, RISK)

    IMPLICIT NONE
    INTEGER :: NH, N
    INTEGER :: I, J
    
    DOUBLE PRECISION :: RISK(N, NH), IAALL(N*N, NH), TRA2ALL(N, NH)
    DOUBLE PRECISION :: IA(N,N), TRA2(N), PILOTF(N), PILOTS

    DO I = 1, NH
        !DO J = 1, N
        !    IA(:, J ) = IAALL( (N*(J-1)+1):(N*J), I) 
        !END DO
        IA = RESHAPE(IAALL(:, I), (/ N, N /))
        TRA2 = TRA2ALL(:, I)
        RISK(:,I) = (MATMUL(IA, PILOTF)**2) + (PILOTS**2) * TRA2(:)
    END DO
    
END SUBROUTINE LOCALRISKALL

SUBROUTINE LOCALRISKALL2 (N, NH, PILOTF, PILOTS, IAALL, TRA2ALL, RISK, BIAS, VARI)

    IMPLICIT NONE
    INTEGER :: NH, N
    INTEGER :: I, J
    
    DOUBLE PRECISION :: RISK(N, NH), IAALL(N*N, NH), TRA2ALL(N, NH), BIAS(N, NH), VARI(N, NH)
    DOUBLE PRECISION :: IA(N,N), TRA2(N), PILOTF(N), PILOTS

    DO I = 1, NH
        !DO J = 1, N
        !    IA(:, J ) = IAALL( (N*(J-1)+1):(N*J), I) 
        !END DO
        IA = RESHAPE(IAALL(:, I), (/ N, N /))
        TRA2 = TRA2ALL(:, I)
        BIAS(:,I) = MATMUL(IA, PILOTF)
        VARI(:,I) = (PILOTS**2) * TRA2(:)
        RISK(:,I) = BIAS(:,I)**2 + VARI(:,I)
    END DO
    
END SUBROUTINE LOCALRISKALL2


SUBROUTINE MOVINGRISK(LRISK, N, NH, X, NEWX, NN, RVAL, MOVRISK, MLOCRISK)

    IMPLICIT NONE
    INTEGER :: N, NN, NH
    DOUBLE PRECISION :: LRISK(N, NH), X(N), NEWX(NN), RVAL, MOVRISK(NN)
    INTEGER :: MLOCRISK(NN)   
    
    INTEGER :: I, J, IND(N), TMPIND(N, NH), IND2
    DOUBLE PRECISION :: TMPVAL(N), TMPRISK(NH)
    INTEGER :: NTH = 4

    DO I = 1, NN
        IND(:) = 0
        TMPVAL(:) = ABS(X(:) - NEWX(I)) - RVAL
        DO J = 1, N
            IF(TMPVAL(J) .LT. 0 ) IND(J) = 1

        END DO

        TMPIND = SPREAD(IND, DIM=2, NCOPIES=NH)
        TMPRISK = SUM(LRISK, DIM=1, MASK= TMPIND .GT. 0)
        
        TMPRISK(:) = TMPRISK(:)/REAL(SUM(IND))
        IND2 = MINLOC(TMPRISK, DIM=1)
        MLOCRISK(I) = IND2
        MOVRISK(I) = TMPRISK(IND2)

    END DO

END SUBROUTINE MOVINGRISK

SUBROUTINE MOVINGRISK3(LRISK, N, NH, X, NEWX, NN, RVAL, MOVRISK, MLOCRISK)

    IMPLICIT NONE
    INTEGER :: N, NN, NH
    DOUBLE PRECISION :: LRISK(N, NH), X(N), NEWX(NN), RVAL, MOVRISK(NN)
    INTEGER :: MLOCRISK(NN)   
    
    INTEGER :: I, J, IND(N), TMPIND(N, NH), IND2
    DOUBLE PRECISION :: TMPRISK(NH)
    INTEGER :: NTH = 4

	DOUBLE PRECISION :: TMPVAL, TMPSUM(NH), FLIPBD, MAXX, MINX
	DOUBLE PRECISION :: NLLEN, NRLEN, LLEN, RLEN
	INTEGER :: ITER
    
	
	MAXX = MAXVAL(X)
	MINX = MINVAL(X)
	FLIPBD = MIN(RVAL, (MAXX - MINX)/4)
	DO I = 1, NN  ! newx index
        TMPSUM(:) = 0.0
		ITER = 0
		! SUM LRISK( , :) * I(|X(J) - NEWX(I)| < RVAL)
		DO J = 1, N ! x index
			TMPVAL = ABS(X(J) - NEWX(I)) - RVAL
			IF (TMPVAL .LT. 0.0) THEN
				TMPSUM = TMPSUM + LRISK(J,:)
				ITER = ITER + 1
			ENDIF
			! Left boundary correct
			NLLEN = NEWX(I) - MINX
			NRLEN = MAXX - NEWX(I)
			LLEN = X(J) - MINX
			RLEN = MAXX - X(J)
			
			IF ((NLLEN .LT. RVAL) .AND. (LLEN .LT. (FLIPBD - NLLEN))) THEN
				TMPSUM = TMPSUM + LRISK(J,:)
				ITER = ITER + 1
			ENDIF
			
			IF ((NRLEN .LT. RVAL) .AND. (RLEN .LT. (FLIPBD - NRLEN))) THEN
				TMPSUM = TMPSUM + LRISK(J,:)
				ITER = ITER + 1
			ENDIF
			
		END DO

		TMPRISK = TMPSUM(:)/ITER
		IND2 = MINLOC(TMPRISK, DIM=1)
		MLOCRISK(I) = IND2
		MOVRISK(I) = TMPRISK(IND2)

    END DO

END SUBROUTINE MOVINGRISK3


SUBROUTINE MOVINGRISK2(LRISK, N, NH, X, NEWX, NN, RVAL, MOVRISK, MLOCRISK)

    IMPLICIT NONE
    INTEGER :: N, NN, NH
    DOUBLE PRECISION :: LRISK(N, NH), X(N, 2), NEWX(NN, 2), RVAL, MOVRISK(NN)
    INTEGER :: MLOCRISK(NN)   
    
    INTEGER :: I, J, IND(N), TMPIND(N, NH), IND2
    DOUBLE PRECISION :: TMPVAL(N), TMPRISK(NH)
    INTEGER :: NTH = 4

    DO I = 1, NN

        IND(:) = 0
        TMPVAL(:) = (X(:,1) - NEWX(I,1))**2 + (X(:,2) - NEWX(I,2))**2 - RVAL**2
        DO J = 1, N
            IF(TMPVAL(J) .LT. 0 ) IND(J) = 1
        END DO

        TMPIND = SPREAD(IND, DIM=2, NCOPIES=NH)
        TMPRISK = SUM(LRISK, DIM=1, MASK= TMPIND .GT. 0)
        TMPRISK(:) = TMPRISK(:)/REAL(SUM(IND))
        IND2 = MINLOC(TMPRISK, DIM=1)
        MLOCRISK(I) = IND2
        MOVRISK(I) = TMPRISK(IND2)
    END DO

END SUBROUTINE MOVINGRISK2