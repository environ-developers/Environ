CC     ----------------------------------------------------------------------
CC     This file contains the LBFGS algorithm and supporting routines
CC
CC     ****************
CC     LBFGS SUBROUTINE
CC     ****************
CC
C      SUBROUTINE LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
CC
C      INTEGER N,M,IPRINT(2),IFLAG
C      DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
C      DOUBLE PRECISION F,EPS,XTOL
C      LOGICAL DIAGCO
CC
CC        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
CC                          JORGE NOCEDAL
CC                        *** July 1990 ***
CC
CC 
CC     This subroutine solves the unconstrained minimization problem
CC 
CC                      min F(x),    x= (x1,x2,...,xN),
CC
CC      using the limited memory BFGS method. The routine is especially
CC      effective on problems involving a large number of variables. In
CC      a typical iteration of this method an approximation Hk to the
CC      inverse of the Hessian is obtained by applying M BFGS updates to
CC      a diagonal matrix Hk0, using information from the previous M steps.
CC      The user specifies the number M, which determines the amount of
CC      storage required by the routine. The user may also provide the
CC      diagonal matrices Hk0 if not satisfied with the default choice.
CC      The algorithm is described in "On the limited memory BFGS method
CC      for large scale optimization", by D. Liu and J. Nocedal,
CC      Mathematical Programming B 45 (1989) 503-528.
CC 
CC      The user is required to calculate the function value F and its
CC      gradient G. In order to allow the user complete control over
CC      these computations, reverse  communication is used. The routine
CC      must be called repeatedly under the control of the parameter
CC      IFLAG. 
CC
CC      The steplength is determined at each iteration by means of the
CC      line search routine MCVSRCH, which is a slight modification of
CC      the routine CSRCH written by More' and Thuente.
CC 
CC      The calling statement is 
CC 
CC          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
CC 
CC      where
CC 
CC     N       is an INTEGER variable that must be set by the user to the
CC             number of variables. It is not altered by the routine.
CC             Restriction: N>0.
CC 
CC     M       is an INTEGER variable that must be set by the user to
CC             the number of corrections used in the BFGS update. It
CC             is not altered by the routine. Values of M less than 3 are
CC             not recommended; large values of M will result in excessive
CC             computing time. 3<= M <=7 is recommended. Restriction: M>0.
CC 
CC     X       is a DOUBLE PRECISION array of length N. On initial entry
CC             it must be set by the user to the values of the initial
CC             estimate of the solution vector. On exit with IFLAG=0, it
CC             contains the values of the variables at the best point
CC             found (usually a solution).
CC 
CC     F       is a DOUBLE PRECISION variable. Before initial entry and on
CC             a re-entry with IFLAG=1, it must be set by the user to
CC             contain the value of the function F at the point X.
CC 
CC     G       is a DOUBLE PRECISION array of length N. Before initial
CC             entry and on a re-entry with IFLAG=1, it must be set by
CC             the user to contain the components of the gradient G at
CC             the point X.
CC 
CC     DIAGCO  is a LOGICAL variable that must be set to .TRUE. if the
CC             user  wishes to provide the diagonal matrix Hk0 at each
CC             iteration. Otherwise it should be set to .FALSE., in which
CC             case  LBFGS will use a default value described below. If
CC             DIAGCO is set to .TRUE. the routine will return at each
CC             iteration of the algorithm with IFLAG=2, and the diagonal
CC              matrix Hk0  must be provided in the array DIAG.
CC 
CC 
CC     DIAG    is a DOUBLE PRECISION array of length N. If DIAGCO=.TRUE.,
CC             then on initial entry or on re-entry with IFLAG=2, DIAG
CC             it must be set by the user to contain the values of the 
CC             diagonal matrix Hk0.  Restriction: all elements of DIAG
CC             must be positive.
CC 
CC     IPRINT  is an INTEGER array of length two which must be set by the
CC             user.
CC 
CC             IPRINT(1) specifies the frequency of the output:
CC                IPRINT(1) < 0 : no output is generated,
CC                IPRINT(1) = 0 : output only at first and last iteration,
CC                IPRINT(1) > 0 : output every IPRINT(1) iterations.
CC 
CC             IPRINT(2) specifies the type of output generated:
CC                IPRINT(2) = 0 : iteration count, number of function 
CC                                evaluations, function value, norm of the
CC                                gradient, and steplength,
CC                IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of
CC                                variables and  gradient vector at the
CC                                initial point,
CC                IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of
CC                                variables,
CC                IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.
CC 
CC 
CC     EPS     is a positive DOUBLE PRECISION variable that must be set by
CC             the user, and determines the accuracy with which the solution
CC             is to be found. The subroutine terminates when
CC
CC                         ||G|| < EPS max(1,||X||),
CC
CC             where ||.|| denotes the Euclidean norm.
CC 
CC     XTOL    is a  positive DOUBLE PRECISION variable that must be set by
CC             the user to an estimate of the machine precision (e.g.
CC             10**(-16) on a SUN station 3/60). The line search routine will
CC             terminate if the relative width of the interval of uncertainty
CC             is less than XTOL.
CC 
CC     W       is a DOUBLE PRECISION array of length N(2M+1)+2M used as
CC             workspace for LBFGS. This array must not be altered by the
CC             user.
CC 
CC     IFLAG   is an INTEGER variable that must be set to 0 on initial entry
CC             to the subroutine. A return with IFLAG<0 indicates an error,
CC             and IFLAG=0 indicates that the routine has terminated without
CC             detecting errors. On a return with IFLAG=1, the user must
CC             evaluate the function F and gradient G. On a return with
CC             IFLAG=2, the user must provide the diagonal matrix Hk0.
CC 
CC             The following negative values of IFLAG, detecting an error,
CC             are possible:
CC 
CC              IFLAG=-1  The line search routine MCSRCH failed. The
CC                        parameter INFO provides more detailed information
CC                        (see also the documentation of MCSRCH):
CC
CC                       INFO = 0  IMPROPER INPUT PARAMETERS.
CC
CC                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
CC                                 UNCERTAINTY IS AT MOST XTOL.
CC
CC                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE
CC                                 REQUIRED AT THE PRESENT ITERATION.
CC
CC                       INFO = 4  THE STEP IS TOO SMALL.
CC
CC                       INFO = 5  THE STEP IS TOO LARGE.
CC
CC                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. 
CC                                 THERE MAY NOT BE A STEP WHICH SATISFIES
CC                                 THE SUFFICIENT DECREASE AND CURVATURE
CC                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL.
CC
CC 
CC              IFLAG=-2  The i-th diagonal element of the diagonal inverse
CC                        Hessian approximation, given in DIAG, is not
CC                        positive.
CC           
CC              IFLAG=-3  Improper input parameters for LBFGS (N or M are
CC                        not positive).
CC 
CC
CC
CC    ON THE DRIVER:
CC
CC    The program that calls LBFGS must contain the declaration:
CC
CC                       EXTERNAL LB2
CC
CC    LB2 is a BLOCK DATA that defines the default values of several
CC    parameters described in the COMMON section. 
CC
CC 
CC 
CC    COMMON:
CC 
CC     The subroutine contains one common area, which the user may wish to
CC    reference:
CC 
C         COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
CC 
CC    MP  is an INTEGER variable with default value 6. It is used as the
CC        unit number for the printing of the monitoring information
CC        controlled by IPRINT.
CC 
CC    LP  is an INTEGER variable with default value 6. It is used as the
CC        unit number for the printing of error messages. This printing
CC        may be suppressed by setting LP to a non-positive value.
CC 
CC    GTOL is a DOUBLE PRECISION variable with default value 0.9, which
CC        controls the accuracy of the line search routine MCSRCH. If the
CC        function and gradient evaluations are inexpensive with respect
CC        to the cost of the iteration (which is sometimes the case when
CC        solving very large problems) it may be advantageous to set GTOL
CC        to a small value. A typical small value is 0.1.  Restriction:
CC        GTOL should be greater than 1.D-04.
CC 
CC    STPMIN and STPMAX are non-negative DOUBLE PRECISION variables which
CC        specify lower and uper bounds for the step in the line search.
CC        Their default values are 1.D-20 and 1.D+20, respectively. These
CC        values need not be modified unless the exponents are too large
CC        for the machine being used, or unless the problem is extremely
CC        badly scaled (in which case the exponents should be increased).
CC 
CC
CC  MACHINE DEPENDENCIES
CC
CC        The only variables that are machine-dependent are XTOL,
CC        STPMIN and STPMAX.
CC 
CC
CC  GENERAL INFORMATION
CC 
CC    Other routines called directly:  DAXPY, DDOT, LB1, MCSRCH
CC 
CC    Input/Output  :  No input; diagnostic messages on unit MP and
CC                     error messages on unit LP.
CC 
CC 
CC     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CC
C      DOUBLE PRECISION GTOL,ONE,ZERO,GNORM,DDOT,STP1,FTOL,STPMIN,
C     .                 STPMAX,STP,YS,YY,SQ,YR,BETA,XNORM
C      INTEGER MP,LP,ITER,NFUN,POINT,ISPT,IYPT,MAXFEV,INFO,
C     .        BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN
C      LOGICAL FINISH
CC
C      SAVE
C      DATA ONE,ZERO/1.0D+0,0.0D+0/
CC
CC     INITIALIZE
CC     ----------
CC
C      IF(IFLAG.EQ.0) GO TO 10
C      GO TO (172,100) IFLAG
C  10  ITER= 0
C      IF(N.LE.0.OR.M.LE.0) GO TO 196
C      IF(GTOL.LE.1.D-04) THEN
C        IF(LP.GT.0) WRITE(LP,245)
C        GTOL=9.D-01
C      ENDIF
C      NFUN= 1
C      POINT= 0
C      FINISH= .FALSE.
C      IF(DIAGCO) THEN
C         DO 30 I=1,N
C 30      IF (DIAG(I).LE.ZERO) GO TO 195
C      ELSE
C         DO 40 I=1,N
C 40      DIAG(I)= 1.0D0
C      ENDIF
CC
CC     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
CC     ---------------------------------------
CC     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
CC         OTHER TEMPORARY INFORMATION.
CC     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
CC     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
CC         IN THE FORMULA THAT COMPUTES H*G.
CC     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
CC         STEPS.
CC     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
CC         GRADIENT DIFFERENCES.
CC
CC     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
CC     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
CC
C      ISPT= N+2*M
C      IYPT= ISPT+N*M     
C      DO 50 I=1,N
C 50   W(ISPT+I)= -G(I)*DIAG(I)
C      GNORM= DSQRT(DDOT(N,G,1,G,1))
C      STP1= ONE/GNORM
CC
CC     PARAMETERS FOR LINE SEARCH ROUTINE
CC     
C      FTOL= 1.0D-4
C      MAXFEV= 20
CC
C      IF(IPRINT(1).GE.0) CALL LB1(IPRINT,ITER,NFUN,
C     *                     GNORM,N,M,X,F,G,STP,FINISH)
CC
CC    --------------------
CC     MAIN ITERATION LOOP
CC    --------------------
CC
C 80   ITER= ITER+1
C      INFO=0
C      BOUND=ITER-1
C      IF(ITER.EQ.1) GO TO 165
C      IF (ITER .GT. M)BOUND=M
CC
C         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
C      IF(.NOT.DIAGCO) THEN
C         YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
C         DO 90 I=1,N
C   90    DIAG(I)= YS/YY
C      ELSE
C         IFLAG=2
C         RETURN
C      ENDIF
C 100  CONTINUE
C      IF(DIAGCO) THEN
C        DO 110 I=1,N
C 110    IF (DIAG(I).LE.ZERO) GO TO 195
C      ENDIF
CC
CC     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
CC     "Updating quasi-Newton matrices with limited storage",
CC     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
CC     ---------------------------------------------------------
CC
C      CP= POINT
C      IF (POINT.EQ.0) CP=M
C      W(N+CP)= ONE/YS
C      DO 112 I=1,N
C 112  W(I)= -G(I)
C      CP= POINT
C      DO 125 I= 1,BOUND
C         CP=CP-1
C         IF (CP.EQ. -1)CP=M-1
C         SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
C         INMC=N+M+CP+1
C         IYCN=IYPT+CP*N
C         W(INMC)= W(N+CP+1)*SQ
C         CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
C 125  CONTINUE
CC
C      DO 130 I=1,N
C 130  W(I)=DIAG(I)*W(I)
CC
C      DO 145 I=1,BOUND
C         YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
C         BETA= W(N+CP+1)*YR
C         INMC=N+M+CP+1
C         BETA= W(INMC)-BETA
C         ISCN=ISPT+CP*N
C         CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
C         CP=CP+1
C         IF (CP.EQ.M)CP=0
C 145  CONTINUE
CC
CC     STORE THE NEW SEARCH DIRECTION
CC     ------------------------------
CC
C       DO 160 I=1,N
C 160   W(ISPT+POINT*N+I)= W(I)
CC
CC     OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION 
CC     BY USING THE LINE SEARCH ROUTINE MCSRCH
CC     ----------------------------------------------------
C 165  NFEV=0
C      STP=ONE
C      IF (ITER.EQ.1) STP=STP1
C      DO 170 I=1,N
C 170  W(I)=G(I)
C 172  CONTINUE
C      CALL MCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP,FTOL,
C     *            XTOL,MAXFEV,INFO,NFEV,DIAG)
C      IF (INFO .EQ. -1) THEN
C        IFLAG=1
C        RETURN
C      ENDIF
C      IF (INFO .NE. 1) GO TO 190
C      NFUN= NFUN + NFEV
CC
CC     COMPUTE THE NEW STEP AND GRADIENT CHANGE 
CC     -----------------------------------------
CC
C      NPT=POINT*N
C      DO 175 I=1,N
C      W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
C 175  W(IYPT+NPT+I)= G(I)-W(I)
C      POINT=POINT+1
C      IF (POINT.EQ.M)POINT=0
CC
CC     TERMINATION TEST
CC     ----------------
CC
C      GNORM= DSQRT(DDOT(N,G,1,G,1))
C      XNORM= DSQRT(DDOT(N,X,1,X,1))
C      XNORM= DMAX1(1.0D0,XNORM)
C      IF (GNORM/XNORM .LE. EPS) FINISH=.TRUE.
CC
C      IF(IPRINT(1).GE.0) CALL LB1(IPRINT,ITER,NFUN,
C     *               GNORM,N,M,X,F,G,STP,FINISH)
C      IF (FINISH) THEN
C         IFLAG=0
C         RETURN
C      ENDIF
C      GO TO 80
CC
CC     ------------------------------------------------------------
CC     END OF MAIN ITERATION LOOP. ERROR EXITS.
CC     ------------------------------------------------------------
CC
C 190  IFLAG=-1
C      IF(LP.GT.0) WRITE(LP,200) INFO
C      RETURN
C 195  IFLAG=-2
C      IF(LP.GT.0) WRITE(LP,235) I
C      RETURN
C 196  IFLAG= -3
C      IF(LP.GT.0) WRITE(LP,240)
CC
CC     FORMATS
CC     -------
CC
C 200  FORMAT(/' IFLAG= -1 ',/' LINE SEARCH FAILED. SEE'
C     .          ' DOCUMENTATION OF ROUTINE MCSRCH',/' ERROR RETURN'
C     .          ' OF LINE SEARCH: INFO= ',I2,/
C     .          ' POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT',/,
C     .          ' OR INCORRECT TOLERANCES')
C 235  FORMAT(/' IFLAG= -2',/' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,
C     .       ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
C 240  FORMAT(/' IFLAG= -3',/' IMPROPER INPUT PARAMETERS (N OR M',
C     .       ' ARE NOT POSITIVE)')
C 245  FORMAT(/'  GTOL IS LESS THAN OR EQUAL TO 1.D-04',
C     .       / ' IT HAS BEEN RESET TO 9.D-01')
C      RETURN
C      END
CC
CC     LAST LINE OF SUBROUTINE LBFGS
CC
CC
C      SUBROUTINE LB1(IPRINT,ITER,NFUN,
C     *                     GNORM,N,M,X,F,G,STP,FINISH)
CC
CC     -------------------------------------------------------------
CC     THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND
CC     AMOUNT OF OUTPUT ARE CONTROLLED BY IPRINT.
CC     -------------------------------------------------------------
CC
C      INTEGER IPRINT(2),ITER,NFUN,LP,MP,N,M
C      DOUBLE PRECISION X(N),G(N),F,GNORM,STP,GTOL,STPMIN,STPMAX
C      LOGICAL FINISH
C      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
CC
C      IF (ITER.EQ.0)THEN
C           WRITE(MP,10)
C           WRITE(MP,20) N,M
C           WRITE(MP,30)F,GNORM
C                 IF (IPRINT(2).GE.1)THEN
C                     WRITE(MP,40)
C                     WRITE(MP,50) (X(I),I=1,N)
C                     WRITE(MP,60)
C                     WRITE(MP,50) (G(I),I=1,N)
C                  ENDIF
C           WRITE(MP,10)
C           WRITE(MP,70)
C      ELSE
C          IF ((IPRINT(1).EQ.0).AND.(ITER.NE.1.AND..NOT.FINISH))RETURN
C              IF (IPRINT(1).NE.0)THEN
C                   IF(MOD(ITER-1,IPRINT(1)).EQ.0.OR.FINISH)THEN
C                         IF(IPRINT(2).GT.1.AND.ITER.GT.1) WRITE(MP,70)
C                         WRITE(MP,80)ITER,NFUN,F,GNORM,STP
C                   ELSE
C                         RETURN
C                   ENDIF
C              ELSE
C                   IF( IPRINT(2).GT.1.AND.FINISH) WRITE(MP,70)
C                   WRITE(MP,80)ITER,NFUN,F,GNORM,STP
C              ENDIF
C              IF (IPRINT(2).EQ.2.OR.IPRINT(2).EQ.3)THEN
C                    IF (FINISH)THEN
C                        WRITE(MP,90)
C                    ELSE
C                        WRITE(MP,40)
C                    ENDIF
C                      WRITE(MP,50)(X(I),I=1,N)
C                  IF (IPRINT(2).EQ.3)THEN
C                      WRITE(MP,60)
C                      WRITE(MP,50)(G(I),I=1,N)
C                  ENDIF
C              ENDIF
C            IF (FINISH) WRITE(MP,100)
C      ENDIF
CC
C 10   FORMAT('*************************************************')
C 20   FORMAT('  N=',I5,'   NUMBER OF CORRECTIONS=',I2,
C     .       /,  '       INITIAL VALUES')
C 30   FORMAT(' F= ',1PD10.3,'   GNORM= ',1PD10.3)
C 40   FORMAT(' VECTOR X= ')
C 50   FORMAT(6(2X,1PD10.3))
C 60   FORMAT(' GRADIENT VECTOR G= ')
C 70   FORMAT(/'   I   NFN',4X,'FUNC',8X,'GNORM',7X,'STEPLENGTH'/)
C 80   FORMAT(2(I4,1X),3X,3(1PD10.3,2X))
C 90   FORMAT(' FINAL POINT X= ')
C 100  FORMAT(/' THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.',
C     .       /' IFLAG = 0')
CC
C      RETURN
C      END
CC     ******
CC
CC
CC   ----------------------------------------------------------
CC     DATA 
CC   ----------------------------------------------------------
CC
C      BLOCK DATA LB2
C      INTEGER LP,MP
C      DOUBLE PRECISION GTOL,STPMIN,STPMAX
C      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
C      DATA MP,LP,GTOL,STPMIN,STPMAX/11,11,9.0D-01,1.0D-20,1.0D+20/
C      END
CC
CC
CC   ----------------------------------------------------------
CC    ------------------------------------------------------------------
CC
CC     **************************
CC     LINE SEARCH ROUTINE MCSRCH
CC     **************************
CC
C      SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL,MAXFEV,INFO,NFEV,WA)
C      INTEGER N,MAXFEV,INFO,NFEV
C      DOUBLE PRECISION F,STP,FTOL,GTOL,XTOL,STPMIN,STPMAX
C      DOUBLE PRECISION X(N),G(N),S(N),WA(N)
C      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
C      SAVE
CC
CC                     SUBROUTINE MCSRCH
CC                
CC     A slight modification of the subroutine CSRCH of More' and Thuente.
CC     The changes are to allow reverse communication, and do not affect
CC     the performance of the routine. 
CC
CC     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
CC     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
CC
CC     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
CC     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
CC     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
CC     MINIMIZER OF THE MODIFIED FUNCTION
CC
CC          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
CC
CC     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
CC     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
CC     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
CC     CONTAINS A MINIMIZER OF F(X+STP*S).
CC
CC     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
CC     THE SUFFICIENT DECREASE CONDITION
CC
CC           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),
CC
CC     AND THE CURVATURE CONDITION
CC
CC           ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).
CC
CC     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
CC     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
CC     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
CC     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
CC     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
CC     SATISFIES THE SUFFICIENT DECREASE CONDITION.
CC
CC     THE SUBROUTINE STATEMENT IS
CC
CC        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
CC     WHERE
CC
CC       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
CC         OF VARIABLES.
CC
CC       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
CC         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
CC         X + STP*S.
CC
CC       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
CC         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
CC
CC       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
CC         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
CC         OF F AT X + STP*S.
CC
CC       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
CC         SEARCH DIRECTION.
CC
CC       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
CC         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
CC         STP CONTAINS THE FINAL ESTIMATE.
CC
CC       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
CC         communication implementation GTOL is defined in a COMMON
CC         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
CC         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
CC         SATISFIED.
CC
CC       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
CC         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
CC         IS AT MOST XTOL.
CC
CC       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
CC         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
CC         communication implementatin they are defined in a COMMON
CC         statement).
CC
CC       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
CC         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
CC         MAXFEV BY THE END OF AN ITERATION.
CC
CC       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
CC
CC         INFO = 0  IMPROPER INPUT PARAMETERS.
CC
CC         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
CC
CC         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
CC                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
CC
CC         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
CC                   IS AT MOST XTOL.
CC
CC         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
CC
CC         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
CC
CC         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
CC
CC         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
CC                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
CC                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
CC                   TOLERANCES MAY BE TOO SMALL.
CC
CC       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
CC         CALLS TO FCN.
CC
CC       WA IS A WORK ARRAY OF LENGTH N.
CC
CC     SUBPROGRAMS CALLED
CC
CC       MCSTEP
CC
CC       FORTRAN-SUPPLIED...ABS,MAX,MIN
CC
CC     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
CC     JORGE J. MORE', DAVID J. THUENTE
CC
CC     **********
C      INTEGER INFOC,J
C      LOGICAL BRACKT,STAGE1
C      DOUBLE PRECISION DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,
C     *       FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY,
C     *       STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO
C      DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/
C      IF(INFO.EQ.-1) GO TO 45
C      INFOC = 1
CC
CC     CHECK THE INPUT PARAMETERS FOR ERRORS.
CC
C      IF (N .LE. 0 .OR. STP .LE. ZERO .OR. FTOL .LT. ZERO .OR.
C     *    GTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. STPMIN .LT. ZERO
C     *    .OR. STPMAX .LT. STPMIN .OR. MAXFEV .LE. 0) RETURN
CC
CC     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
CC     AND CHECK THAT S IS A DESCENT DIRECTION.
CC
C      DGINIT = ZERO
C      DO 10 J = 1, N
C         DGINIT = DGINIT + G(J)*S(J)
C   10    CONTINUE
C      IF (DGINIT .GE. ZERO) then
C         write(LP,15)
C   15    FORMAT(/'  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION')
C         RETURN
C         ENDIF
CC
CC     INITIALIZE LOCAL VARIABLES.
CC
C      BRACKT = .FALSE.
C      STAGE1 = .TRUE.
C      NFEV = 0
C      FINIT = F
C      DGTEST = FTOL*DGINIT
C      WIDTH = STPMAX - STPMIN
C      WIDTH1 = WIDTH/P5
C      DO 20 J = 1, N
C         WA(J) = X(J)
C   20    CONTINUE
CC
CC     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
CC     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
CC     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
CC     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
CC     THE INTERVAL OF UNCERTAINTY.
CC     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
CC     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
CC
C      STX = ZERO
C      FX = FINIT
C      DGX = DGINIT
C      STY = ZERO
C      FY = FINIT
C      DGY = DGINIT
CC
CC     START OF ITERATION.
CC
C   30 CONTINUE
CC
CC        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
CC        TO THE PRESENT INTERVAL OF UNCERTAINTY.
CC
C         IF (BRACKT) THEN
C            STMIN = MIN(STX,STY)
C            STMAX = MAX(STX,STY)
C         ELSE
C            STMIN = STX
C            STMAX = STP + XTRAPF*(STP - STX)
C            END IF
CC
CC        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
CC
C         STP = MAX(STP,STPMIN)
C         STP = MIN(STP,STPMAX)
CC
CC        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
CC        STP BE THE LOWEST POINT OBTAINED SO FAR.
CC
C         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
C     *      .OR. NFEV .GE. MAXFEV-1 .OR. INFOC .EQ. 0
C     *      .OR. (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)) STP = STX
CC
CC        EVALUATE THE FUNCTION AND GRADIENT AT STP
CC        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
CC        We return to main program to obtain F and G.
CC
C         DO 40 J = 1, N
C            X(J) = WA(J) + STP*S(J)
C   40       CONTINUE
C         INFO=-1
C         RETURN
CC
C   45    INFO=0
C         NFEV = NFEV + 1
C         DG = ZERO
C         DO 50 J = 1, N
C            DG = DG + G(J)*S(J)
C   50       CONTINUE
C         FTEST1 = FINIT + STP*DGTEST
CC
CC        TEST FOR CONVERGENCE.
CC
C         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
C     *      .OR. INFOC .EQ. 0) INFO = 6
C         IF (STP .EQ. STPMAX .AND.
C     *       F .LE. FTEST1 .AND. DG .LE. DGTEST) INFO = 5
C         IF (STP .EQ. STPMIN .AND.
C     *       (F .GT. FTEST1 .OR. DG .GE. DGTEST)) INFO = 4
C         IF (NFEV .GE. MAXFEV) INFO = 3
C         IF (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX) INFO = 2
C         IF (F .LE. FTEST1 .AND. ABS(DG) .LE. GTOL*(-DGINIT)) INFO = 1
CC
CC        CHECK FOR TERMINATION.
CC
C         IF (INFO .NE. 0) RETURN
CC
CC        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
CC        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
CC
C         IF (STAGE1 .AND. F .LE. FTEST1 .AND.
C     *       DG .GE. MIN(FTOL,GTOL)*DGINIT) STAGE1 = .FALSE.
CC
CC        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
CC        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
CC        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
CC        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
CC        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
CC
C         IF (STAGE1 .AND. F .LE. FX .AND. F .GT. FTEST1) THEN
CC
CC           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
CC
C            FM = F - STP*DGTEST
C            FXM = FX - STX*DGTEST
C            FYM = FY - STY*DGTEST
C            DGM = DG - DGTEST
C            DGXM = DGX - DGTEST
C            DGYM = DGY - DGTEST
CC
CC           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
CC           AND TO COMPUTE THE NEW STEP.
CC
C            CALL MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,
C     *                 BRACKT,STMIN,STMAX,INFOC)
CC
CC           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
CC
C            FX = FXM + STX*DGTEST
C            FY = FYM + STY*DGTEST
C            DGX = DGXM + DGTEST
C            DGY = DGYM + DGTEST
C         ELSE
CC
CC           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
CC           AND TO COMPUTE THE NEW STEP.
CC
C            CALL MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,
C     *                 BRACKT,STMIN,STMAX,INFOC)
C            END IF
CC
CC        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
CC        INTERVAL OF UNCERTAINTY.
CC
C         IF (BRACKT) THEN
C            IF (ABS(STY-STX) .GE. P66*WIDTH1)
C     *         STP = STX + P5*(STY - STX)
C            WIDTH1 = WIDTH
C            WIDTH = ABS(STY-STX)
C            END IF
CC
CC        END OF ITERATION.
CC
C         GO TO 30
CC
CC     LAST LINE OF SUBROUTINE MCSRCH.
CC
C      END
C      SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
C     *                 STPMIN,STPMAX,INFO)
C      INTEGER INFO
C      DOUBLE PRECISION STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
C      LOGICAL BRACKT,BOUND
CC
CC     SUBROUTINE MCSTEP
CC
CC     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
CC     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
CC     A MINIMIZER OF THE FUNCTION.
CC
CC     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
CC     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
CC     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
CC     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
CC     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
CC     WITH ENDPOINTS STX AND STY.
CC
CC     THE SUBROUTINE STATEMENT IS
CC
CC       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
CC                        STPMIN,STPMAX,INFO)
CC
CC     WHERE
CC
CC       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
CC         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
CC         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
CC         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
CC         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
CC
CC       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
CC         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
CC         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
CC         UPDATED APPROPRIATELY.
CC
CC       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
CC         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
CC         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
CC         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
CC
CC       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
CC         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
CC         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
CC         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
CC
CC       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
CC         AND UPPER BOUNDS FOR THE STEP.
CC
CC       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
CC         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
CC         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
CC         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
CC
CC     SUBPROGRAMS CALLED
CC
CC       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
CC
CC     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
CC     JORGE J. MORE', DAVID J. THUENTE
CC
C      DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
C      INFO = 0
CC
CC     CHECK THE INPUT PARAMETERS FOR ERRORS.
CC
C      IF ((BRACKT .AND. (STP .LE. MIN(STX,STY) .OR.
C     *    STP .GE. MAX(STX,STY))) .OR.
C     *    DX*(STP-STX) .GE. 0.0 .OR. STPMAX .LT. STPMIN) RETURN
CC
CC     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
CC
C      SGND = DP*(DX/ABS(DX))
CC
CC     FIRST CASE. A HIGHER FUNCTION VALUE.
CC     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
CC     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
CC     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
CC
C      IF (FP .GT. FX) THEN
C         INFO = 1
C         BOUND = .TRUE.
C         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
C         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
C         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
C         IF (STP .LT. STX) GAMMA = -GAMMA
C         P = (GAMMA - DX) + THETA
C         Q = ((GAMMA - DX) + GAMMA) + DP
C         R = P/Q
C         STPC = STX + R*(STP - STX)
C         STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)
C         IF (ABS(STPC-STX) .LT. ABS(STPQ-STX)) THEN
C            STPF = STPC
C         ELSE
C           STPF = STPC + (STPQ - STPC)/2
C           END IF
C         BRACKT = .TRUE.
CC
CC     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
CC     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
CC     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
CC     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
CC
C      ELSE IF (SGND .LT. 0.0) THEN
C         INFO = 2
C         BOUND = .FALSE.
C         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
C         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
C         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
C         IF (STP .GT. STX) GAMMA = -GAMMA
C         P = (GAMMA - DP) + THETA
C         Q = ((GAMMA - DP) + GAMMA) + DX
C         R = P/Q
C         STPC = STP + R*(STX - STP)
C         STPQ = STP + (DP/(DP-DX))*(STX - STP)
C         IF (ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN
C            STPF = STPC
C         ELSE
C            STPF = STPQ
C            END IF
C         BRACKT = .TRUE.
CC
CC     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
CC     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
CC     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
CC     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
CC     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
CC     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
CC     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
CC     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
CC
C      ELSE IF (ABS(DP) .LT. ABS(DX)) THEN
C         INFO = 3
C         BOUND = .TRUE.
C         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
C         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
CC
CC        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
CC        TO INFINITY IN THE DIRECTION OF THE STEP.
CC
C         GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
C         IF (STP .GT. STX) GAMMA = -GAMMA
C         P = (GAMMA - DP) + THETA
C         Q = (GAMMA + (DX - DP)) + GAMMA
C         R = P/Q
C         IF (R .LT. 0.0 .AND. GAMMA .NE. 0.0) THEN
C            STPC = STP + R*(STX - STP)
C         ELSE IF (STP .GT. STX) THEN
C            STPC = STPMAX
C         ELSE
C            STPC = STPMIN
C            END IF
C         STPQ = STP + (DP/(DP-DX))*(STX - STP)
C         IF (BRACKT) THEN
C            IF (ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN
C               STPF = STPC
C            ELSE
C               STPF = STPQ
C               END IF
C         ELSE
C            IF (ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN
C               STPF = STPC
C            ELSE
C               STPF = STPQ
C               END IF
C            END IF
CC
CC     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
CC     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
CC     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
CC     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
CC
C      ELSE
C         INFO = 4
C         BOUND = .FALSE.
C         IF (BRACKT) THEN
C            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
C            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
C            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
C            IF (STP .GT. STY) GAMMA = -GAMMA
C            P = (GAMMA - DP) + THETA
C            Q = ((GAMMA - DP) + GAMMA) + DY
C            R = P/Q
C            STPC = STP + R*(STY - STP)
C            STPF = STPC
C         ELSE IF (STP .GT. STX) THEN
C            STPF = STPMAX
C         ELSE
C            STPF = STPMIN
C            END IF
C         END IF
CC
CC     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
CC     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
CC
C      IF (FP .GT. FX) THEN
C         STY = STP
C         FY = FP
C         DY = DP
C      ELSE
C         IF (SGND .LT. 0.0) THEN
C            STY = STX
C            FY = FX
C            DY = DX
C            END IF
C         STX = STP
C         FX = FP
C         DX = DP
C         END IF
CC
CC     COMPUTE THE NEW STEP AND SAFEGUARD IT.
CC
C      STPF = MIN(STPMAX,STPF)
C      STPF = MAX(STPMIN,STPF)
C      STP = STPF
C      IF (BRACKT .AND. BOUND) THEN
C         IF (STY .GT. STX) THEN
C            STP = MIN(STX+0.66*(STY-STX),STP)
C         ELSE
C            STP = MAX(STX+0.66*(STY-STX),STP)
C            END IF
C         END IF
C      RETURN
CC
CC     LAST LINE OF SUBROUTINE MCSTEP.
CC
C      END
C
