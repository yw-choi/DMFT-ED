module ed_projection

contains

    function func(x)
        implicit none
        double precision :: x()
        double precision :: func
    end function func

    subroutine dfunc(x,df)
        implicit none
        double precision :: x,df
    end subroutine dfunc

!     The following routines are from the numerical recipes.
!     Given a starting point p that is a vector of length n, Fletcher-Reeves-Polak-Ribiere minimization
!     is performed on a function func, using its gradient as calculated by a routine dfunc.
!     The convergence tolerance on the function value is input as ftol. Returned quantities are p 
!     (the location of the minimum), iter ( the # of iterations that were performed), and fret (the minimum value of function)
!     The routine linmin is called to perform line minimizations
!     parameter: nmax is the maximum anticipated value of n; itmax is the maximum allowd number of iterations;
!     eps is a small number of rectify special case of converging to exactly zero function value.
      SUBROUTINE FRPRMN(P,N,FTOL,ITER,FRET)
      PARAMETER (NMAX=50,ITMAX=200,EPS=1.E-10)
      DOUBLE PRECISION P(N),G(NMAX),H(NMAX),XI(NMAX)
      FP=FUNC(P)
      CALL DFUNC(P,XI)
      DO 11 J=1,N
        G(J)=-XI(J)
        H(J)=G(J)
        XI(J)=H(J)
11    CONTINUE
      DO 14 ITS=1,ITMAX
        ITER=ITS
        CALL LINMIN(P,XI,N,FRET)
        IF(2.*ABS(FRET-FP).LE.FTOL*(ABS(FRET)+ABS(FP)+EPS))RETURN
        FP=FUNC(P)
        CALL DFUNC(P,XI)
        GG=0.
        DGG=0.
        DO 12 J=1,N
          GG=GG+G(J)**2
!         DGG=DGG+XI(J)**2
          DGG=DGG+(XI(J)+G(J))*XI(J)
12      CONTINUE
        IF(GG.EQ.0.)RETURN
        GAM=DGG/GG
        DO 13 J=1,N
          G(J)=-XI(J)
          H(J)=G(J)+GAM*H(J)
          XI(J)=H(J)
13      CONTINUE
14    CONTINUE
      PAUSE 'FRPR maximum iterations exceeded'
      RETURN
      END SUBROUTINE FRPRMN

      SUBROUTINE LINMIN(P,XI,N,FRET)
      PARAMETER (NMAX=50,TOL=1.E-4)
      ! EXTERNAL F1DIM
      DOUBLE PRECISION P(N),XI(N)
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)
      NCOM=N
      DO 15 J=1,N
        PCOM(J)=P(J)
        XICOM(J)=XI(J)
15    CONTINUE
      AX=0.
      XX=1.
      BX=2.
      CALL MNBRAK(AX,XX,BX,FA,FX,FB,F1DIM)
      FRET=BRENT(AX,XX,BX,F1DIM,TOL,XMIN)
      DO 16 J=1,N
        XI(J)=XMIN*XI(J)
        P(J)=P(J)+XI(J)
16    CONTINUE
      RETURN
      END SUBROUTINE LINMIN

      SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUNC)
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.E-20)
      FA=FUNC(AX)
      FB=FUNC(BX)
      IF(FB.GT.FA)THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BX+GOLD*(BX-AX)
      FC=FUNC(CX)
17    IF(FB.GE.FC)THEN
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            GO TO 17
          ELSE IF(FU.GT.FB)THEN
            CX=U
            FC=FU
            GO TO 17
          ENDIF
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
            FU=FUNC(U)
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
          U=ULIM
          FU=FUNC(U)
        ELSE
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
        GO TO 17
      ENDIF
      RETURN
      END SUBROUTINE MNBRAK

      FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN)
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0E-10)
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      FX=F(X)
      FV=FX
      FW=FX
      DO 18 ITER=1,ITMAX
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 23
        IF(ABS(E).GT.TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.*(Q-R)
          IF(Q.GT.0.) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR.P.GE.Q*(B-X)) GOTO 25
          D=P/Q
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 24
        ENDIF
25      IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
24      IF(ABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=F(U)
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
18    CONTINUE
      PAUSE 'Brent exceed maximum iterations.'
23    XMIN=X
      BRENT=FX
      RETURN
      END FUNCTION BRENT

      FUNCTION F1DIM(X)
      PARAMETER (NMAX=50)
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)
      DOUBLE PRECISION XT(NMAX)
      DO 30 J=1,NCOM
        XT(J)=PCOM(J)+X*XICOM(J)
30    CONTINUE
      F1DIM=FUNC(XT)
      RETURN
      END FUNCTION F1DIM
end module ed_projection

