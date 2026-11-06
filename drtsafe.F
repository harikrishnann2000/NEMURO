c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /www/html/ipsl/OCMIP/phase2/simulations/Abiotic/Cchem/RCS/drtsafe.f,v $ 
c_ $Revision: 1.1 $
c_ $Date: 1999/04/03 22:00:42 $   ;  $State: Exp $
c_ $Author: orr $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: drtsafe.f,v $
c_ Revision 1.1  1999/04/03 22:00:42  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_ 
      REAL FUNCTION DRTSAFE(FUNCD,X1,X2,XACC)
C
C	File taken from Numerical Recipes. Modified  R.M.Key 4/94
C
      MAXIT=100
      CALL FUNCD(X1,FL,DF)
      CALL FUNCD(X2,FH,DF)
      IF(FL .LT. 0.0) THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        SWAP=FL
        FL=FH
        FH=SWAP
      END IF
      DRTSAFE=.5*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL FUNCD(DRTSAFE,F,DF)
      DO 100, J=1,MAXIT
        IF(((DRTSAFE-XH)*DF-F)*((DRTSAFE-XL)*DF-F) .GE. 0.0 .OR.
     &	      ABS(2.0*F) .GT. ABS(DXOLD*DF)) THEN
          DXOLD=DX
          DX=0.5*(XH-XL)
          DRTSAFE=XL+DX
          IF(XL .EQ. DRTSAFE)RETURN
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=DRTSAFE
          DRTSAFE=DRTSAFE-DX
          IF(TEMP .EQ. DRTSAFE)RETURN
	END IF
        IF(ABS(DX) .LT. XACC)RETURN
        CALL FUNCD(DRTSAFE,F,DF)
        IF(F .LT. 0.0) THEN
          XL=DRTSAFE
          FL=F
        ELSE
          XH=DRTSAFE
          FH=F
        END IF
  100  CONTINUE
      RETURN
      END
