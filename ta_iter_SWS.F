c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /home/orr/TeX/reports/noces/year2/howto/Cchem/mouchet/RCS/ta_iter_SWS.f,v $ 
c_ $Revision: 1.4 $
c_ $Date: 2004/06/16 13:18:29 $   ;  $State: Exp $
c_ $Author: orr $ ;  $Locker: orr $
c_
c_ ---------------------------------------------------------------------
c_ $Log: ta_iter_SWS.f,v $
c_ Revision 1.4  2004/06/16 13:18:29  orr
c_ Initial Revision
c_
c_ ---------------------------------------------------------------------
c_ 
	subroutine ta_iter_SWS(x,fn,df)
	real k12,k12p,k123p
        real k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
	common /const/k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,ff,hSWS
	common /species/bt,st,ft,sit,pt,dic,ta
C
C Modified from ta_iter_1.f (RCS version 1.2, OCMIP-2)
C - by A. Mouchet, 2004:
C Fixed Problems w/ version of ta_iter_1.f used in OCMIP-2 (vers. 1.2)
C  1) fixed errors in signs, parenthesis and coefficient c in derivative
C  2) changed from Total to Seawater Scale 
C     * c defined for seawater H scale; 
C     * fn and df adapted to KF on free H scale
C     * comments have been adapted
C

C
C This routine expresses TA as a function of DIC, hSWS and constants.
C It also calculates the derivative of this function with respect to 
C hSWS. It is used in the iterative solution for hSWS. In the call
C "x" is the input value for hSWS, "fn" is the calculated value for TA
C and "df" is the value for dTA/dhSWS
C
	x2=x*x
	x3=x2*x
	k12 = k1*k2
	k12p = k1p*k2p
	k123p = k12p*k3p
	c = 1.0 + st/ks + ft/kf
	a = x3 + k1p*x2 + k12p*x + k123p
	a2=a*a
	da = 3.0*x2 + 2.0*k1p*x + k12p
	b = x2 + k1*x + k12
	b2=b*b
	db = 2.0*x + k1
C
C	fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
C===========================================================================
C
	fn = k1*x*dic/b +
     &	     2.0*dic*k12/b +
     &	     bt/(1.0 + x/kb) +
     &	     kw/x +
     &	     pt*k12p*x/a +
     &	     2.0*pt*k123p/a +
     &	     sit/(1.0 + x/ksi) -
     &	     x/c -
     &	     st/(1.0 + ks/(x/c)) -
     &	     ft/(1.0 + kf/(x/c)) -
     &	     pt*x3/a -
     &	     ta
C
C	df = dfn/dx
C
	df = ((k1*dic*b) - k1*x*dic*db)/b2 -
     &	     2.0*dic*k12*db/b2 -
     &	     bt/kb/(1.0+x/kb)**2. -
     &	     kw/x2 +
     &	     (pt*k12p*(a - x*da))/a2 -
     &	     2.0*pt*k123p*da/a2 -
     &	     sit/ksi/(1.0+x/ksi)**2. -
     &	     1.0/c -
     &	     st *(1.0 + ks/(x/c))**(-2.0) *(ks*c/x2) -
     &	     ft*(1.0 + kf/(x/c))**(-2.0) *(kf*c/x2) -
     &	     pt*x2*(3.0*a-x*da)/a2
C
	return
	end
