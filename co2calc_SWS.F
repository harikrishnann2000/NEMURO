c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /home/orr/TeX/reports/noces/year2/howto/Cchem/mouchet/RCS/co2calc_SWS.f,v $ 
c_ $Revision: 1.5 $
c_ $Date: 2004/07/02 15:29:01 $   ;  $State: Exp $
c_ $Author: orr $ ;  $Locker: orr $
c_
c_ ---------------------------------------------------------------------
c_ $Log: co2calc_SWS.f,v $
c_ Revision 1.5  2004/07/02 15:29:01  orr
c_ Added missing parenthesis at end of "kf" equation (bug identified by Anne Mouchet)
c_
c_ Revision 1.4  2004/06/16 13:18:29  orr
c_ Initial Revision
c_
c_ ---------------------------------------------------------------------
c_ 
      subroutine co2calc_SWS(t,s,dic_in,ta_in,pt_in,sit_in, rho_in,
     &                   phlo,phhi,ph,xco2_in,atmpres,
     &                   co2star,dco2star,pCO2surf,dpco2,rff)
	real t,s,dic_in,ta_in,pt_in,sit_in, rho_in
	real phlo,phhi,ph,xco2_in,atmpres
	real co2star,dco2star,pCO2surf,dpco2
C
C-------------------------------------------------------------------------
C
C Modified from co2calc.f (RCS version 1.8, OCMIP-2) 
C - by A. Mouchet, 2004:
C
C NOW; "All" constants are given on seawater H scale (hSWS) 
C - HOWEVER, dissociation constants for S and F are on the 'free' H scale 
C            (this is necessary to work with hSWS)
C - this routine corrects for inconsistencies in the previous version.
C
C - Other changes:
C   * use ta_iter_SWS instead of ta_iter_1;
C   * hSWS replaces htotal;
C   * moved upward the calculation of bt, st and ft 
C     (needed in evaluation of kb);
C   * added various comments and references.
C
C
C SUBROUTINE CO2CALC_SWS
C
C PURPOSE
C	Calculate delta co2* from total alkalinity and total CO2 at
C temperature (t), salinity (s) and "atmpres" atmosphere total pressure. 
C
C USAGE
C       call co2calc_SWS(t,s,dic_in,ta_in,pt_in,sit_in
C    &                  ,phlo,phhi,ph,xco2_in,atmpres
C    &                  ,co2star,dco2star,pCO2surf,dpco2)
C
C INPUT
C	dic_in = total inorganic carbon (mol/m^3) 
C                where 1 T = 1 metric ton = 1000 kg
C	ta_in  = total alkalinity (eq/m^3) 
C	pt_in  = inorganic phosphate (mol/m^3) 
C	sit_in = inorganic silicate (mol/m^3) 
C	t      = temperature (degrees C)
C	s      = salinity (PSU)
C	phlo   = lower limit of pH range
C	phhi   = upper limit of pH range
C	xco2_in=atmospheric mole fraction CO2 in dry air (ppmv) 
C	atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
C
C       Note: arguments dic_in, ta_in, pt_in, sit_in, and xco2_in are 
C             used to initialize variables dic, ta, pt, sit, and xco2.
C             * Variables dic, ta, pt, and sit are in the common block 
C               "species".
C             * Variable xco2 is a local variable.
C             * Variables with "_in" suffix have different units 
C               than those without.

C OUTPUT
C	co2star  = CO2*water (mol/m^3)
C	dco2star = delta CO2 (mol/m^3)
c       pco2surf = oceanic pCO2 (ppmv)
c       dpco2    = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)
C
C IMPORTANT: Some words about units - (JCO, 4/4/1999)
c     - Models carry tracers in mol/m^3 (on a per volume basis)
c     - Conversely, this routine, which was written by observationalists 
c       (C. Sabine and R. Key), passes input arguments in umol/kg  
c       (i.e., on a per mass basis)
c     - I have changed things slightly so that input arguments are in mol/m^3,
c     - Thus, all input concentrations (dic_in, ta_in, pt_in, and st_in) 
c       should be given in mol/m^3; output arguments "co2star" and "dco2star"  
c       are likewise in mol/m^3.

C FILES and PROGRAMS NEEDED
C	drtsafe
C	ta_iter_SWS
C
C--------------------------------------------------------------------------
C
	real invtk,is,is2
        real k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi, rff
	real gamm, co2s, hg, cag, bohg, h3po4g, h2po4g, hpo4g
	real po4g, siooh3g, denom, dummy, fg, hSWS, hSWS2
	common /const/k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,ff,hSWS
	common /species/bt,st,ft,sit,pt,dic,ta
	external ta_iter_SWS
C

c       ---------------------------------------------------------------------
C       Change units from the input of mol/m^3 -> mol/kg:
c       (1 mol/m^3)  x (1 m^3/1024.5 kg)
c       where the oceans mean surface density is 1024.5 kg/m^3
c       Note: mol/kg are actually what the body of this routine uses 
c       for calculations.  
c       ---------------------------------------------------------------------
	permil = 1.0 / rho_in !1024.5
c       To convert input in mol/m^3 -> mol/kg 
	pt=pt_in*permil
	sit=sit_in*permil
	ta=ta_in*permil
	dic=dic_in*permil

c       ---------------------------------------------------------------------
C       Change units from uatm to atm. That is, atm is what the body of 
c       this routine uses for calculations.
c       ---------------------------------------------------------------------
	permeg=1.e-6
c       To convert input in uatm -> atm
	xco2=xco2_in*permeg
c       ---------------------------------------------------------------------
C
C*************************************************************************
C Calculate all constants needed to convert between various measured
C carbon species. References for each equation are noted in the code. 
C Once calculated, the constants are
C stored and passed in the common block "const". The original version of this
C code was based on the code by Dickson in Version 2 of "Handbook of Methods
C for the Analysis of the Various Parameters of the Carbon Dioxide System
C in Seawater", DOE, 1994 (SOP No. 3, p25-26). 
C
C Derive simple terms used more than once
C
	tk = 273.15 + t
	tk100 = tk/100.0
	tk1002=tk100*tk100
	invtk=1.0/tk
	dlogtk=log(tk)
	is=19.924*s/(1000.-1.005*s)
	is2=is*is
	sqrtis=sqrt(is)
	s2=s*s
	sqrts=sqrt(s)
	s15=s**1.5
	scl=s/1.80655
C
C------------------------------------------------------------------------
C Calculate concentrations for borate, sulfate, and fluoride
C
C Uppstrom (1974)
	bt = 0.000232 * scl/10.811
C Morris & Riley (1966)
	st = 0.14 * scl/96.062
C Riley (1965)
	ft = 0.000067 * scl/18.9984
C
C------------------------------------------------------------------------
C f = k0(1-pH2O)*correction term for non-ideality
C
C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
C
	ff = exp(-162.8301 + 218.2968/tk100  +
     &		90.9241*log(tk100) - 1.47696*tk1002 +
     &		s * (.025695 - .025225*tk100 + 
     &		0.0049867*tk1002))
C
C K0 from Weiss 1974
C
	k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) +
     &		s * (.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))

C
C------------------------------------------------------------------------
C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]     on hSWS
C
C Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale 
C (Original reference: Dickson and Millero, DSR, 1987)
C
	k1=10**(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk -
     &		0.0118 * s + 0.000116*s2))
C
	k2=10**(-1*(1394.7*invtk + 4.777 - 
     &		0.0184*s + 0.000118*s2))
C
C------------------------------------------------------------------------
C k1p = [H][H2PO4]/[H3PO4] on hSWS
C
C Millero p.670 (1995)
C
	k1p = exp(-4576.752*invtk + 115.540 - 18.453 * dlogtk +
     &		(-106.736*invtk + 0.69171) * sqrts +
     &		(-0.65643*invtk - 0.01844) * s)
C
C------------------------------------------------------------------------
C k2p = [H][HPO4]/[H2PO4] on hSWS
C
C Millero p.670 (1995)
C
	k2p = exp(-8814.715*invtk + 172.1033 - 27.927 * dlogtk +
     &		(-160.340*invtk + 1.3566) * sqrts +
     &		(0.37335*invtk - 0.05778) * s)
C
C------------------------------------------------------------------------
C k3p = [H][PO4]/[HPO4] on hSWS
C
C Millero p.670 (1995)
C
	k3p = exp(-3070.75*invtk - 18.126 + 
     &		(17.27039*invtk + 2.81197) *
     &		sqrts + (-44.99486*invtk - 0.09984) * s)
C
C------------------------------------------------------------------------
C ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
C
C Millero p.671 (1995) using data from Yao and Millero (1995)
C change to (mol/ kg soln)
C
	ksi = exp(-8904.2*invtk + 117.400 - 19.334 * dlogtk +
     &		(-458.79*invtk + 3.5913) * sqrtis +
     &		(188.74*invtk - 1.5998) * is +
     &		(-12.1652*invtk + 0.07871) * is2 +
     &		log(1.0-0.001005*s))
C
C------------------------------------------------------------------------
C kw = [H][OH] on hSWS
C
C Millero p.670 (1995) using composite data
C
	kw = exp(-13847.26*invtk + 148.9802 - 23.6521 * dlogtk +
     &		(118.67*invtk - 5.977 + 1.0495 * dlogtk) *
     &		sqrts - 0.01615 * s)
C
C------------------------------------------------------------------------
C ks = [H][SO4]/[HSO4] on free H scale
C
C Dickson (1990, J. chem. Thermodynamics 22, 113)
C change to (mol/ kg soln)
C
	ks=exp(-4276.1*invtk + 141.328 - 23.093*dlogtk +
     &		(-13856*invtk + 324.57 - 47.986*dlogtk) * sqrtis +
     &		(35474*invtk - 771.54 + 114.723*dlogtk) * is -
     &		2698*invtk*is**1.5 + 1776*invtk*is2 +
     &		log(1.0 - 0.001005*s))
C
C------------------------------------------------------------------------
C kf = [H][F]/[HF] on free H scale
C
C Dickson and Riley (1979)
C change to (mol/ kg soln)
C
	kf=exp(1590.2*invtk - 12.641 + 1.525*sqrtis +
     &		log(1.0 - 0.001005*s)) 
C
C------------------------------------------------------------------------
C kb = [H][BO2]/[HBO2] on hSWS
C
C Dickson p.673 (1990)
C change from htotal to hSWS
C
	kb=exp( (-8966.90 - 2890.53*sqrts - 77.942*s +
     &		1.728*s15 - 0.0996*s2)*invtk +
     &		(148.0248 + 137.1942*sqrts + 1.62142*s) +
     &		(-24.4344 - 25.085*sqrts - 0.2474*s) *
     &		dlogtk + 0.053105*sqrts*tk +
     &          log((1+(st/ks)+(ft/kf))/(1+(st/ks))) )
C
C*************************************************************************
C
C Calculate [H+] SWS when DIC and TA are known at T, S and 1 atm.
C The solution converges to err of xacc. The solution must be within
C the range x1 to x2.
C
C If DIC and TA are known then either a root finding or iterative method
C must be used to calculate hSWS. In this case we use the Newton-Raphson
C "safe" method taken from "Numerical Recipes" (function "rtsafe.f" with
C error trapping removed).
C
C As currently set, this procedure iterates about 12 times. The x1 and x2
C values set below will accomodate ANY oceanographic values. If an initial
C guess of the pH is known, then the number of iterations can be reduced to
C about 5 by narrowing the gap between x1 and x2. It is recommended that
C the first few time steps be run with x1 and x2 set as below. After that,
C set x1 and x2 to the previous value of the pH +/- ~0.5. The current
C setting of xacc will result in co2star accurate to 3 significant figures
C (xx.y). Making xacc bigger will result in faster convergence also, but this
C is not recommended (xacc of 10**-9 drops precision to 2 significant figures).
c ===========================================================================
c NEW CODE according to Follows et al., (2006), Ocean Modelling, 
c added on 21-Jan-2015

c First guess of [H+]: from last timestep *OR* fixed for cold start
	hg = 10.0**(-ph)
c estimate contributions to total alk from borate, silicate, phosphate
	bohg = bt*kb/(hg + kb)
	siooh3g = sit*ksi/(ksi + hg)
	denom = hg*hg*hg + (k1p*hg*hg) + k1p*k2p*hg + k1p*k2p*k3p
	h3po4g = (pt*hg*hg*hg)/denom
	h2po4g = (pt*k1p*hg*hg)/denom
	hpo4g = (pt*k1p*k2p*hg)/denom
	po4g = (pt*k1p*k2p*k3p)/denom
c estimate carbonate alkalinity
	fg = -bohg - (kw/hg) + hg - hpo4g - 2.0*po4g + h3po4g - siooh3g
	cag = ta + fg
c improved estimate of hydrogen ion conc
	gamm = dic/cag
	dummy = (1.0-gamm)*(1.0-gamm)*k1*k1 - 4.0*k1*k2*(1.0-2.0*gamm)
	hSWS = 0.5*((gamm-1.0)*k1 + sqrt(dummy))
	if ( hSWS .ne. hSWS) hSWS = 0.0
c ===========================================================================
C
C Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
C ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
C
C
C Parentheses added around negative exponents (Keith Lindsay)
C
	x1 = 10.0**(-phhi)
	x2 = 10.0**(-phlo)
	xacc = 1.e-10
	hSWS = drtsafe(ta_iter_SWS,x1,x2,xacc)
	if ( hSWS .ne. hSWS) hSWS = 0.0

	hSWS2=hSWS*hSWS
	rtemp = (hSWS2 + k1*hSWS + k1*k2)
	if (rtemp .eq. 0.0) rtemp = 1.0e0
c	co2star=dic*hSWS2/(hSWS2 + k1*hSWS + k1*k2)
	co2star=dic*hSWS2/rtemp
	co2starair=xco2*ff*atmpres
c	co2starair=xco2*k0*atmpres
	dco2star=co2starair-co2star
	if (rtemp .eq. 1.0) dco2star = 0.0 ! Date: 16-Mar-2015
	ph=-log10(hSWS)


c
c       ---------------------------------------------------------------
cc      Add two output arguments for storing pCO2surf
cc      Should we be using K0 or ff for the solubility here?
c       ---------------------------------------------------------------
        pCO2surf = co2star / ff
        dpCO2    = pCO2surf - xco2*atmpres
!! pCO2 to dic convertor; DIC = [ pCO2  ]*[rtemp*ff/hSWS2]. 27-Jan-2024
	if (hSWS2 .ne. 0.0) then
        rff = rtemp*ff/hSWS2
	else
	rff = 0.0
	endif
        
C
C  Convert units of output arguments
c      Note: co2star and dco2star are calculated in mol/kg within this routine 
c      Thus Convert now from mol/kg -> mol/m^3
       co2star  = co2star / permil
       dco2star = dco2star / permil

c      Note: pCO2surf and dpCO2 are calculated in atm above. 
c      Thus convert now to uatm
       pCO2surf = pCO2surf / permeg
       dpCO2    = dpCO2 / permeg
C
	return
	end
