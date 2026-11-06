c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /www/html/ipsl/OCMIP/phase2/simulations/Biotic/boundcond/RCS/co2flux.f,v $ 
c_ $Revision: 1.6 $
c_ $Date: 1999/11/25 09:07:10 $   ;  $State: Exp $
c_ $Author: orr $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: co2flux.f,v $
c_ Revision 1.6  1999/11/25 09:07:10  orr
c_ Global change pp0 -> ppo (for consistency with other Biotic routines)
c_
c_ Revision 1.5  1999/11/25 08:52:16  orr
c_ Fixed error with ppo vs pp0 inconsistency (global change: ppo -> pp0)
c_ Added IMPLICIT NONE and declared all vars to avoid future problems
c_
c_ Revision 1.4  1999/04/06 13:19:34  orr
c_ Added two arguments to co2calc.f: pCO2surf and dpCO2 (see HOWTO)
c_
c_ Revision 1.3  1999/04/05 15:57:42  orr
c_ Added comments concerning input/output units
c_ Input for co2calc is now mol/m^3.
c_
c_ Revision 1.2  1999/04/04 02:46:19  orr
c_ Fixed erroneous CO2 piston veocity calculation;
c_ Updated arguments for call co2calc, based on new release of that
c_ routine from Sabine and Key.
c_
c_ Revision 1.1  1999/03/22 12:57:47  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_ 
      subroutine co2flux(t,s,kw660,ppo,dic,alk,po4,si,rhoin,xco2,dz1,
     &   co2ex,ph,pCO2surf,ff)
c
c**********************************************************************
c
c  Computes the time rate of change of DIC in the surface
c  layer due to air-sea gas exchange in mol/m^3/s.
c
c  Inputs:
c    t        model surface temperature (deg C)
c    s        model surface salinity (permil)
c    kw660    gas transfer velocity at a Schmidt number of 660, accounting
c               for sea ice fraction (m/s)
c    ppo      surface pressure divided by 1 atm
c    dic      surface DIC concentration (mol/m^3)
c    alk      surface alkalinity (eq/m^3)
c    po4      surface phosphate concentration (mol/m^3)
c    si       surface silicate concentration (mol/m^3)
c    xco2     atmospheric CO2 mixing ratio (ppm)
c    dz1      surface grid box thickness (m)
c  Output:
c    co2ex    time rate of change of DIC in the surface layer due
c               to air-sea exchange (mol/m^3/s)
c
c  Two functions are called:
c    scco2    Schmidt number of CO2
cc   co2sato  CO2 saturation concentration at 1 atm (mol/m^3)
c  One subroutine is called:
c    co2calc  Computes actual aqueous CO2 concentration in surface
c             layer (mol/m^3), the saturation concentration at a
c             given atmospheric pressure, and there difference. 
c
c  Numbers in brackets refer to equation numbers in simulation design
c  document.
c
c  Ray Najjar 1/29/99
c  Modifications: James Orr 4/4/99 
c                 (in particular, see changes near "cc" prefix)
c
c**********************************************************************
c
      IMPLICIT none

c     Declare arguments to this routine
      REAL t, s, kw660, ppo, dic, alk, po4, si, xco2, dz1, co2ex
      REAL rhoin, ff

c     Declare FUNCTION type
      REAL ScCO2

c     Declare output arguments to co2calc
      REAL ph, co2star, dco2star, pCO2surf, dpco2

c     Declare local variables
      real kwco2
      REAL phlo, phhi
      REAL dicin, alkin, po4in, siin

c  Compute the transfer velocity for CO2 in m/s [4]
c
cc    Somethings seems to be missing here, and units change is unneccessary
cc    because OCMIP gas exchange fields (i.e., xKw) is already in m/s!
cc    kwco2 = (scco2(t)/660)**-0.5*0.01/3600.0
      kwco2 =          kw660   !* (660/ScCO2(t))**0.5   ! The Sc-Number is already included
c
c  ---------------------------------------------------------------------------
c  >>> RGN: Do the carbonate equilibrium calculations.  Pick values of
c  >>>      ph range.  You may want to change these per suggestions of Sabine
c  >>>      and Key.  Note also that co2calc subroutine expects 
c  >>>      concentrations in mol/kg and eq/kg, so average surface density 
c  >>>      is used to convert.
C
c  >>> JCO: Yes, but one can pass arguments in mol/m^3 directly,
c  >>>      IF input and output arguments are consistent (i.e., always in
c  >>>      per volume basis, not per mass). See my comments in modified 
c  >>>      co2calc.f.  Related errors are small.

      phlo = 6.0
      phhi = 9.0

cc    dicin = dic/rhoavg
cc    alkin = alk/rhoavg
cc    po4in = po4/rhoavg
cc    siin = si/rhoavg
      dicin = dic
      alkin = alk
      po4in = po4
      siin = si
c  ---------------------------------------------------------------------------

cc    call co2calc(t,s,dicin,alkin,po4in,siin,phlo,phi,ph,co2star)
        phlo = 5.0
        phhi = 9.0
      CALL co2calc_SWS(t,s,dicin,alkin,po4in,siin,rhoin,phlo,phhi,
     &            ph,xco2,ppo,
     &            co2star,dco2star,pCO2surf,dpco2,ff)

c  Compute the saturation concentration for CO2 [3]
c
cc    co2sat = co2sato(t,s,xco2)*ppo
c
c  Compute time rate of change of CO2 due to gas exchange [1]
c
cc    co2ex = kwco2*(co2sat-co2star)/dz1
      co2ex = kwco2*(dco2star)/dz1

cc    Further handling of output arguments (1) pCO2surf (pCO2ocn) and 
cc    (2) dpCO2 (pCO2ocn - pCO2atm) are left to the discretion of the 
cc    the modeler.  Both arguments are needed for standard OCMIP-2 output. 
c
      return
      END







