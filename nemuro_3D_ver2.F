C     NEMURO 3-D biogeochemistry model.
	subroutine nemuro_3D (
     I  PS,PL,ZS,ZL,ZP,NO3,NH4,POM,DOM,SIOH4,Opal,Ca,CaCO3,TCO2,TALK,
     I  T0,S0,UVEL0,VVEL0, loop, dt, rlight,k,km,pomin,opalin,cain, caco3in,
     I  no3in,sioh4in,ppindx,
     I  wvelocity,tempin,saltin,uvelin,vvelin,t0_grad, s0_grad,
     I  uvel0_grad, vvel0_grad,no3_grad,POM_F_depth,Opal_F_depth,
     I  CaCO3_F_depth,no3_obss,sioh4_obss,
     I  RnewS,RnewL,NO3_lim,NH4_lim,N_lim,Si_lim,diff_N_and_Si_lim,
     I  no3_relax_obs,sioh4_relax_obs,sfx,
     I  temp_in, salt_in, kw660, ppo, dic_in, alk_in, po4, si, rhoin, xco2,
     I  dz1 , co2ex,ph,pCO2surf,ff,pco2Z,phZ,tco2in_rstr,talkin_rstr,
     I  rlxindx, ZLSUM,Onto)

#include "nemuro_new.h"

	real dz
	real(8) Monod, expPsi
	real(8) rkappa, rlight, rsum, sfx, day_night
	real(8) rtemp1, rtemp2
        real w1,w2,no31,no32,sioh41,sioh42,t01,t02,s01,s02,uvel01,uvel02
        real vvel01,vvel02
        real pom_grad(51),opal_grad(51),caco3_grad(51),no3_grad(51)
        real sioh4_grad(51)
        real t0_grad(51),s0_grad(51),uvel0_grad(51),vvel0_grad(51)
        real pomin(51),opalin(51),cain(51), caco3in(51),no3in(51),sioh4in(51)
        real pomin_sm(51),opalin_sm(51),caco3in_sm(51),no3in_sm(51), cain_sm(51)
        real sioh4in_sm(51)
        real wvelocity(51)
        real tempin_sm(51),saltin_sm(51),uvelin_sm(51),vvelin_sm(51)
        real tempin(51),saltin(51),uvelin(51),vvelin(51)
        real POM_F_depth,Opal_F_depth,CaCO3_F_depth,a
        real POM_F_Z(51),Opal_F_Z(51),CaCO3_F_Z(51),R
        real r_CP,d
        real w_diff,wb
        real delta,grad_POM, grad_Opal, grad_CaCO3
        real POM_ave
        real no3_obss,sioh4_obss
        real relaxation_NO3,relaxation_SIOH4
        real NO3_lim,NH4_lim,N_lim,Si_lim,diff_N_and_Si_lim
        real no3_relax_obs,sioh4_relax_obs
	integer loop, dt
        integer k,km,kmm,Z_c, kzc, kzcp

	real(8) PS_f, PL_f, ZS_f, ZL_f, ZP_f, NO3_f, NH4_f, POM_f, DOM_f,
     &             SIOH4_f, Opal_f, Ca_f, CaCO3_f, TCO2_f, TALK_f
        real(8) Ca_fTCO2, ZLSUM
        real    temp_in, salt_in, kw660, dic_in, alk_in
        real    si, rhoin, xco2, ppo, po4
        real    dz1 , co2ex,ph,pCO2surf,ff, ppindx(51)
        real talkin_rstr(51), tco2in_rstr(51), rlxindx(51)
        real pco2Z(51), phZ(51),Onto(51)

C     Grid resolution and model constants.
        kmm = km -3
        dz = dz1
        a   = 0.9

        R= 0.08
        r_CP = 117.0
        d = 3500.0

C     Primary production and nutrient limitation.
	PS_Photosynthesis =
     &  VmaxS * (Monod (NO3, KNO3S)*expPsi(PsiS, NH4)   +
     &           Monod (NH4, KNH4S))                 *
     &  exp(KS*T0)*(rlight/IoptS)*exp(1.0d0 - rlight/IoptS) * PS

	RnewS = Monod(Monod(NO3, KNO3S)*
     &   expPsi(PsiS, NH4) ,Monod(NH4, KNH4S))

	rtemp1 = Monod (NO3, KNO3L)*expPsi(PsiL, NH4)   +
     &               (Monod (NH4, KNH4L))

	rtemp2 = Monod(SIOH4,KSIL)

	PL_Photosynthesis = VmaxL *
     &  min (rtemp1, rtemp2)*
     &  exp(KL*T0)*(rlight/IoptL)*exp(1.0 - rlight/IoptL) * PL

        NO3_lim =  Monod (NO3, KNO3L)*expPsi(PsiL, NH4)
        NH4_lim= Monod (NH4, KNH4L)
        N_lim= rtemp1

        NO3_lim =  Monod (NO3, KNO3L)*expPsi(PsiL, NH4)
        NH4_lim= Monod (NH4, KNH4L)
        N_lim= rtemp1
        Si_lim = rtemp2

        NO3_lim =  Monod (NO3, KNO3L)*expPsi(PsiL, NH4)
        NH4_lim= Monod (NH4, KNH4L)
        N_lim= rtemp1
        Si_lim = rtemp2
        diff_N_and_Si_lim = rtemp1-rtemp2

	RnewL = Monod(Monod(NO3, KNO3L)*
     & expPsi(PsiL, NH4) ,Monod(NH4, KNH4L))

C     Respiration, mortality, grazing, and excretion.
	PS_Respiration = RPS0 * exp(KRS*T0) * PS

	PL_Respiration = RPL0 * exp(KRL*T0) * PL

	PS_Extracellular_Excretion = GammaS * PS_Photosynthesis

	PL_Extracellular_Excretion = GammaL * PL_Photosynthesis

	PS_Mortality = MPS0 * exp(KMPS*T0) * PS**2

	PL_Mortality = MPL0 * exp(KMPL*T0) * PL**2

	ZS_Mortality = MZS0 * exp(KZS*T0) * ZS**2

	ZL_Mortality = MZL0 * exp(KZL*T0) * ZL**2

	ZP_Mortality = MZP0 * exp(KZP*T0) * ZP**2

	rtemp1 = 1.0 - exp(LamdaS*(PSZS - PS))
	PS_Grazing_by_ZS = GRmaxS *
     &       max (0.0d0, rtemp1 ) *
     &       exp (KGS*T0) * ZS

	rtemp1 = 1.0 - exp(LamdaL*(PSZL - PS))
	PS_Grazing_by_ZL = GRmaxLPS *
     &       max (0.0d0, rtemp1 ) *
     &       exp (KGL*T0) * ZL

	rtemp1 = 1.0 - exp(LamdaL*(PLZL - PL))
	PL_Grazing_by_ZL =  GRmaxLPL *
     &       max (0.0d0, rtemp1 ) *
     &       exp (KGL*T0) * ZL

	rtemp1 = 1.0 - exp(LamdaL*(ZSZL - ZS))
	ZS_Predating_by_ZL = GRmaxLZS *
     &       max (0.0d0, rtemp1 ) *
     &       exp (KGL*T0) * ZL

	rtemp1 = 1.0 - exp(LamdaP*(PLZP - PL))
	PL_Grazing_by_ZP =   GRmaxPPL *
     &       max (0.0d0, rtemp1 ) *
     &       expPsi( PusaiPL, ZS+ZL)*exp(KGP*T0)*ZP

	rtemp1 = 1.0 - exp(LamdaP*(ZSZP - ZS))
	ZS_Predating_by_ZP = GRmaxPZS *
     &       max (0.0d0, rtemp1 ) *
     &       expPsi( PusaiZS, ZL)*exp(KGP*T0)*ZP

	rtemp1 = 1.0 - exp(LamdaP*(ZLZP - ZL))
	ZL_Predating_by_ZP = GRmaxPZL *
     &       max (0.0d0, rtemp1 ) *
     &       exp (KGP*T0) * ZP

        Ontogenetic_Migration = ZLSUM / (86400.0)

	ZS_Excretion = (GammaZS - BettaZS)* PS_Grazing_by_ZS

	ZL_Excretion = (GammaZL - BettaZL)* (PS_Grazing_by_ZL +
     &                 PL_Grazing_by_ZL + ZS_Predating_by_ZL )

	ZP_Excretion = (GammaZP - BettaZP)* (PL_Grazing_by_ZP +
     &                 ZS_Predating_by_ZP + ZL_Predating_by_ZP )

	ZS_Egestion = (1.0d0 - GammaZS) * PS_Grazing_by_ZS

	ZL_Egestion = (1.0d0 - GammaZL) * (PS_Grazing_by_ZL +
     &                 PL_Grazing_by_ZL + ZS_Predating_by_ZL)

	ZP_Egestion = (1.0d0 - GammaZP) * (PL_Grazing_by_ZP +
     &                 ZS_Predating_by_ZP + ZL_Predating_by_ZP)

C     Remineralization and elemental transformations.
	POM_Remineralization =  VPA0*exp(KPA*T0)*POM

	POM_Decomposition_to_DOM = VPD0*exp(KPD*T0)*POM

	DOM_Remineralization = VDA0*exp(KDA*T0)*DOM

	Opal_Decomposition = VOpal*exp(KOpal*T0)*Opal

	CaCO3_Decomposition =VCaCO3*exp(kCaCO3*T0)*CaCO3

	Nitrification = NNit0 * exp(KNit*T0)*NH4

	Si_PL_Shell_Formation = (PL_Photosynthesis - PL_Respiration
     &            - PL_Extracellular_Excretion)*RSiN

	Si_PL_Mortality = PL_Mortality * RSiN

	Si_ZL_Egestion = PL_Grazing_by_ZL * RSiN

	Si_ZP_Egestion = PL_Grazing_by_ZP * RSiN

        Ca_PS_Shell_Formation =( PS_Photosynthesis - PS_Respiration
     &       - PS_Extracellular_Excretion) * RCN * Rcoco * RCco / (RCN*RCco)

        Ca_ZS_Shell_Formation = (PS_Grazing_by_ZS) * BettaZS * RCN
     &       * Rfora * RCfo / (RCN * RCfo)

        Ca_PS_Mortality = (PS_Mortality) * RCN * Rcoco * RCco

        Ca_ZS_Mortality = (ZS_Mortality) * RCN * Rfora * RCfo

        Ca_ZS_Egestion = (PS_Grazing_by_ZS) * RCN *  Rcoco * RCco

        Ca_ZL_Egestion = (PS_Grazing_by_ZL) * RCN * Rcoco * RCco
     &                 + (ZS_Predating_by_ZL) * RCN * Rfora * RCfo

        Ca_ZP_Egestion = (ZS_Predating_by_ZP) * RCN * Rfora * RCfo
C     Ecosystem and nutrient tendency terms.
	PS_f = PS_Photosynthesis - PS_Respiration -
     &         PS_Extracellular_Excretion - PS_Mortality -
     &         PS_Grazing_by_ZS  - PS_Grazing_by_ZL

	PL_f = PL_Photosynthesis - PL_Respiration -
     &         PL_Extracellular_Excretion - PL_Mortality -
     &         PL_Grazing_by_ZL - PL_Grazing_by_ZP

	ZS_f = PS_Grazing_by_ZS - ZS_Excretion -
     &         ZS_Egestion - ZS_Mortality -
     &         ZS_Predating_by_ZL - ZS_Predating_by_ZP

	ZL_f = PS_Grazing_by_ZL + PL_Grazing_by_ZL +
     &         ZS_Predating_by_ZL - ZL_Excretion -
     &         ZL_Egestion - ZL_Mortality -
     &         ZL_Predating_by_ZP

	ZP_f = PL_Grazing_by_ZP + ZS_Predating_by_ZP +
     &         ZL_Predating_by_ZP - ZP_Excretion -
     &         ZP_Egestion - ZP_Mortality

	NO3_f = Nitrification   -
     &          (PS_Photosynthesis - PS_Respiration)*RnewS -
     &          (PL_Photosynthesis - PL_Respiration)*RnewL
	NH4_f =  (ZS_Excretion + ZL_Excretion + ZP_Excretion) *1.0    +
     &           (DOM_Remineralization + POM_Remineralization)*1.0 -
     &  Nitrification -
     & (PS_Photosynthesis - PS_Respiration)*(1.0d0 - RnewS) -
     & (PL_Photosynthesis - PL_Respiration)*(1.0d0 - RnewL)

	POM_f = PS_Mortality + PL_Mortality + ZS_Mortality +
     &          ZL_Mortality + ZP_Mortality + ZS_Egestion +
     &          ZL_Egestion + ZP_Egestion - POM_Remineralization -
     &          POM_Decomposition_to_DOM

	DOM_f = PS_Extracellular_Excretion + PL_Extracellular_Excretion +
     &          POM_Decomposition_to_DOM - DOM_Remineralization

	SIOH4_f = Opal_Decomposition - Si_PL_Shell_Formation

	Opal_f = Si_ZL_Egestion + Si_ZP_Egestion + Si_PL_Mortality -
     &           Opal_Decomposition

        Ca_f = CaCO3_Decomposition
     &  - Ca_ZS_Shell_Formation
     &  - Ca_PS_Shell_Formation

        Ca_fTCO2 = CaCO3_Decomposition
     &  - Ca_ZS_Shell_Formation * Rcoco * RCco * RCN
     &  - Ca_PS_Shell_Formation * Rfora * RCfo * RCN

        CaCO3_f = Ca_ZS_Egestion + Ca_ZL_Egestion + Ca_ZP_Egestion +
     &  Ca_PS_Mortality + Ca_ZS_Mortality - CaCO3_Decomposition

C     Carbonate chemistry and air-sea exchange.
        call co2flux(temp_in,salt_in,kw660,ppo,dic_in,alk_in,
     I   po4,si,rhoin,xco2,dz1,
     &   co2ex,ph,pCO2surf,ff)
        pco2Z(k) = pCO2surf
        phZ(k)   = ph

        if (k .eq. 1) then
	CO2_Air_Sea_Gas_Exchange = co2ex / 1025.0
        else
        CO2_Air_Sea_Gas_Exchange = 0.0
        endif
        TCO2_f = ((NO3_f + NH4_f))*RCN*ppindx(k) + Ca_fTCO2
     &          + CO2_Air_Sea_Gas_Exchange
     &          + (tco2in_rstr(k) - TCO2)/(365*86400*rlxindx(k))

        TALK_f = (2.0 * Ca_fTCO2) - (NO3_f - NH4_f)*ppindx(k)
     &          + (talkin_rstr(k) - TALK)/(365*86400*rlxindx(k))

C     Daily smoothing of forcing and transport profiles.
        do kk = k-1, k+1

        if (mod(loop,86400/dt).eq.0)then

        if ( k.ge.4 .and. k.le.kmm ) then

        pomin_sm(kk) =(1.0/16.0) * ((-1*pomin(kk-2))+(4*pomin(kk-1))
     &          +(10*pomin(kk))+(4*pomin(kk+1))+(-1*pomin(kk+2)))
        opalin_sm(kk)=(1.0/16.0) * ((-1*opalin(kk-2))+(4*opalin(kk-1))
     &        +(10*opalin(kk))+(4*opalin(kk+1))+(-1*opalin(kk+2)))

        caco3in_sm(kk)=(1.0/16.0)*((-1*caco3in(kk-2))+(4*caco3in(kk-1))
     &     +(10*caco3in(kk))+(4*caco3in(kk+1))+(-1*caco3in(kk+2)))

        cain_sm(kk)=(1.0/16.0)*((-1*cain(kk-2))+(4*cain(kk-1))
     &     +(10*cain(kk))+(4*cain(kk+1))+(-1*cain(kk+2)))

        no3in_sm(kk) =(1.0/16.0) *((-1*no3in(kk-2))+(4*no3in(kk-1))
     &     +(10*no3in(kk))+(4*no3in(kk+1))+(-1*no3in(kk+2)))

        sioh4in_sm(kk) =(1.0/16.0)*((-1*sioh4in(kk-2))+(4*sioh4in(kk-1))
     &    +(10*sioh4in(kk))+(4*sioh4in(kk+1))+(-1*sioh4in(kk+2)))

        tempin_sm(kk)=(1.0/16.0) * ((-1*tempin(kk-2))+(4*tempin(kk-1))
     &        +(10*tempin(kk))+(4*tempin(kk+1))+(-1*tempin(kk+2)))

        saltin_sm(kk)=(1.0/16.0) * ((-1*saltin(kk-2))+(4*saltin(kk-1))
     &        +(10*saltin(kk))+(4*saltin(kk+1))+(-1*saltin(kk+2)))

        uvelin_sm(kk)=(1.0/16.0) * ((-1*uvelin(kk-2))+(4*uvelin(kk-1))
     &        +(10*uvelin(kk))+(4*uvelin(kk+1))+(-1*uvelin(kk+2)))

        vvelin_sm(kk)=(1.0/16.0) * ((-1*vvelin(kk-2))+(4*vvelin(kk-1))
     &        +(10*vvelin(kk))+(4*vvelin(kk+1))+(-1*vvelin(kk+2)))

        pomin(kk)=pomin_sm(kk)
        opalin(kk)=opalin_sm(kk)
        caco3in(kk)=caco3in_sm(kk)
        cain(kk)=cain_sm(kk)
        no3in(kk) =no3in_sm(kk)
        sioh4in(kk)= sioh4in_sm(kk)
        tempin(kk) = tempin_sm(kk)
        saltin(kk) = saltin_sm(kk)
        uvelin(kk) = uvelin_sm(kk)
        vvelin(kk) = vvelin_sm(kk)

        end if
        end if
        end do

        if (mod(loop,86400/dt).eq.0)then
        if ( k.ge.4 .and. k.le.kmm ) then

        POM = pomin_sm(k)
        Opal= opalin_sm(k)
        CaCO3 = caco3in_sm(k)
        NO3 = no3in_sm(k)
        SIOH4 = sioh4in_sm(k)
        T0    = tempin_sm(k)
        S0    = saltin_sm(k)
        UVEL0 = uvelin_sm(k)
        VVEL0 = vvelin_sm(k)

        end if
        end if

C     Vertical gradients and transport terms.
        if (k.eq.1) then

        w1 = (wvelocity(k)+wvelocity(k))/2
        w2 = (wvelocity(k)+wvelocity(k+1))/2

        no31 = (no3in(k)+no3in(k))/2
        no32 = (no3in(k)+no3in(k+1))/2
        sioh41 = (sioh4in(k)+sioh4in(k))/2
        sioh42 = (sioh4in(k)+sioh4in(k+1))/2
        t01    = (tempin(k)+tempin(k))/2
        t02    = (tempin(k)+tempin(k+1))/2
        s01    = (saltin(k)+saltin(k))/2
        s02    = (saltin(k)+saltin(k+1))/2

        pom_grad(k) = (((pomin(k)+pomin(k))/2)
     &       -((pomin(k)+pomin(k+1))/2)) /(dz)

        opal_grad(k)=(((opalin(k)+opalin(k))/2)
     &       -((opalin(k)+opalin(k+1))/2))/(dz)

        caco3_grad(k)=(((caco3in(k)+caco3in(k))/2)
     &     -((caco3in(k)+caco3in(k+1))/2))/(dz)

        no3_grad(k) =  ((w1+w2)/2) * (no31 - no32)/(dz)
        sioh4_grad(k) = ((w1+w2)/2) * (sioh41 - sioh42)/(dz)
        t0_grad(k)    = ((w1+w2)/2) * (t01 - t02)/(dz)
        s0_grad(k)    = ((w1+w2)/2) * (s01 - s02)/(dz)

        else if ((k.ge.2) .and. (k.lt.km)) then
        w1 = (wvelocity(k)+wvelocity(k-1))/2
        w2 = (wvelocity(k)+wvelocity(k+1))/2

        no31 = (no3in(k)+no3in(k))/2
        no31 = (no3in(k)+no3in(k-1))/2
        no32 = (no3in(k)+no3in(k+1))/2
        sioh41 = (sioh4in(k)+sioh4in(k-1))/2
        sioh42 = (sioh4in(k)+sioh4in(k+1))/2
        t01    = (tempin(k)+tempin(k-1))/2
        t02    = (tempin(k)+tempin(k+1))/2
        s01    = (saltin(k)+saltin(k-1))/2
        s02    = (saltin(k)+saltin(k+1))/2

        pom_grad(k) = (((pomin(k) + pomin(k-1))/2)
     &      - ((pomin(k)+pomin(k+1))/2))/(dz)

        opal_grad(k) =(((opalin(k) + opalin(k-1))/2)
     &     - ((opalin(k)+opalin(k+1))/2))/(dz)
        caco3_grad(k) = (((caco3in(k) + caco3in(k-1))/2)
     &        -((caco3in(k)+caco3in(k+1))/2))/(dz)

        no3_grad(k) =  ((w1+w2)/2) * (no31 - no32)/(dz)
        sioh4_grad(k) = ((w1+w2)/2) * (sioh41 - sioh42)/(dz)
        t0_grad(k)    = ((w1+w2)/2)  * (t01 - t02)/(dz)
        s0_grad(k)    = ((w1+w2)/2) * (s01 - s02)/(dz)

        else
        w1 = (wvelocity(k-1)+wvelocity(k))/2
        w2 = (wvelocity(k)+wvelocity(k))/2

        no31 = (no3in(k)+no3in(k))/2
        no31 = (no3in(k-1)+no3in(k))/2
        no32 = (no3in(k)+no3in(k))/2
        sioh41 = (sioh4in(k-1)+sioh4in(k))/2
        sioh42 = (sioh4in(k)+sioh4in(k))/2
        t01   = (tempin(k-1)+tempin(k))/2
        t02   = (tempin(k)+tempin(k))/2
        s01   = (saltin(k-1)+saltin(k))/2
        s02   = (saltin(k)+saltin(k))/2

        pom_grad(k) = (((pomin(k-1)+pomin(k))/2)
     &      -((pomin(k)+ pomin(k))/2))/(dz)
        opal_grad(k) = (((opalin(k-1)+opalin(k))/2)
     &         - ((opalin(k)+opalin(k))/2))/(dz)
        caco3_grad(k) = (((caco3in(k-1)+caco3in(k))/2)
     &    - ((caco3in(k) +caco3in(k))/2) )/(dz)

        no3_grad(k) =  ((w1+w2)/2) * (no31 - no32)/(dz)
        sioh4_grad(k) = ((w1+w2)/2) * (sioh41 - sioh42)/(dz)
        t0_grad(k)    = ((w1+w2)/2)  * (t01 - t02)/(dz)
        s0_grad(k)    = ((w1+w2)/2) * (s01 - s02)/(dz)

        end if
        ExpPON   = SPOM  * pom_grad(k)
        ExpOpal  = SOpal * opal_grad(k)
        ExpCaCO3 = SCaCO3  * caco3_grad(k)

        ExcNO3   = no3_grad(k)

        ExcSiOH4 = sioh4_grad(k)
        ExcT0    = t0_grad(k)
        ExcS0    = s0_grad(k)

C     Depth-dependent particle export fluxes.
         Z_c = 55.0
         kzc = 51
         kzcp = kzc + 1

         if(k .le. kzc) then

         grad_POM = pomin(k)/dz
         grad_Opal = opalin(k)/dz
         grad_CaCO3 = caco3in(k)/dz

         POM_f = POM_f - (grad_POM * SPOM)
         Opal_f = Opal_f - (grad_Opal * SOpal )
         CaCO3_f = CaCO3_f - (caco3_grad(k) * SCaCO3)

         POM_F_depth = POM_F_depth  + max(0.0, POM) * dz
         Opal_F_depth = Opal_F_depth + max(0.0,Opal) * dz

         POM_F_Z(k)   = 0.0
         Opal_F_Z(k)  = 0.0
         endif

         POM_F_depth = max(0.0, POM_F_depth)
         Opal_F_depth = max(0.0, Opal_F_depth)

         if((k.ge.kzc).and.(k.le.km)) then
         POM_F_Z(k) = POM_F_depth *  ((k*dz)/Z_c)**(-0.9)

         Opal_F_Z(k) = Opal_F_depth *  ((k*dz)/Z_c)**(-0.9)

         end if

         if(k.ge.kzcp) then
         POM_f = POM_f    + ((POM_F_Z(k-1)-POM_F_Z(k))/dz)
         Opal_f = Opal_f  + ((Opal_F_Z(k-1)-Opal_F_Z(k))/dz)

         end if

C     Nutrient relaxation and tendency corrections.
         if (k.le.8) then
         relaxation_NO3 =(NO3 - no3_obss)/86400.0/180.0
         relaxation_SIOH4=(SIOH4 - sioh4_obss)/86400.0/180.0
         end if

         if ((k.gt.8).and.(k.le.16)) then
         relaxation_NO3 = 5*(NO3 - no3_obss)/86400.0/90.0
         relaxation_SIOH4 = 5*(SIOH4 - sioh4_obss)/86400.0/90.0
         end if
         if (k.gt.16) then
         relaxation_NO3 = 5*(NO3 - no3_obss)/86400.0/30.0
         relaxation_SIOH4 =5*(SIOH4 - sioh4_obss)/86400.0/30.0
         end if

         NO3_f   = NO3_f - (ExcNO3*1.0)  - relaxation_NO3*1.0
         SIOH4_f = SIOH4_f - (ExcSiOH4*1.0) -  relaxation_SIOH4*1.0

         POM_f   = POM_f   - ExpPON *0.0
         Opal_f  = Opal_f  - ExpOpal *0.0
	 CaCO3_f = CaCO3_f - ExpCaCO3 *0.0

C     Update ecosystem state variables.
	PS = max(0.0, PS + PS_f*dt)
	PL = max(0.0, PL + PL_f*dt)
	ZS = max(0.0, ZS + ZS_f*dt)
	ZL = max(0.0, ZL + ZL_f*dt)
	ZP = max(0.0, ZP + ZP_f*dt)
	NO3 = max (0.0 , NO3 + NO3_f*dt)
	NH4 = max (0.0, NH4 + NH4_f*dt)
	POM = max (0.0, POM + POM_f*dt)
	DOM = max(0.0, DOM + DOM_f*dt)
	SIOH4 = max(0.0, SIOH4 + SIOH4_f*dt)
	Opal = max(0.0, Opal + Opal_f*dt)
        Ca = max(0.0, Ca + Ca_f*dt)
        CaCO3 = max(0.0, CaCO3 + CaCO3_f * dt)
        TCO2 = max(0.0, TCO2 + TCO2_f * dt)
        TALK = max(0.0, TALK + TALK_f * dt)
        T0 = T0 - (ExcT0 * dt)
        S0 = S0 - (ExcS0 * dt)

	return
	end

C     Supporting response functions.
	function Monod (a, b)
	real(8) a, b, Monod
	Monod = a/(a+b)
	end

	function expPsi (a, b)
	real(8) expPsi, a, b
	expPsi = exp(-1.0*a*b)
	end
