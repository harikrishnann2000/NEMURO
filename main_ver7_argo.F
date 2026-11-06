	program main
!
!	module to test nemuro_3D.F on subroutine mode
!
!

#include "size.h"
	real(8) PS, PL, ZS, ZL, ZP
        real(8) NO3, NH4, POM, DOM, SIOH4, Opal
        real(8) Ca, CaCO3, TCO2, TALK                               
        real(8) T0,day_night,sfx, rlight, rsum1
        real(8) S0,UVEL0,VVEL0
        real(8) co2 , ZLOnt(0:52), ZLSUM !NEWWW


	real tempin(51), saltin(51), uvelin(51), vvelin(51)
	real psin(51),plin(51),zsin(51),zlin(51),zpin(51)
	real no3in(51),nh4in(51),pomin(51),domin(51)
	real sioh4in(51),opalin(51)
	real cain(51),caco3in(51)
	real tco2in(51),talkin(51)
	real tco2in_rstr(51),talkin_rstr(51)
	real SMout(51), SHout(51)
        real(4)ExcT0_var(51),ExcS0_var(51),ExcUVEL0_var(51)
        real(4)ExcVVEL0_var(51),ExcNO3_var(51) 
        real t0_grad(51),s0_grad(51),uvel0_grad(51),vvel0_grad(51)
        real no3_grad(51)
        real rlightk(51),Pa, qsat, ql, Lv, RO, cp_a, qout
        real(4)tauxin,tauyin !doubt want to define or not
        real(4)Tair1,Tair2,Tair,uwnd1,uwnd2,uwnd,SST_in,Onto1(51),Onto2(51),Onto(51)
        real(4)vwnd1,vwnd2,vwnd,sfx_in1,sfx_in2,sfx_in
        real(4)co2_in1,co2_in2,co2_in


        real(4)pme,ppt,evp,ppt1,ppt2,curl_tau1,curl_tau2,curl_tau,evp1
        real(4)evp2,spe_hum1,spe_hum2,hum1,hum2,hum
        real(4)velocity(51),velocitya(51),velocityb(51)
        real(4)cp_w
        real(4)rho_a,cs,cl,emiss,sigma_const,c1,C,vel         
        real(4)qs,dq,lw_term1,lw_term2,lw_term3,lw           
        real(4)Tair_Ke,salt
        real(4)omega 
        real(4)wvel(615,51),wvel1(51),wvel2(51)!wvel(122,51),wvel1(51),wvel2(51)
        real(4)wvelocity(51)
        real(4)temp_var(51),salt_var(51),ps_var(51),pl_var(51)
        real(4)rlight_var(51)
        real(4)zs_var(51),zl_var(51),zp_var(51),no3_var(51)
        real(4)nh4_var(51),pom_var(51),dom_var(51),sioh4_var(51)
        real(4)opal_var(51),ca_var(51),caco3_var(51),tco2_var(51)
        real(4)talk_var(51)
        real(4)uvel_var(51),vvel_var(51)
        real(4)wvelocity_var(51), ph_var(51)
        real(4)co2ex_var(51),pCO2surf_var(51)
        real(4)xco2_var(51)

        real(4)pme_offset,pme_offset1
        real(4)POM_F_depth,Opal_F_depth,CaCO3_F_depth
        real(4) sswin(51)
        real(4) lw_sst,lw_air,alpha
        real(4) long,year,pme_sum,year_no
        real(4) pme_int
        real(4)sal_clim,sal_clim1,sal_clim2,no3_obs(51),no3_obss
        real(4)sioh4_obs(51),sioh4_obss

        !WHY REAL 4 GIVING nan VALUES
        real(8)RnewL,RnewS,RnewL_d(51),RnewS_d(51),RnewS_var(51)
        real(8)RnewL_var(51)
        !*******************************88
        real NO3_lim,NH4_lim,N_lim,Si_lim,diff_N_and_Si_lim
        real NO3_lim_d(51),NH4_lim_d(51),N_lim_d(51),Si_lim_d(51)
        real diff_N_and_Si_lim_d(51)
        real NO3_lim_var(51),NH4_lim_var(51),N_lim_var(51)
        real Si_lim_var(51),diff_N_and_Si_lim_var(51)
        real no3_relax_obs(51),sioh4_relax_obs(51)
        real   date
        real(8) PSCChla(51),PLCChla(51),PSCChla_var(51), PLCChla_var(51)
        real(8) C_atm_mass, Red_rat, C_chla_rat,Chloro
        real(8) KWA,KCA,KDA,KAA,mlight
	integer dt, taun 
        integer lpd,lsd,lkd,lmd,irec,lco2,lOnt
        integer day_no,d,lo,lo_n 
        real  temp_in, salt_in, kw660, dic_in, alk_in
        real  si, rhoin, xco2, ppo, po4
        real  dz1 , co2ex1,ph1,pCO2surf,ff, ppindx(51),nonlocal(51),rlxindx(51)
        real pco2Z(51), phZ(51)
	
#include "param.h"
	
	i = 1
	j = 1
	k = 1

	month = 1
	month_start_snap = month

	do k=1,km
	zdepth(k) = -1.0*(2.5*k)
        dz(k) = 5.0
	enddo

        open (11,file='ppindx.txt')
        do l=1,km
        read (11,*) ppindx(l)
        enddo
        close (11)

        open (11,file='nonlocal.txt')
        do l=1,km
        read (11,*) nonlocal(l)
        enddo
        close (11)

        open (11,file='rlxindx.txt')
        do l=1,km
        read (11,*) rlxindx(l)
        enddo
        close (11)

	taun = 1
	open(20,file='INPUT/details_argo.txt')
        read(20,*)long
        read(20,*)rLATin
        close(20)
        
!**********DATA************************************************************
!        temp and salinity from levitus_climatology depth interpolated 

	open (1,file='INPUT/temp_in_data.txt')
        open (8,file='INPUT/salt_in_data.txt')        
 
        do k=1,km
        read (1,*) tempin(k)
        enddo
        close (1)
        
        do k=1,km
        read (8,*) saltin(k)
        enddo
        close (8)
        
        open(20,file='INPUT/NO3_in_data.DAT')
        open(21,file='INPUT/SIOH4_in_data.DAT') 
        open(22,file='INPUT/TCO2_in_data_for_argo11.DAT')
        open(23,file='INPUT/TALK_in_data_for_argo11.DAT')
        do k=1,km
        read(20,*) no3_obs(k)
        read(21,*) sioh4_obs(k)
        read(22,*) tco2in(k)
        read(23,*) talkin(k)
        tco2in_rstr(k) = tco2in(k)
        talkin_rstr(k) = talkin(k)
        end do
        
        close(20)
        close(21)
        close(22)
        close(23)
!**************************************************************************
	!initial condition        
        do k=1,km

	u(i,k,j,taun) = 0.0
	v(i,k,j,taun) = 0.0
	t(i,k,j,1,taun) = tempin(k)                   ! TEMP,20.0-0.01*(M-1)
	t(i,k,j,2,taun) = saltin(k)                   !SALT,33.0 + 0.05*k 
	t(i,k,j,3,taun) = 0.1e-6                      !0.1e-6  ! PS !molN/L
	t(i,k,j,4,taun) = 0.1e-6                      !0.1e-6 ! PL
	t(i,k,j,5,taun) = 0.1e-6                      !0.1e-6 ! ZS
	t(i,k,j,6,taun) = 0.1e-6                      !0.1e-6 ! ZL
	t(i,k,j,7,taun) = 0.1e-6                      !0.1e-6 ! ZP
	t(i,k,j,8,taun) = no3_obs(k)*1e-6             !5.0e-6  ! NO3
	t(i,k,j,9,taun) = 0.5e-6                      !0.5e-6 ! NH4
	t(i,k,j,10,taun) = 0.1e-6                     ! POM
	t(i,k,j,11,taun) = 0.1e-6                     ! DOM
	t(i,k,j,12,taun) = sioh4_obs(k)*1e-6          !5.0e-6  !SIOH4 !check this 3e-6 or 5e-6
	t(i,k,j,13,taun) = 0.1e-6                     ! Opal
        t(i,k,j,14,taun) = 0.06e-6                    ! Ca 
        t(i,k,j,15,taun) =1.0e-6                      !CaCO3
        t(i,k,j,16,taun) = tco2in(k)                  !TCO2 (UPDATED) 
        t(i,k,j,17,taun) = talkin(k)                  !TALK  (UPDATED)
          
	enddo
 
        pme_offset = 0.0
        
        year_no = 2 
        do year_n = 1,year_no
        days  = 365*5                                 ! change later1460!365!360 
        dt = 864
        deg2rad = 3.14/180.0
        loop_total = days*86400/dt
        pme_sum = 0.0
        pme_int = 0.0
         
!***********TO INTERPOLATE THE DATA IN PROGRAME***********
 
	
        open(2,file='INPUT/air_temp_data.txt')
        open(3,file='INPUT/uwnd_data.txt') 
        open(4,file='INPUT/vwnd_data.txt') 
        open(6,file='INPUT/sfx_in_data.txt')                  !!sfx and ppt data are in
                                                        !data_for_argo2009 inside ECMWF directory
 
        open(10,file='INPUT/ppt_data.txt')
        open(13,file='INPUT/specific_hum_data.txt')
        open(unit=9,
     &file=
     &'INPUT/wvel_data.bin'
     &,status='old',
     &access='direct',form='unformatted',recl=51)

        open(unit=15,file='INPUT/salinity_clim_data.txt')
        

        open(unit=7,file='INPUT/pCO2_in_data.txt') !NEWWW
        
        open(unit=8,file='INPUT/ontofact_data.txt') !NEWWW
        
        read(2,*)Tair1
        read(2,*)Tair2
        read(3,*)uwnd1
        read(3,*)uwnd2
        read(4,*)vwnd1
        read(4,*)vwnd2
!        read(5,*)hum1
!        read(5,*)hum2
        read(6,*)sfx_in1  
        read(6,*)sfx_in2 
        read(10,*)ppt1   
        read(10,*)ppt2   
!        read(11,*)evp1  
!        read(11,*)evp2 
        read(13,*)spe_hum1
        read(13,*)spe_hum2
        read(15,*)sal_clim1
        read(15,*)sal_clim2 

        read(7,*)co2_in1 !NEWWW 
        read(7,*)co2_in2 !NEWWW

        do kk=1,51
        read(8,*)Onto1(kk)
        read(8,*)Onto2(kk)
        enddo
        
! **********READING A BINARY FILE******************
        
        
        irec = 1

        do ii = 1,615  !122                                !jj=1 is jan 2 2009,jj=2 jan 5 2009,jj=122 31 dec 2009 !615 
        read(9,rec=irec) (wvel(ii,jj),jj=1,km)
        irec=irec+1
       

                if (ii.eq.1) then
                  do k=1,km
                  if (wvel(ii,k).lt. -999) wvel(ii,k) = 0.0
                  wvel1(k)=wvel(ii,k) 
                  enddo
                end if
                
                if (ii.eq.2) then 
                   do k=1,km
                   if (wvel(ii,k).lt. -999) wvel(ii,k) = 0.0
                   wvel2(k)=wvel(ii,k)
                   enddo
                end if
        enddo
        
 
        close(2)
        close(3)
        close(4)
!        close(5)
        close(6)
        close(10)
!        close(11)
        close(13)
        close(9)
        close(15) 
       
        close(7)!NEWWW

        lpd =1.0
        lsd =1.0      
        lkd = 1.0
        lmd = 1.0 
        lco2= 1.0 !NEWWW
        lOnt= 1.0 !NEWWW
        day_no =0.0
        do loop = 1,loop_total
        
        !DATA
        open(2,file='INPUT/air_temp_data.txt')
        open(3,file='INPUT/uwnd_data.txt')
        open(4,file='INPUT/vwnd_data.txt')
        open(6,file='INPUT/sfx_in_data.txt')
        open(10,file='INPUT/ppt_data.txt')
        open(13,file='INPUT/specific_hum_data.txt') 
        open(unit=9,
     &file=
     &'INPUT/wvel_data.bin'
     &,status='old',
     &access='direct',form='unformatted',recl=51)
        open(unit=15,file='INPUT/salinity_clim_data.txt')
        
        open(7,file='INPUT/pCO2_in_data.txt')!NEWWW
 
        
        call relax_nutrient(loop,no3_relax_obs,sioh4_relax_obs) 
         
        if(loop.eq.1) then
        Tair = Tair1
        uwnd = uwnd1
        vwnd = vwnd1
!        hum = hum1
        sfx_in = sfx_in1
        ppt = ppt1
!        evp = evp1
        spe_hum = spe_hum1
        sal_clim = sal_clim1
        co2_in = co2_in1 !NEWWW

        do k=1,km
        Onto(k) = Onto1(k)
        wvelocity(k)=wvel1(k)
        enddo
        
        else


        Tair= Tair + (( (Tair2 - Tair1) * 864.0) / (21600.0))!864 and 21600
        uwnd = uwnd + (((uwnd2 - uwnd1) *864.0) /(21600.0))!change in
        vwnd = vwnd + (((vwnd2 - vwnd1) *864.0) /(21600.0))!terms of dt.
!        hum = hum + (((hum2 - hum1) *864.0) /(21600.0))
        sfx_in = sfx_in + (((sfx_in2 - sfx_in1 ) *864.0) / (43200.0))
        ppt = ppt + (((ppt2 - ppt1 ) * 864.0) / (43200.0))
!        evp = evp + (((evp2 - evp1 ) * 864.0) / (21600.0))
        spe_hum = spe_hum + (((spe_hum2- spe_hum1)* 864.0)/(21600.0))
        !sal_clim = sal_clim + (((sal_clim2- sal_clim1)* 864.0)/(2628000.0)) 
        sal_clim = sal_clim + (((sal_clim2- sal_clim1)* 864.0)/(259200)) 
        
        co2_in = co2_in + (((co2_in2 - co2_in1 ) *864.0) /(86400.0))
        !NEWWW

  
        
        do k=1,km
        Onto(k)=Onto(k) + (((Onto2(k)- Onto1(k))*864.0)
     &     /(259200))
        wvelocity(k)=wvelocity(k)+(((wvel2(k)- wvel1(k))*864.0)
     &     /(259200))
        enddo
        
        end if
        
        close(2)
        close(3)
        close(4)
!        close(5)
        close(6)
        close(10)
!        close(11)
        close(13)
        close(9)
        close(15)
        close(7)!NEWWW
        
        

        if(mod(loop,25).eq.0) then
        lpd = lpd+1
        
        call read_variable (lpd, Tair1,
     &   Tair2,uwnd1,uwnd2,vwnd1,vwnd2,
     &   spe_hum1,spe_hum2)
      
        Tair = Tair1
        uwnd = uwnd1
        vwnd = vwnd1
!        hum = hum1
!        sfx_in = sfx_in1
!        ppt = ppt1
!        evp = evp1
        spe_hum = spe_hum1
        end if
        
        
        if(mod(loop,50).eq.0) then
        lkd = lkd+1
        
        

        
        call read_variable_settwo (lkd,
     &  sfx_in1,sfx_in2,ppt1,ppt2)

        sfx_in = sfx_in1
        ppt = ppt1
        
        end if

        if(mod(loop,300).eq.0) then 
        lsd = lsd+1
        
        call read_variable_wvelocity(lsd,wvel,wvel1,wvel2,km,
     &  wvelocity) 
        
        do k=1,km
        wvelocity(k)=wvel1(k)
        enddo 
        end if
       

!        if(mod(loop,3041).eq.0)then !3042 0r 3041 !CHECK MUSTTTTTTTTTTTTTTT
        if(mod(loop,300).eq.0)then !3042 0r 3041 !CHECK MUSTTTTTTTTTTTTTTT
        lmd = lmd + 1
       
        call read_variable_sal_clim (lmd,sal_clim1,sal_clim2)
        sal_clim = sal_clim1
        end if


        if(mod(loop,100).eq.0) then !NEWWW
        lco2 = lco2+1
         call read_variable_pCO2(lco2,co2_in1,co2_in2)
         co2_in = co2_in1
         end if

        if(mod(loop,100).eq.0) then !NEWWW
        lOnt = lOnt + 1
         call read_variable_Onto(lOnt,Onto1,Onto2)
         Onto(:)=Onto1(:)
        endif


        
! ALL ARGOS 13 and 16*********************************        
        if (mod(loop,146000).eq.0) then !146000, was 36500
        lpd = 1.0
        lsd = 1.0
        lkd = 1.0
        end if
       
!!!        if(mod(loop,36500).eq.0) then
!!!        lsd=1.0
!!!        end if
!*******************************************************************

 
! FOR ARGO 12 and ARGO 14 and ARGO15***********************************************

!        if (year_n .eq. 1) then
!       
!                if (mod(loop,146000).eq.0) then !146000, was 36500
!                lpd = 1.0
!               lkd = 1.0
!               end if
!                if(mod(loop,36500).eq.0) then
!                lsd=1.0
!               end if
!        else
!                if (mod(loop,146000).eq.0) then !146000, was 36500
!               lpd = 1.0
!                lsd = 1.0
!                lkd = 1.0
!                end if
!        end if
!**************************************************************************


!DELETE AND ADD LSD IN ABOVE LOOP
!        if(mod(loop,36500).eq.0) then
!        lsd=1.0
!        end if
!!DELETE
        if(mod(loop,146000).eq.0)then  !36500
        lmd=1.0
        endif
        
!**************************************************now************************
        do k=1,km
        pomin(k)  = t(i,k,j,10,taun)!dbt in this 1e6 needed or not
        opalin(k) = t(i,k,j,13,taun) !molN/L
        caco3in(k)= t(i,k,j,15,taun)
        no3in(k)  = t(i,k,j,8,taun)
        sioh4in(k)= t(i,k,j,12,taun)
        tempin(k) = t(i,k,j,1,taun) 
        saltin(k) = t(i,k,j,2,taun)
        uvelin(k) = u(i,k,j,taun)  
        vvelin(k) = v(i,k,j,taun)  
        enddo
!********************************************************************now****        
        
        POM_F_depth = 0.0  !downward flux of POM at a depth 
        Opal_F_depth = 0.0
        CaCO3_F_depth = 0.0

        ZLOnt(:) = 0.0
        do k = 1,km
        ZLOnt(k) = t(i,k,j,6,taun)
        enddo
        do k=1,km
	T0    =  t(i,k,j,1,taun)
        S0    =  t(i,k,j,2,taun) 
        UVEL0 =  u(i,k,j,taun)   
        VVEL0 =  v(i,k,j,taun)   
        PS    =  t(i,k,j,3,taun) !molN/L
        PL    =  t(i,k,j,4,taun)
        ZS    =  t(i,k,j,5,taun)
        ZL    =  t(i,k,j,6,taun)
        ZP    =  t(i,k,j,7,taun)
        NO3   =  t(i,k,j,8,taun)
        NH4   =  t(i,k,j,9,taun)
        POM   =  t(i,k,j,10,taun)
        DOM   =  t(i,k,j,11,taun)
        SIOH4 =  t(i,k,j,12,taun)
        Opal  =  t(i,k,j,13,taun)
        Ca    =  t(i,k,j,14,taun)
        CaCO3 = t(i,k,j,15,taun)
        TCO2  = t(i,k,j,16,taun)
        TALK  = t(i,k,j,17,taun)
        no3_obss=no3_relax_obs(k)*1e-6       !no3_obs(k)*1e-6
        sioh4_obss=sioh4_relax_obs(k)*1e-6   !sioh4_obs(k)*1e-6
        if (k .gt. 1 .and. k .lt. 30) then
        ZLSUM = -0.5*ZLOnt(k) + 0.5*ZLOnt(k+int(Onto(k)))
        else
        ZLSUM = 0.0
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!treatement for surface light
        Alpha1    = 0.025  !0.03   light dissipation co-ef       !0.035(from org. Nemuro)  to 2.5*0.04 
        Alpha2    = 0.025/1e-6 !0.018/1e-6  self shading co-ef    !0.0281/1e-6(original nemuro.) 
        
        rsum = 0.0
        day_night = sin((2.0*3.14159265/86400.0)*loop*dt)
        day_night = max(0.0d0 , day_night)
        sfx = sfx_in * 0.45
        
        do kk = 1,k
        rkappa = Alpha1 + Alpha2 * (t(i,kk,j,3,taun) + t(i,kk,j,4,taun))
        rsum = rsum +   rkappa * dz(kk)
        enddo

        rlight = sfx * exp(-1.0d0 * rsum)
        rlightk(k) = rlight
       
!*******************************************************************************        
        !to get profile of C:Chla ratio
        
        C_atm_mass = 12.0  
        Red_rat = 117.0/16.0 ! Redfld ratio
        C_chla_rat = 40.0 !C:Chla ratio
        
        Chloro = (Red_rat * (t(i,kk,j,3,taun) + t(i,kk,j,4,taun)) *
     &            C_atm_mass * 1e6 ) / (C_chla_rat)
                
       
        KWA = 0.028 ! attenuation coeff. of  water
        KCA = 0.058 ! attenuation coeff. of chlorophyll
        KDA = 0.008 ! attenuation coeff. of DOM
        
        KAA = KWA + (KCA*Chloro) + (KDA*DOM*1e6)  ! DOM in mmol N/m3
        
        mlight = sfx * exp(-1*KAA * dz(k)*k)
        
!*********************************************************************************8        
!*********************************************************************************************        
        rtemp = T0
        xconv = (1.0/3.6e5)
        ra = 0.337
        xkw_temp = (uwnd**2 + vwnd**2)*xconv*ra
        rsc = scco2 (rtemp)
        rkw = xkw_temp * (660/rsc)**0.5

        temp_in = T0
        salt_in = S0
        kw660   = rkw
        ppo     = 1.0 + (dz(k) - 5.0) * 0.1
        po4     = (NO3/16.0) * 1025.0  ! mole/m3
        dic_in  = TCO2 * 1025.0  ! mole/m3
        alk_in  = TALK * 1025.0  ! mole/m3
        si      = SIOH4* 1025.0  ! mole/m3
        rhoin   = 1025.0
        xco2    = co2_in  !NEWWW
        dz1     = 5.0



	call nemuro_3D (
     I  PS,PL,ZS,ZL,ZP,NO3,NH4,POM,DOM,SIOH4,Opal,Ca, CaCO3, TCO2, TALK,
     I  T0,S0,UVEL0,VVEL0, loop, dt, rlight,k,km,pomin,opalin,
     I  cain,caco3in,no3in,sioh4in,ppindx,
     I  wvelocity,tempin,saltin,uvelin,vvelin, t0_grad, s0_grad, 
     I  uvel0_grad, vvel0_grad,no3_grad,POM_F_depth,Opal_F_depth,
     I  CaCO3_F_depth,no3_obss,sioh4_obss,
     I  RnewL,RnewS,NO3_lim,NH4_lim,N_lim,Si_lim,diff_N_and_Si_lim,
     I  no3_relax_obs,sioh4_relax_obs,sfx,
     I  temp_in, salt_in, kw660, ppo, dic_in, alk_in, po4, si, rhoin, xco2,
     I  dz1 , co2ex1,ph1,pCO2surf,ff,pco2Z,phZ, tco2in_rstr,
     I  talkin_rstr, rlxindx,ZLSUM,Onto)

        t(i,k,j,1,taun) = T0 
        t(i,k,j,2,taun) = S0 
        u(i,k,j,taun)   = UVEL0 
        v(i,k,j,taun)   = VVEL0 
        t(i,k,j,3,taun) = PS !molN/L
        t(i,k,j,4,taun) = PL
        t(i,k,j,5,taun) = ZS
        t(i,k,j,6,taun) = ZL
        t(i,k,j,7,taun) = ZP
        t(i,k,j,8,taun) = NO3
        t(i,k,j,9,taun) = NH4
        t(i,k,j,10,taun)= POM
        t(i,k,j,11,taun)= DOM
        t(i,k,j,12,taun)= SIOH4
        t(i,k,j,13,taun)= Opal
        t(i,k,j,14,taun)= Ca
        t(i,k,j,15,taun)= CaCO3
        t(i,k,j,16,taun)= TCO2
        t(i,k,j,17,taun)= TALK
        
        RnewS_d(k)             = RnewS
        RnewL_d(k)             = RnewL
        NO3_lim_d(k)           = NO3_lim
        NH4_lim_d(k)           = NH4_lim
        N_lim_d(k)             = N_lim
        Si_lim_d(k)            = Si_lim
        diff_N_and_Si_lim_d(k) = diff_N_and_Si_lim

	enddo
        

		
	do k=1,km
	tempin(k) = t(i,k,j,1,taun)
	saltin(k) = t(i,k,j,2,taun)
	uvelin(k) = u(i,k,j,taun)
	vvelin(k) = v(i,k,j,taun)
	psin(k) = t(i,k,j,3,taun)*1e6 !micromolN/L
	plin(k) = t(i,k,j,4,taun)*1e6
	zsin(k) = t(i,k,j,5,taun)*1e6
	zlin(k) = t(i,k,j,6,taun)*1e6
	zpin(k) = t(i,k,j,7,taun)*1e6
	no3in(k) = t(i,k,j,8,taun)*1e6
	nh4in(k) = t(i,k,j,9,taun)*1e6
	pomin(k) = t(i,k,j,10,taun)*1e6
	domin(k) = t(i,k,j,11,taun)*1e6
	sioh4in(k) = t(i,k,j,12,taun)*1e6
	opalin(k) = t(i,k,j,13,taun)*1e6
	cain(k) = t(i,k,j,14,taun)*1e6
	caco3in(k) = t(i,k,j,15,taun)*1e6
	tco2in(k) = t(i,k,j,16,taun)*1e6
	talkin(k) = t(i,k,j,17,taun)*1e6
        
        !treatment of light taken as equation from the book "Small scale
!processes in Geophysical fluid flows" page no :124
        
	enddo
        
        call light_treatment_by_band (sfx_in,km,sswin,psin,plin)
        SST_in =  tempin(1)
        
c        call stresstc (uwnd, vwnd, Tair,SST_in,
c     &  tauxin, tauyin, vel)

        
                        
	dt_rg = dt
!======================================================================================
        cp_a  = 1004.0                  ! (J/kg degree C), spec.h eat of air
        cp_w  = 4200.0                  ! J/kg degree C), spec.heat water
        rho_a = 1.175                   !or 1.22 ! kg/m3, air density
        RO    = 1025.0                  !ocean density. kg/m3 ! 
!        cs    = 1.0e-3                 ! empirical coefficient
!        cl    = 1.5e-3  
!        cd    = 1.0e-3               ! empirical coefficient
        Lv    = 2.44e6                  ! J/kg, latnt heat of evaporatn.
        Pa    = 1000.0                  ! atm. pressre. dobt in unt milibar??
        emiss = 0.97                    ! emissivity of sea surface
        sigma_const = 5.670367e-8       ! Stefan-Boltzman constant. check it.
        c1    = 0.55                    ! cloud coefficient.
        C     = 0.1                     ! Fractional cloudiness., it can have values between 0 and 1

        
        
        !DRAG VALUES FROM SMITH (1988), JGR, TABLE-2. FOR T > 20 DEGREES
        ! AND T < 20
        !ref :/home/user/Desktop/PHd/my_work/coefficints to improve SST.odt
 
        vel = sqrt(uwnd**2+vwnd**2)
!	write (*,*) uwnd,vwnd,vel

        if(SST_in.ge.20.0) then
        if((vel.ge.0.0).and.(vel.lt.2.0))cs =2.07e-3            !cs -empirical coeff.
        if((vel.ge.2.0).and.(vel.lt.5.0))cs =1.48e-3
        if((vel.ge.5.0).and.(vel.lt.10.0))cs =1.20e-3
        if((vel.ge.10.0).and.(vel.lt.15.0))cs =1.10e-3
        if((vel.ge.15.0).and.(vel.lt.20.0))cs =1.06e-3
        if((vel.ge.20.0).and.(vel.lt.25.0))cs =1.04e-3
        cl=cs*1.20                                              !cl - empirical coefficnt.
        end if

        if(SST_in.lt.20.0) then
        if((vel.ge.0.0).and.(vel.lt.2.0))cs =1.96e-3
        if((vel.ge.2.0).and.(vel.lt.5.0))cs =1.42e-3
        if((vel.ge.5.0).and.(vel.lt.10.0))cs =1.17e-3
        if((vel.ge.10.0).and.(vel.lt.15.0))cs =1.08e-3
        if((vel.ge.15.0).and.(vel.lt.20.0))cs =1.05e-3
        if((vel.ge.20.0).and.(vel.lt.25.0))cs =1.03e-3
        cl=cs*1.20
        end if

        if (SST_in.ge.20) then
        if((vel.ge.0.0).and.(vel.lt.2.0))cd = 1.68e-3
        if((vel.ge.2.0).and.(vel.lt.5.0))cd = 1.45e-3
        if((vel.ge.5.0).and.(vel.lt.10.0))cd = 1.57e-3
        if ((vel.ge.10.0).and.(vel.lt.15.0))cd = 1.74e-3
        if ((vel.ge.15.0).and.(vel.lt.20.0))cd = 1.93e-3
        if ((vel.ge.20.0).and.(vel.lt.25.0))cd = 2.14e-3
        endif

        if (SST_in.lt.20) then
        if ((vel.ge.0.0).and.(vel.lt.2.0)) cd = 1.61e-3
        if ((vel.ge.2.0).and.(vel.lt.5.0)) cd = 1.40e-3
        if ((vel.ge.5.0).and.(vel.lt.10.0)) cd = 1.52e-3
        if ((vel.ge.10.0).and.(vel.lt.15.0)) cd = 1.71e-3
        if ((vel.ge.15.0).and.(vel.lt.20.0)) cd = 1.91e-3
        if ((vel.ge.20.0).and.(vel.lt.25.0)) cd = 2.12e-3
        endif

        tauxin= rho_a*cd*uwnd*abs(uwnd)
        tauyin= rho_a*cd*vwnd*abs(vwnd)
        

!sensible heat flux         
        qs = rho_a*cp_a*cs*vel*(Tair-SST_in)

!latent heat flux

!        dq=(hum/100.0)*qsat(Tair, Pa) - 0.98*qsat(SST_in, Pa) ! when
!        using relative hum.
        dq = spe_hum - 0.98*qsat(SST_in, Pa) ! when using spec. hum
        ql= rho_a * cl * Lv * vel * dq

!	print *,spe_hum
!************************************************************************        
        !NOTE
!        if dq is +ve ql is becomes +ve, in that case evp is becomes
!        +ve. that is it turn in to ppt. so make evp =0.0 for +ve ql. 
        !calculate evp frm ql / (1000 * Lv). 
!MODEL STARTS FROM 5.30 mrng. bec data starts frm 3.00 UTC. 
!***********************************************************************        
        if(ql.gt.0)then
        evp = 0.0 
        else
        evp = ql/(Lv*1024.5)
        end if
!KEEP   
        pme = ((ppt) + evp)*(-1.0) - pme_offset
        pme_sum = pme_sum + pme
        
        if(mod(loop,3000).eq.0)then
        pme_offset1 = pme_sum/(3000)
        end if

!        change 3500 to 3000
        
        if(mod(loop,3000).eq.0)then
        pme_offset = 0.0
        else 
        pme_offset = pme_offset1 
        end if 
      
        
! KEEP

!        pme = ppt+evp - pme_offset
!
!        pme_sum = pme_sum + ppt+evp!pme 
!        
!        if(mod(loop,36500).eq.0)then
!        pme_offset1 = pme_sum/(36500)
!        end if
!        
!        if(year_n.le.1)then
!        pme_offset = 0.0
!        else
!        pme_offset = pme_offset1
!        end if

! KEEP


!**********************************DELETE*****************************
!        pme_offset = 9.531E-08

!       pme =   ppt - (-1*evp)+pme_offset
!        pme_offset = 0.0
!************************************************************************* 
 
        salt = pme * 35.0 * (1.0/(dz(1)))   
                                                !salt should be less than 0.
!*******************************************************************************        
                                  
        
                                   
!*********long wave radiation frm stephen-boltzmanns law***********
                Tair_Ke = Tair+273.15                               ! temp in kelvin
                                                                    !refe. kantha clayson page:425
        alpha = 0.045
        lw_sst = -1* emiss*sigma_const * ((SST_in+273.15)**4)  

        lw_air = sigma_const*(Tair_Ke**4)*(1-(0.261*exp(-7.77e-4*
     &   (273.15-Tair_Ke)**2)))
        
        lw = (lw_sst) + (1-alpha)*(lw_air)                          !in buk 1-alpha is
                                                                    !0.045  and -ve sign is not there in the eqation.
!***************************************************************************      

        emmis=0.9 !0.97
                                        ! KANTHA AND CLAYSON, Small scale processes, Page-425
        lw=-1.0*emmis*sigma_const*(SST_in+273.15)**4 + (1)*
     &      sigma_const*((Tair+273.15)**4)* (1.0-0.261*exp((-7.77e-4)*
     &      (273.15 - (Tair+273.15))**2))
                                        ! net outward flux = sensible heat + latent heat - longwave radiation
        qout=qs+ql+lw

        qnet =(qout)/(RO*cp_w*(dz(1)))   ! qnet in dT/dt unit.
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
!	do kkk_light = 2,km
!	sswin(kkk_light) = rlightk(kkk_light)
!	enddo
	call MT1 (tempin, saltin, uvelin, vvelin, 
     &  psin, plin, zsin,zlin,zpin,no3in,nh4in,pomin,domin,sioh4in,
     &  opalin,cain,caco3in,tco2in,talkin,
     &  tauxin,tauyin,
     &  qnet, pme, rLATin,dt_rg,SMout,SHout,sswin,salt,cp_w,sal_clim,nonlocal*mod(loop,2)) !,rlightk) !,salt,cp_w) 
 
        taup = taun   ! ONLY TEMPORORY FOR OUTPUTS, MUST CHANGE LATER
	do k=1,km
	t(i,k,j,1,taun) = tempin(k) 
	t(i,k,j,2,taun) = saltin(k) 
	u(i,k,j,taun) = uvelin(k) 
	v(i,k,j,taun) = vvelin(k)
	t(i,k,j,3,taun) = psin(k) *1e-6 !molN/L 
	t(i,k,j,4,taun) = plin(k) *1e-6
	t(i,k,j,5,taun) = zsin(k) *1e-6
	t(i,k,j,6,taun) = zlin(k) *1e-6
	t(i,k,j,7,taun) = zpin(k) *1e-6
	t(i,k,j,8,taun) = no3in(k) *1e-6
	t(i,k,j,9,taun) = nh4in(k) *1e-6
	t(i,k,j,10,taun) = pomin(k) *1e-6
	t(i,k,j,11,taun) = domin(k) *1e-6
	t(i,k,j,12,taun) = sioh4in(k) *1e-6
	t(i,k,j,13,taun) = opalin(k) *1e-6
	t(i,k,j,14,taun) = cain(k) * 1e-6
	t(i,k,j,15,taun) = caco3in(k) *1e-6
	t(i,k,j,16,taun) = tco2in(k) *1e-6 
	t(i,k,j,17,taun) = talkin(k) *1e-6
	enddo
	t(i,1,j,16,taun) = tco2in(1) *1e-6   + pme * dt_rg * 1.9/(5.0*1024.5)
	t(i,1,j,17,taun) = talkin(1) *1e-6   + pme * dt_rg * 2.2/(5.0*1024.5)
        
        
!**********************************************************       
       !TO GET DAILY AVERAGE VALUES************ 
       
        do k=1,km
        temp_var(k)  = tempin(k)
        salt_var(k)  = saltin(k)
        ps_var(k)    = psin(k) * 1e-6 !molN/L 
        pl_var(k)    = plin(k) * 1e-6
        zs_var(k)    = zsin(k) * 1e-6
        zl_var(k)    = zlin(k) * 1e-6
        zp_var(k)    = zpin(k) * 1e-6
        no3_var(k)   = no3in(k) * 1e-6
        nh4_var(k)   = nh4in(k) * 1e-6
        pom_var(k)   = pomin(k) * 1e-6
        dom_var(k)   = domin(k) * 1e-6
        sioh4_var(k) = sioh4in(k) * 1e-6
        opal_var(k)  = opalin(k) * 1e-6
        ca_var(k)    = cain(k) * 1e-6
        caco3_var(k) = caco3in(k)*1e-6
        tco2_var(k)  = tco2in(k) * 1e-6
        talk_var(k)  = talkin(k) * 1e-6
        uvel_var(k)  = uvelin(k)
        vvel_var(k)  = vvelin(k)
        ExcT0_var(k) = t0_grad(k)
        ExcS0_var(k) = s0_grad(k)
        !ExcUVEL0_var(k)=uvel0_grad(k)
        !ExcVVEL0_var(k)=vvel0_grad(k)
        ExcNO3_var(k) = no3_grad(k)
        wvelocity_var(k) = wvelocity(k)
        rlight_var(k)  = rlightk(k)
        RnewS_var(k) = RnewS_d(k)
        RnewL_var(k) = RnewL_d(k)
        NO3_lim_var(k) =NO3_lim_d(k) 
        NH4_lim_var(k) =  NH4_lim_d(k)
        N_lim_var(k) = N_lim_d(k)
        Si_lim_var(k) =Si_lim_d(k)
        diff_N_and_Si_lim_var(k)= diff_N_and_Si_lim_d(k)
        PSCChla_var(k)    =  PSCChla(k)
        PLCChla_var(k)    =  PLCChla(k)
        ph_var(k) = phZ(k)
        co2ex_var(k) = co2ex1! NEW ADDITION
        pCO2surf_var(k) = pCO2Z(k) ! NEW ADDITION
        xco2_var(k)     = xco2!NEWWW
        enddo

        
        
!        if(year_n.eq.7) then
        call daily_ave(temp_var,salt_var,ps_var,pl_var,zs_var,
     &  zl_var,zp_var,no3_var,nh4_var,pom_var,dom_var,
     &  sioh4_var,opal_var,ca_var,caco3_var,tco2_var,talk_var,
     &  rlight_var,SHout,SMout,km,loop,day_no,qnet,sfx_in,ppt,evp,pme,
     &  uvel_var,vvel_var,ExcT0_var,ExcS0_var,
     &  ExcNO3_var,wvelocity_var,Tair,uwnd,vwnd,spe_hum,qs,ql,lw,sswin,
     &  RnewS_var,RnewL_var,NO3_lim_var,NH4_lim_var,N_lim_var,
     &  Si_lim_var,diff_N_and_Si_lim_var,PSCChla_var,
     &  PLCChla_var,ph_var,co2ex_var,pCO2surf_var,xco2_var)!NEW ADDITION
!        end if
! hum is going to daily, but spe_hum is using        
	
c	write(*,*) loop, 'loop', ps_var

        end do 

!***********************************************************************************
        
	loop_ind = 0
	if (mod (loop,100) .eq. 0 ) then
         loop_ind = loop_ind + 1
        endif
        enddo

        end
        
        subroutine read_variable (lpd,Tair1, Tair2,uwnd1,uwnd2,
     & vwnd1,vwnd2,spe_hum1,spe_hum2)
        real(4)  Tair1, Tair2, uwnd1, uwnd2, vwnd1, vwnd2, 
     &  spe_hum1,spe_hum2
        integer lpd 
        
        !DATA
        open(2,file='INPUT/air_temp_data.txt')
        open(3,file='INPUT/uwnd_data.txt')
        open(4,file='INPUT/vwnd_data.txt')
!        open(5,file='relative_hum_data.txt')
!        open(6,file='sfx_in_data.txt')
!        open(10,file='ppt_data.txt')
!        open(11,file='evp_data.txt')
        open(13,file='INPUT/specific_hum_data.txt')


        do i= 1,lpd-1
        read(2,*)Tair1
        read(3,*)uwnd1
        read(4,*)vwnd1
!        read(5,*)hum1
!        read(6,*)sfx_in1
!        read(10,*)ppt1
!        read(11,*)evp1
        read(13,*)spe_hum1
        enddo
        

        read(2,*)Tair1
        read(2,*)Tair2
        read(3,*)uwnd1
        read(3,*)uwnd2
        read(4,*)vwnd1
        read(4,*)vwnd2
!        read(5,*)hum1
!        read(5,*)hum2
!        read(6,*)sfx_in1
!        read(6,*)sfx_in2
!        read(10,*),ppt1
!        read(10,*)ppt2
!        read(11,*)evp1
!        read(11,*)evp2
        read(13,*)spe_hum1
        read(13,*)spe_hum2
!*************************************        
                
        
        Tair=Tair1
        uwnd=uwnd1
        vwnd=vwnd1
!        hum=hum1
!        sfx_in=sfx_in1
!        ppt = ppt1
!        evp = evp1
        spe_hum = spe_hum1
          

        close(2)
        close(3)
        close(4)
!        close(5)
!        close(6)
!        close(10)
!        close(11)
        close(13)

        end subroutine read_variable  
        
         subroutine read_variable_settwo (lkd,
     & sfx_in1,sfx_in2,ppt1,ppt2)
        real(4) sfx_in1,sfx_in2,ppt1,ppt2
        integer lkd
        
        !DATA
        open(6,file='INPUT/sfx_in_data.txt')
        open(10,file='INPUT/ppt_data.txt')
        
         do i= 1,lkd-1
        read(6,*)sfx_in1
        read(10,*)ppt1
        enddo
        
        read(6,*)sfx_in1
        read(6,*)sfx_in2
        read(10,*)ppt1
        read(10,*)ppt2

        sfx_in=sfx_in1
        ppt = ppt1

        close(6)
        close(10)
                   
        end subroutine read_variable_settwo


        subroutine read_variable_wvelocity(lsd,wvel,wvel1,wvel2,km,
     &  wvelocity)
        real(4) wvel(615,51),wvel1(51),wvel2(51),wvelocity(51)
        integer lsd
       
        
        do j=1,km
        wvel1(j) = wvel(lsd+0,j)
        if (wvel1(j).lt. -999) wvel1(j) = 0.0 
        
        wvel2(j) = wvel(lsd+1,j)
        if (wvel2(j).lt. -999) wvel2(j) = 0.0 
        enddo
       
        do k=1,km
        wvelocity(k)=wvel1(k)
        enddo
        
        end subroutine read_variable_wvelocity
        
        subroutine read_variable_sal_clim(lmd,sal_clim1,sal_clim2)
        real(4) sal_clim1,sal_clim2
        integer lmd
        
         !DATA
        open(15,file='INPUT/salinity_clim_data.txt')

         do i= 1,lmd-1
        read(15,*)sal_clim1
        enddo

        read(15,*)sal_clim1
        read(15,*)sal_clim2

        sal_clim=sal_clim1

        close(15)

        end subroutine read_variable_sal_clim

        subroutine read_variable_Onto(lpd, Onto1, Onto2)
        real Onto1(51), Onto2(51)
        integer lpd

        open(11,file='INPUT/ontofact_data.txt')

        do i = 1,lpd-1
        do kk=1,51
        read(11,*)Onto1(kk)
        enddo
        enddo

        do kk=1,51
        read(11,*)Onto1(kk)
        enddo
        do kk=1,51
        read(11,*)Onto2(kk)
        enddo
        close(11)

        end subroutine

        function qsat(Ta,Pa)
        real Ta
        real Pa,ew,qsat
        ew = 6.1121*(1.0007+3.46e-6*Pa)*exp((17.502*Ta)/(240.97+Ta))!in mb
!doubt, check it for more information to the woodshol routine.
        qsat = 0.62197*(ew/(Pa-0.378*ew))  !% mb -> kg/kg 
 
        end function qsat

        subroutine read_variable_pCO2 (lpd, co2_in1, co2_in2)
        real co2_in1, co2_in2

        open (10, file='INPUT/pCO2_in_data.txt')

        do i= 1, lpd-1
        read (10,*) co2_in1
        enddo

        read (10,*) co2_in1
        read (10,*) co2_in2

        close (11)
        end subroutine read_variable_pCO2

