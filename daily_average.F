        subroutine daily_ave(temp_var,salt_var,ps_var,pl_var,zs_var,
     &  zl_var,zp_var,no3_var,nh4_var,pom_var,dom_var,
     &  sioh4_var,opal_var,ca_var,caco3_var,tco2_var,talk_var,
     &  rlight_var,SHout,SMout,km,loop,day_no,qnet,sfx_in,ppt,evp,pme,
     &  uvel_var,vvel_var,ExcT0_var,ExcS0_var,
     &  ExcNO3_var,wvelocity_var,Tair,uwnd,vwnd,spe_hum,qs,ql,lw,
     &  sswin,RnewS_var,RnewL_var,NO3_lim_var,NH4_lim_var,N_lim_var,
     &  Si_lim_var,diff_N_and_Si_lim_var,PSCChla_var,PLCChla_var,
     &  ph_var,co2ex_var,pCO2surf_var,xco2_var)! NEW ADDITION !NEWWW
        !add loop_int and remove PSCChla_var,PLCChla_var) while running along ship track 
!     &  Si_lim_var,diff_N_and_Si_lim_var,loop_int)!,PSCChla_var,PLCChla_var) add loop_int and remove PSCChla_var,PLCChla_var) while running along ship track 

!add ,PSCChla_var,PLCChla_var) and remove loop_int while running along
!argo tracks
                                          
       

!NOTE FOR RUN along ship include receive loop_int in above routine to
!control the output.

 
        integer loop_int  !define only along ship track
        real(8) PS_factor,PL_factor,PS_Chla(51),PL_Chla(51)
        real Tair,uwnd,vwnd,hum,qs,ql,lw 
        real ExcT0_var(51),ExcT0_c(100,51),ExcT0_sum,ExcT0_total!full
        real ExcT0_mean(51)!full
        real ExcS0_var(51),ExcS0_c(100,51),ExcS0_sum,ExcS0_total!full
        real ExcS0_mean(51)!full
        real ExcUVEL0_var(51),ExcUVEL0_c(100,51),ExcUVEL0_sum
        real ExcUVEL0_total,ExcUVEL0_mean(51)!full
        real ExcVVEL0_var(51),ExcVVEL0_c(100,51),ExcVVEL0_sum
        real ExcVVEL0_total,ExcVVEL0_mean(51)!full
        real ExcNO3_var(51),ExcNO3_c(100,51),ExcNO3_sum,ExcNO3_total!fu
        real ExcNO3_mean(51)!full


!***********************ENABLE ALONG ARGO*******************************        
!        real(8)PSCChla_var(51),PLCChla_var(51),PSCChla_total     
!        real(8)PLCChla_total,PSCChla_mean(51),PLCChla_mean(51)
!**************************************************************************
        
       !WHY RnewS and RnewL average is coming wrong, is it bec of
       !declaration as real(8) 
        real(8)RnewS_var(51),RnewL_var(51),RnewL_c(100,51)
        real(8)RnewS_c(100,51),RnewS_sum,RnewL_sum
        real(8)RnewS_total,RnewL_total,RnewL_mean(51),RnewS_mean(51)

        !real RnewS_var(51),RnewL_var(51),RnewL_c(100,51)
        !real RnewS_c(100,51),RnewS_sum,RnewL_sum
        !real RnewS_total,RnewL_total,RnewL_mean(51),RnewS_mean(51)



        real temp_var(51),salt_var(51),ps_var(51),pl_var(51)
        real zs_var(51),zl_var(51),zp_var(51),no3_var(51)
        real nh4_var(51),pom_var(51)
        real dom_var(51),sioh4_var(51),opal_var(51),ca_var(51)
        real caco3_var(51)
        real tco2_var(51),talk_var(51)
        real  uvel_var(51),vvel_var(51)! 15/02/2017-uvel_var onwards
        real rlight_var(51),SHout(51),SMout(51)
        real wvelocity_var(51), ph_var(51)
        real co2ex_var(51),pCO2surf_var(51)
        real xco2_var(51)!NEWWW
        

                            !now

        real NO3_lim_var(51),NH4_lim_var(51),N_lim_var(51)
        real Si_lim_var(51),diff_N_and_Si_lim_var(51)
        real NO3_lim_c(100,51),NH4_lim_c(100,51),N_lim_c(100,51)
        real Si_lim_c(100,51),diff_N_and_Si_lim_c(100,51)
        real NO3_lim_sum,NH4_lim_sum,N_lim_sum,Si_lim_sum
        real diff_N_and_Si_lim_sum,NO3_lim_total,NH4_lim_total
        real N_lim_total,Si_lim_total,diff_N_and_Si_lim_total   
        real NO3_lim_mean(51),NH4_lim_mean(51), N_lim_mean(51)
        real Si_lim_mean(51),diff_N_and_Si_lim_mean(51)
 
        real temp_c(100,51),salt_c(100,51),ps_c(100,51)
        real pl_c(100,51),zs_c(100,51),zl_c(100,51),zp_c(100,51)
        real no3_c(100,51),nh4_c(100,51),pom_c(100,51),dom_c(100,51)
        real sioh4_c(100,51),opal_c(100,51),ca_c(100,51),caco3_c(100,51)
        real tco2_c(100,51),talk_c(100,51),rlight_c(100,51)
        real SHout_c(100,51),SMout_c(100,51)
        real uvel_c(100,51),vvel_c(100,51)!15/02/2017-full line
        real wvelocity_c(100,51)
        real ph_c(100,51)!NEW ADDITION
        real co2ex_c(100,51)
        real pCO2surf_c(100,51)
        real xco2_c(100,51) !NEWWW


!ENABLE ALONG
!ARGO*******************************************************
!        real(8)PSCChla_c(100,51),PLCChla_c(100,51)
!**********************************************************************
        real temp_sum,salt_sum,ps_sum,pl_sum,zs_sum,zl_sum,zp_sum
        real no3_sum,nh4_sum,pom_sum,dom_sum,sioh4_sum,opal_sum
        real ca_sum,caco3_sum,tco2_sum,talk_sum,rlight_sum
        real SHout_sum,SMout_sum
        real wvelocity_sum
        real ph_sum !NEW ADDITION
        real co2ex_sum
        real pCO2surf_sum
        real xco2_sum!NEWWW

!ENABLE ALONG ARGO**************************************
!        real(8)PSCChla_sum,PLCChla_sum 
!*******************************************************

        real uvel_sum,vvel_sum !15/02/2017-full line
        real temp_total,salt_total,ps_total,pl_total,zs_total
        real zl_total,zp_total,no3_total,nh4_total,pom_total,dom_total
        real sioh4_total,opal_total,ca_total,caco3_total,tco2_total
        real talk_total,rlight_total,SHout_total,SMout_total
        real uvel_total,vvel_total !15/02/2017-full line
        real wvelocity_total
        real ph_total!NEW ADDITION
        real co2ex_total
        real pCO2surf_total
        real xco2_total !NEWWW


        real temp_mean(51),salt_mean(51),ps_mean(51),pl_mean(51)
        real zs_mean(51),zl_mean(51),zp_mean(51),no3_mean(51)
        real nh4_mean(51),pom_mean(51),dom_mean(51),sioh4_mean(51)
        real opal_mean(51),ca_mean(51),caco3_mean(51),tco2_mean(51)
        real talk_mean(51),rlight_mean(51),SHout_mean(51)
        real uvel_mean(51),vvel_mean(51)!15/02/2017-full line
        real ph_mean(51)!NEW ADDITION
        real co2ex_mean(51)
        real pCO2surf_mean(51)
        real xco2_mean(51)!NEWWW

        real SMout_mean(51),qnet,qnet_c(100),qnet_sum,qnet_total
        real light_mean(51),sfx_in,sfx_in_c(100),sfx_in_sum
        real ppt_c(100),evp_c(100),ppt_sum,evp_sum,ppt_total
        real wvelocity_mean(51)
        real Tair_c(100),uwnd_c(100),vwnd_c(100),hum_c(100),qs_c(100)
        real ql_c(100),lw_c(100)
        real Tair_sum,uwnd_sum,vwnd_sum,hum_sum,qs_sum,ql_sum,lw_sum
        real Tair_total,uwnd_total,vwnd_total,hum_total,qs_total
        real ql_total,lw_total
        
        real sswin(51),sswin_c(100,51),sswin_sum,sswin_total
        real sswin_mean(51) 
        real evp_total,ppt_mean,evp_mean
        real Tair_mean,uwnd_mean,vwnd_mean,hum_mean,qs_mean,ql_mean
        real lw_mean
        real qnet_mean,qnet_mean_org,sfx_in_total,sfx_in_mean
        real ppt,evp
        real pme,pme_c(100),pme_sum,pme_total,pme_mean
        real(8) chla_PS(51),chla_PL(51),fact, CChla_rat_PS, CChla_rat_PL 
        integer km,loop,d,day_no
        
        
        !do i=1,51 
        !print*,PLCChla_var(i),i,loop !,PLCChla_var,k,loop 
        !end do

!        do k=1,51
!        if((loop.gt.100).and.(loop.le.200).and.(k.eq.1))print*
!     &   ,RnewL_var(k),loop,k
!        end do

        !print*,loop_int 
       d = loop - (day_no*100)

!	write (*,*) "loop:",loop,"day_no:",day_no,"d:",d
                temp_c(d,1:51)   =  temp_var(1:51)
                salt_c(d,1:51)   =  salt_var(1:51)
                ps_c(d,1:51)     =  ps_var(1:51)
                pl_c(d,1:51)     =  pl_var(1:51)
                zs_c(d,1:51)     =  zs_var(1:51)
                zl_c(d,1:51)     =  zl_var(1:51)
                zp_c(d,1:51)     =  zp_var(1:51)
                no3_c(d,1:51)    =  no3_var(1:51)
                nh4_c(d,1:51)    =  nh4_var(1:51)
                pom_c(d,1:51)    =  pom_var(1:51)
                dom_c(d,1:51)    =  dom_var(1:51)
                sioh4_c(d,1:51)  =  sioh4_var(1:51)
                opal_c(d,1:51)   =  opal_var(1:51)
                ca_c(d,1:51)     =  ca_var(1:51)
                caco3_c(d,1:51)  =  caco3_var(1:51)
                tco2_c(d,1:51)   =  tco2_var(1:51)
                talk_c(d,1:51)   =  talk_var(1:51)
                rlight_c(d,1:51) =  rlight_var(1:51)
                SHout_c(d,1:51)  =  SHout(1:51)
                SMout_c(d,1:51)  =  SMout(1:51)
                uvel_c(d,1:51)   =  uvel_var(1:51)!15/02/2017
                vvel_c(d,1:51)   =  vvel_var(1:51)!15/02/2017
                ph_c(d,1:51)     =  ph_var(1:51)!NEW ADDITION
                co2ex_c(d,1:51)  =  co2ex_var(1:51)
                pCO2surf_c(d,1:51) =  pCO2surf_var(1:51)
                xco2_c(d,1:51)   =  xco2_var(1:51)!NEWWW

                qnet_c(d)        =  qnet
                sfx_in_c(d)      = sfx_in
                ppt_c(d)         = ppt
                evp_c(d)         = evp
                pme_c(d)         = pme
                Tair_c(d)        = Tair

                uwnd_c(d)        = uwnd
                vwnd_c(d)        = vwnd
                hum_c(d)         = spe_hum
                qs_c(d)          = qs
                ql_c(d)          = ql
                lw_c(d)          = lw
                ExcT0_c(d,1:51)  = ExcT0_var(1:51)
                ExcS0_c(d,1:51)  = ExcS0_var(1:51)
!                ExcUVEL0_c(d,1:51)=ExcUVEL0_var(1:51)
!                ExcVVEL0_c(d,1:51)=ExcVVEL0_var(1:51)
                ExcNO3_c(d,1:51)  =ExcNO3_var(1:51)
                wvelocity_c(d,1:51)   =  wvelocity_var(1:51)
                sswin_c(d,1:51)   = sswin(1:51)
                RnewL_c(d,1:51)   = RnewL_var(1:51)
                RnewS_c(d,1:51)   = RnewS_var(1:51)
                NO3_lim_c(d,1:51) = NO3_lim_var(1:51)      
                NH4_lim_c(d,1:51) = NH4_lim_var(1:51)
                N_lim_c(d,1:51)   = N_lim_var(1:51)
                Si_lim_c(d,1:51)  = Si_lim_var(1:51)
                diff_N_and_Si_lim_c(d,1:51)=diff_N_and_Si_lim_var(1:51)       
 
!ENABLE ALONG ARGO*********************************************               
!                PSCChla_c(d,1:51) = PSCChla_var(1:51)
!                PLCChla_c(d,1:51) = PLCChla_var(1:51)
!**************************************************************                 

        if(mod(loop,100).eq.0)then

         do k=1,km
          temp_sum     =  0.0
          salt_sum     =  0.0
          ps_sum       =  0.0
          pl_sum       =  0.0
          zs_sum       =  0.0
          zl_sum       =  0.0
          zp_sum       =  0.0
          no3_sum      =  0.0
          nh4_sum      =  0.0
          pom_sum      =  0.0
          dom_sum      =  0.0
          sioh4_sum    =  0.0
          opal_sum     =  0.0
          ca_sum       =  0.0
          caco3_sum    =  0.0
          tco2_sum     =  0.0
          talk_sum     =  0.0
          rlight_sum  =  0.0
          SHout_sum    =  0.0
          SMout_sum    =  0.0
          qnet_sum     =  0.0
          sfx_in_sum   =  0.0
          ppt_sum      =  0.0
          evp_sum      =  0.0
          pme_sum      =  0.0
          Tair_sum     =  0.0
          uwnd_sum     =  0.0
          vwnd_sum     =  0.0
          hum_sum      =  0.0
          qs_sum       =  0.0
          ql_sum       =  0.0
          lw_sum       =  0.0
          ph_sum       =  0.0!NEW ADDITION
          co2ex_sum    =  0.0
          pCO2surf_sum =  0.0
          xco2_sum     =  0.0!NEWWW


          uvel_sum     =  0.0 !15/02/2017
          vvel_sum     =  0.0 !15/02/2017
          ExcT0_sum    =  0.0
          ExcS0_sum    =  0.0
!          ExcUVEL0_sum =  0.0
!          ExcVVEL0_sum =  0.0
          ExcNO3_sum   =  0.0
          wvelocity_sum = 0.0
          sswin_sum     = 0.0
          RnewS_sum     = 0.0
          RnewL_sum     = 0.0
          NO3_lim_sum   = 0.0
          NH4_lim_sum   = 0.0
          N_lim_sum     = 0.0
          Si_lim_sum    = 0.0
          diff_N_and_Si_lim_sum = 0.0

!ENABLE ALONG ARGO*******************************************        
!          PSCChla_sum =0.0
!          PLCChla_sum =0.0
!************************************************************

        do ii=1,100

           temp_total     = temp_sum + temp_c(ii,k)
           temp_sum       = temp_total
                
           salt_total     = salt_sum + salt_c(ii,k)
           salt_sum       = salt_total
        
           ps_total       = ps_sum + ps_c(ii,k)
           ps_sum         = ps_total

           pl_total       = pl_sum + pl_c(ii,k)
           pl_sum         = pl_total
                
           zs_total       = zs_sum + zs_c(ii,k)
           zs_sum         = zs_total
        
           zl_total       = zl_sum + zl_c(ii,k)
           zl_sum         = zl_total
        
           zp_total       = zp_sum + zp_c(ii,k)
           zp_sum         = zp_total
        
           no3_total      = no3_sum + no3_c(ii,k)
           no3_sum        = no3_total
        
           nh4_total      = nh4_sum + nh4_c(ii,k)
           nh4_sum        = nh4_total
        
           pom_total      = pom_sum + pom_c(ii,k)
           pom_sum        = pom_total
        
           dom_total      = dom_sum + dom_c(ii,k)
           dom_sum        = dom_total
                
           sioh4_total    = sioh4_sum + sioh4_c(ii,k)
           sioh4_sum      = sioh4_total
        
           opal_total     = opal_sum + opal_c(ii,k)
           opal_sum       = opal_total
        
           ca_total       = ca_sum + ca_c(ii,k)
           ca_sum         = ca_total
        
           caco3_total    = caco3_sum + caco3_c(ii,k)
           caco3_sum      = caco3_total
        
           tco2_total     = tco2_sum + tco2_c(ii,k)
           tco2_sum       = tco2_total

           talk_total     = talk_sum + talk_c(ii,k)
           talk_sum       = talk_total
        
           rlight_total  = rlight_sum + rlight_c(ii,k)
           rlight_sum    = rlight_total
!           print*,loop,k,rlight_sum 
        
           SHout_total    = SHout_sum + SHout_c(ii,k)
           SHout_sum      = SHout_total
        
           SMout_total    = SMout_sum + SMout_c(ii,k)
           SMout_sum      = SMout_total
        
           uvel_total     = uvel_sum + uvel_c(ii,k)!15/02/2017
           uvel_sum       = uvel_total             !15/02/2017
           
           vvel_total     = vvel_sum + vvel_c(ii,k)!15/02/2017
           vvel_sum       = vvel_total             !15/02/2017
        
           qnet_total     = qnet_sum + qnet_c(ii)
           qnet_sum       = qnet_total
        
           sfx_in_total   = sfx_in_sum + sfx_in_c(ii)
           sfx_in_sum     = sfx_in_total
                
           ppt_total      = ppt_sum + ppt_c(ii)
           ppt_sum        = ppt_total
          
           evp_total      = evp_sum + evp_c(ii)
           evp_sum        = evp_total
        
           pme_total      = pme_sum + pme_c(ii)
           pme_sum        = pme_total
        
           Tair_total      = Tair_sum + Tair_c(ii)
           Tair_sum        = Tair_total
        
           uwnd_total      = uwnd_sum + uwnd_c(ii)
           uwnd_sum        = uwnd_total
        
           vwnd_total      = vwnd_sum + vwnd_c(ii)
           vwnd_sum        = vwnd_total

           hum_total      = hum_sum + hum_c(ii)
           hum_sum        = hum_total

        
           qs_total      = qs_sum + qs_c(ii)
           qs_sum        = qs_total

        
           ql_total      = ql_sum + ql_c(ii)
           ql_sum        = ql_total

        
           lw_total      = lw_sum + lw_c(ii)
           lw_sum        = lw_total

        
           ExcT0_total    = ExcT0_sum + ExcT0_c(ii,k)
           ExcT0_sum      = ExcT0_total
        
           ExcS0_total    = ExcS0_sum + ExcS0_c(ii,k)
           ExcS0_sum      = ExcS0_total


           ph_total       = ph_sum + ph_c(ii,k)!NEW ADDITION
           ph_sum         = ph_total


           co2ex_total    = co2ex_sum + co2ex_c(ii,k)!NEW ADDITION
           co2ex_sum         = co2ex_total


           pCO2surf_total = pCO2surf_sum + pCO2surf_c(ii,k)!NEW ADDITION
           pCO2surf_sum   = pCO2surf_total


           
           xco2_total    = xco2_sum + xco2_c(ii,k)!NEWWW
           xco2_sum         = xco2_total




!           ExcUVEL0_total = ExcUVEL0_sum + ExcUVEL0_c(ii,k)
!           ExcUVEL0_sum   = ExcUVEL0_total

!           ExcVVEL0_total = ExcVVEL0_sum + ExcVVEL0_c(ii,k)
!           ExcVVEL0_sum   = ExcVVEL0_total
          
           ExcNO3_total   = ExcNO3_sum + ExcNO3_c(ii,k)
           ExcNO3_sum     = ExcNO3_total
        
           wvelocity_total     = wvelocity_sum + wvelocity_c(ii,k)
           wvelocity_sum       = wvelocity_total
        
           sswin_total     = sswin_sum + sswin_c(ii,k)
           sswin_sum       = sswin_total
        
           RnewL_total     = RnewL_sum + RnewL_c(ii,k)
           RnewL_sum       = RnewL_total
        
           RnewS_total     = RnewS_sum + RnewS_c(ii,k)
           RnewS_sum       = RnewS_total
        
           NO3_lim_total  = NO3_lim_sum + NO3_lim_c(ii,k)
           NO3_lim_sum    = NO3_lim_total
        
           NH4_lim_total  = NH4_lim_sum + NH4_lim_c(ii,k)
           NH4_lim_sum    = NH4_lim_total
                 
           N_lim_total  = N_lim_sum + N_lim_c(ii,k)
           N_lim_sum    = N_lim_total
        
           Si_lim_total  = Si_lim_sum + Si_lim_c(ii,k)
           Si_lim_sum    = Si_lim_total
        
           diff_N_and_Si_lim_total  = diff_N_and_Si_lim_sum + 
     &  diff_N_and_Si_lim_c(ii,k)
           diff_N_and_Si_lim_sum    = diff_N_and_Si_lim_total

!ENABLE ALONG ARGO***********************************************
!        PSCChla_total  = PSCChla_sum + PSCChla_c (ii,k)
!        PSCChla_sum    = PSCChla_total
        
!        PLCChla_total  = PLCChla_sum + PLCChla_c (ii,k)
!        PLCChla_sum    = PLCChla_total
!***************************************************************        
        
        end do
        
        temp_mean(k)       =  temp_total/100.0 !100 no of loops for 1 day
        salt_mean(k)       =  salt_total/100.0 !give it in terms of d if 
        ps_mean(k)         =  ps_total/100.0   !possible.
        pl_mean(k)         =  pl_total/100.0
        zs_mean(k)         =  zs_total/100.0
        zl_mean(k)         =  zl_total/100.0
        zp_mean(k)         =  zp_total/100.0
        no3_mean(k)        =  no3_total/100.0
        nh4_mean(k)        =  nh4_total/100.0
        pom_mean(k)        =  pom_total/100.0
        dom_mean(k)        =  dom_total/100.0
        sioh4_mean(k)      =  sioh4_total/100.0
        opal_mean(k)       =  opal_total/100.0
        ca_mean(k)         =  ca_total/100.0
        caco3_mean(k)      =  caco3_total/100.0
        tco2_mean(k)       =  tco2_total/100.0
        talk_mean(k)       =  talk_total/100.0
        rlight_mean(k)    =  rlight_total/100.0
        SHout_mean(k)      =  SHout_total/100.0
        SMout_mean(k)      =  SMout_total/100.0
        qnet_mean          =  qnet_total/100.0
        sfx_in_mean        = sfx_in_total /100.0
        ppt_mean           = ppt_total/100.0
        evp_mean           = evp_total/100.0
        pme_mean           = pme_total/100.0
        Tair_mean          = Tair_total/100.0
        uwnd_mean          = uwnd_total/100.0
        vwnd_mean          = vwnd_total/100.0
        hum_mean           = hum_total/100.0
        qs_mean            = qs_total/100.0
        ql_mean            = ql_total/100.0
        lw_mean            = lw_total/100.0

        ph_mean(k)         =  ph_total/100.0!NEW ADDITION
        co2ex_mean(k)      =  co2ex_total/100.0
        pCO2surf_mean(k)   =  pCO2surf_total/100.0
        xco2_mean(k)      =  xco2_total/100.0



        uvel_mean(k)       = uvel_total/100.0!15/02/2017
        vvel_mean(k)       = vvel_total/100.0!15/02/2017 
        ExcT0_mean(k)      = ExcT0_total/100.0
        ExcS0_mean(k)      = ExcS0_total/100.0
!        ExcUVEL0_mean(k)   = ExcUVEL0_total/100.0
!        ExcVVEL0_mean(k)   = ExcVVEL0_total/100.0
        ExcNO3_mean(k)     = ExcNO3_total/100.0
        wvelocity_mean(k)       =  wvelocity_total/100.0
        sswin_mean(k)   = sswin_total/100.0
      
        RnewL_mean(k)       =  RnewL_total/100.0
        RnewS_mean(k)       = RnewS_total/100.0!100 no of loops for 1 day
       
        NO3_lim_mean(k)     = NO3_lim_total/100.0
        NH4_lim_mean(k)     = NH4_lim_total/100.0
        N_lim_mean(k)       = N_lim_total/100.0
        Si_lim_mean(k)     = Si_lim_total/100.0
        diff_N_and_Si_lim_mean(k)     =  diff_N_and_Si_lim_total/100.0

!*****************************************ENABLE ALONG
!ARGO***************        
!        PSCChla_mean(k)    =PSCChla_total/100.0
!        PLCChla_mean(k)    =PLCChla_total/100.0
        
        
!        if(PSCChla_mean(k).le.0)then
!        PSCChla_mean(k)=1.0
!        else
!        PSCChla_mean(k)=PSCChla_mean(k)
!        end if
        
!        if(PLCChla_mean(k) .le.0)then
!        PLCChla_mean(k) =1.0
!        else
!        PLCChla_mean(k) =PLCChla_mean(k) 
!        end if
!        PS_factor = (117.0*12.0*1e6)/(16.0*PSCChla_mean(k))
!        PS_Chla(k) = PS_mean(k)*PS_factor
        
!        PL_factor = (117.0*12.0*1e6)/(16.0*PLCChla_mean(k))
!        PL_Chla(k) =PL_mean(k)*PL_factor
        
!*****************************************************************************        
        
        
!JUST FOR TEST RUN DELETE THESE TERMS LATER        
        qnet_mean_org = qnet_mean * 1000.0 * 4200.0* 5 !(it is in w/m2
!not in unit of temp)
!!************************************************************************

! ENABLE ALONG ARGO********************************
        
!         fact = (117.0 *12.0*1e6)/16.0
!          
!        if(k.le.8)then
!        CChla_rat_PS = 60.0
!        CChla_rat_PL = 40.0
!        else 
!        CChla_rat_PS = 40.0
!        CChla_rat_PL = 20.0
!        end if

!        chla_PS(k)=ps_mean(k)* fact/50.0
!        chla_PL(k)=pl_mean(k)* fact/30.0
        
!***************************************************        
        
        
        
!        if((mod(loop,100).eq.0).and.(loop_int.gt.149100))then !DISABLE
!ABOVE IF LINE FOR ARGO RUNS
!!only enable run along ship track 
	write(*,*) temp_mean(k),salt_mean(k),ps_mean(k),pl_mean(k),
     & zs_mean(k),zl_mean(k),zp_mean(k),no3_mean(k),nh4_mean(k),
     & pom_mean(k),sioh4_mean(k), ca_mean(k) ,
     &tco2_mean(k),caco3_mean(k),talk_mean(k),dom_mean(k),
     &ph_mean(k),co2ex_mean(k),pCO2surf_mean(k),xco2_mean(k) ! NEWADDITION
        
!	write(*,*) NO3_lim_mean(k),NH4_lim_mean(k),Si_lim_mean(k),
!     & N_lim_mean(k),diff_N_and_Si_lim_mean(k)

!         write(*,*) temp_mean(k),salt_mean(k),ps_mean(k),pl_mean(k),
!     & zs_mean(k),zl_mean(k),zp_mean(k),no3_mean(k),nh4_mean(k),
!     & hum_mean,sfx_in_mean,ppt_mean,pme_mean,
!     & Tair_mean,uwnd_mean,vwnd_mean,wvelocity_mean(k)

!        end if

! remember here qnet_mean is actualy qnet_mean*rho*cp*depth
! (qnet_mean*1000*4200*5)
        end do

        day_no=day_no+1
        end if

        return
        end
