        subroutine daily_ave(temp_var,salt_var,ps_var,pl_var,zs_var,
     &  zl_var,zp_var,no3_var,nh4_var,pom_var,dom_var,
     &  sioh4_var,opal_var,ca_var,caco3_var,tco2_var,talk_var,
     &  rlight_var,SHout,SMout,km,loop,day_no,qnet,sfx_in,ppt,evp,pme,
     &  uvel_var,vvel_var,ExcT0_var,ExcS0_var,
     &  ExcNO3_var,wvelocity_var,Tair,uwnd,vwnd,spe_hum,qs,ql,lw,
     &  sswin,RnewS_var,RnewL_var,NO3_lim_var,NH4_lim_var,N_lim_var,
     &  Si_lim_var,diff_N_and_Si_lim_var,PSCChla_var,PLCChla_var,
     &  ph_var,co2ex_var,pCO2surf_var,xco2_var,dz,dt)

C     Local variables and daily-history arrays
        integer loop_int
        real(8) PS_factor,PL_factor,PS_Chla(51),PL_Chla(51)
        real Tair,uwnd,vwnd,hum,qs,ql,lw
        real ExcT0_var(51),ExcT0_sum,ExcT0_total
        real ExcT0_mean(51)
        real ExcS0_var(51),ExcS0_sum,ExcS0_total
        real ExcS0_mean(51)
        real ExcUVEL0_var(51),ExcUVEL0_sum
        real ExcUVEL0_total,ExcUVEL0_mean(51)
        real ExcVVEL0_var(51),ExcVVEL0_sum
        real ExcVVEL0_total,ExcVVEL0_mean(51)
        real ExcNO3_var(51),ExcNO3_sum,ExcNO3_total
        real ExcNO3_mean(51)

        real(8)RnewS_var(51),RnewL_var(51)
        real(8)RnewS_sum,RnewL_sum
        real(8)RnewS_total,RnewL_total,RnewL_mean(51),RnewS_mean(51)

        real temp_var(51),salt_var(51),ps_var(51),pl_var(51)
        real zs_var(51),zl_var(51),zp_var(51),no3_var(51)
        real nh4_var(51),pom_var(51)
        real dom_var(51),sioh4_var(51),opal_var(51),ca_var(51)
        real caco3_var(51)
        real tco2_var(51),talk_var(51)
        real  uvel_var(51),vvel_var(51)
        real rlight_var(51),SHout(51),SMout(51)
        real wvelocity_var(51), ph_var(51)
        real co2ex_var(51),pCO2surf_var(51)
        real xco2_var(51)

        real NO3_lim_var(51),NH4_lim_var(51),N_lim_var(51)
        real Si_lim_var(51),diff_N_and_Si_lim_var(51)
        real NO3_lim_sum,NH4_lim_sum,N_lim_sum
        real Si_lim_sum,diff_N_and_Si_lim_sum
        real NO3_lim_total,NH4_lim_total
        real N_lim_total,Si_lim_total,diff_N_and_Si_lim_total
        real NO3_lim_mean(51),NH4_lim_mean(51), N_lim_mean(51)
        real Si_lim_mean(51),diff_N_and_Si_lim_mean(51)

        real temp_sum,salt_sum,ps_sum,pl_sum,zs_sum,zl_sum,zp_sum
        real no3_sum,nh4_sum,pom_sum,dom_sum,sioh4_sum,opal_sum
        real ca_sum,caco3_sum,tco2_sum,talk_sum,rlight_sum
        real SHout_sum,SMout_sum
        real wvelocity_sum
        real ph_sum
        real co2ex_sum
        real pCO2surf_sum
        real xco2_sum

        real uvel_sum,vvel_sum
        real temp_total,salt_total,ps_total,pl_total,zs_total
        real zl_total,zp_total,no3_total,nh4_total,pom_total,dom_total
        real sioh4_total,opal_total,ca_total,caco3_total,tco2_total
        real talk_total,rlight_total,SHout_total,SMout_total
        real uvel_total,vvel_total
        real wvelocity_total
        real ph_total
        real co2ex_total
        real pCO2surf_total
        real xco2_total

        real temp_mean(51),salt_mean(51),ps_mean(51),pl_mean(51)
        real zs_mean(51),zl_mean(51),zp_mean(51),no3_mean(51)
        real nh4_mean(51),pom_mean(51),dom_mean(51),sioh4_mean(51)
        real opal_mean(51),ca_mean(51),caco3_mean(51),tco2_mean(51)
        real talk_mean(51),rlight_mean(51),SHout_mean(51)
        real uvel_mean(51),vvel_mean(51)
        real ph_mean(51)
        real co2ex_mean(51)
        real pCO2surf_mean(51)
        real xco2_mean(51)

        real SMout_mean(51),qnet,qnet_sum,qnet_total
        real light_mean(51),sfx_in,sfx_in_sum
        real ppt_sum,evp_sum,ppt_total
        real wvelocity_mean(51)

        real Tair_sum,uwnd_sum,vwnd_sum,hum_sum,qs_sum,ql_sum,lw_sum
        real Tair_total,uwnd_total,vwnd_total,hum_total,qs_total
        real ql_total,lw_total

        real sswin(51),sswin_sum,sswin_total
        real sswin_mean(51)
        real evp_total,ppt_mean,evp_mean
        real Tair_mean,uwnd_mean,vwnd_mean,hum_mean,qs_mean,ql_mean
        real lw_mean
        real qnet_mean,qnet_mean_org,sfx_in_total,sfx_in_mean
        real ppt,evp
        real pme,pme_sum,pme_total,pme_mean
        real(8) chla_PS(51),chla_PL(51),fact, CChla_rat_PS, CChla_rat_PL
        integer km,loop,d,day_no,nsteps_day,ii,dt
        real dz
        integer, save :: nsteps_allocated = 0

        real, allocatable, save :: ExcT0_c(:,:),ExcS0_c(:,:)
        real, allocatable, save :: ExcUVEL0_c(:,:),ExcVVEL0_c(:,:)
        real, allocatable, save :: ExcNO3_c(:,:)
        real(8), allocatable, save :: RnewS_c(:,:),RnewL_c(:,:)
        real, allocatable, save :: NO3_lim_c(:,:),NH4_lim_c(:,:)
        real, allocatable, save :: N_lim_c(:,:),Si_lim_c(:,:)
        real, allocatable, save :: diff_N_and_Si_lim_c(:,:)
        real, allocatable, save :: temp_c(:,:),salt_c(:,:),ps_c(:,:)
        real, allocatable, save :: pl_c(:,:),zs_c(:,:),zl_c(:,:)
        real, allocatable, save :: zp_c(:,:),no3_c(:,:),nh4_c(:,:)
        real, allocatable, save :: pom_c(:,:),dom_c(:,:),sioh4_c(:,:)
        real, allocatable, save :: opal_c(:,:),ca_c(:,:),caco3_c(:,:)
        real, allocatable, save :: tco2_c(:,:),talk_c(:,:),rlight_c(:,:)
        real, allocatable, save :: SHout_c(:,:),SMout_c(:,:)
        real, allocatable, save :: uvel_c(:,:),vvel_c(:,:)
        real, allocatable, save :: wvelocity_c(:,:),ph_c(:,:)
        real, allocatable, save :: co2ex_c(:,:),pCO2surf_c(:,:)
        real, allocatable, save :: xco2_c(:,:),sswin_c(:,:)
        real, allocatable, save :: qnet_c(:),sfx_in_c(:)
        real, allocatable, save :: ppt_c(:),evp_c(:),pme_c(:)
        real, allocatable, save :: Tair_c(:),uwnd_c(:),vwnd_c(:)
        real, allocatable, save :: hum_c(:),qs_c(:),ql_c(:),lw_c(:)

C     Validate resolution settings
        if (dt .le. 0) then
          write(*,*) 'ERROR: dt must be positive (seconds).'
          stop 1
        endif
        if (dz .le. 0.0) then
          write(*,*) 'ERROR: dz must be positive (metres).'
          stop 1
        endif
        nsteps_day = nint(86400.0/real(dt))
        if (nsteps_day .lt. 1 .or.
     &      abs(real(nsteps_day)*real(dt)-86400.0) .gt. 0.01) then
          write(*,*) 'ERROR: dt must divide 86400 seconds exactly.'
          write(*,*) 'dt = ',dt
          stop 1
        endif

C     Allocate history arrays for one model day
        if (nsteps_allocated .eq. 0) then
          allocate(ExcT0_c(nsteps_day,51),ExcS0_c(nsteps_day,51),
     &             ExcUVEL0_c(nsteps_day,51),
     &             ExcVVEL0_c(nsteps_day,51),
     &             ExcNO3_c(nsteps_day,51),
     &             RnewS_c(nsteps_day,51),RnewL_c(nsteps_day,51),
     &             NO3_lim_c(nsteps_day,51),
     &             NH4_lim_c(nsteps_day,51),N_lim_c(nsteps_day,51),
     &             Si_lim_c(nsteps_day,51),
     &             diff_N_and_Si_lim_c(nsteps_day,51),
     &             temp_c(nsteps_day,51),salt_c(nsteps_day,51),
     &             ps_c(nsteps_day,51),pl_c(nsteps_day,51),
     &             zs_c(nsteps_day,51),zl_c(nsteps_day,51),
     &             zp_c(nsteps_day,51),no3_c(nsteps_day,51),
     &             nh4_c(nsteps_day,51),pom_c(nsteps_day,51),
     &             dom_c(nsteps_day,51),sioh4_c(nsteps_day,51),
     &             opal_c(nsteps_day,51),ca_c(nsteps_day,51),
     &             caco3_c(nsteps_day,51),tco2_c(nsteps_day,51),
     &             talk_c(nsteps_day,51),rlight_c(nsteps_day,51),
     &             SHout_c(nsteps_day,51),SMout_c(nsteps_day,51),
     &             uvel_c(nsteps_day,51),vvel_c(nsteps_day,51),
     &             wvelocity_c(nsteps_day,51),ph_c(nsteps_day,51),
     &             co2ex_c(nsteps_day,51),
     &             pCO2surf_c(nsteps_day,51),
     &             xco2_c(nsteps_day,51),sswin_c(nsteps_day,51),
     &             qnet_c(nsteps_day),sfx_in_c(nsteps_day),
     &             ppt_c(nsteps_day),evp_c(nsteps_day),
     &             pme_c(nsteps_day),Tair_c(nsteps_day),
     &             uwnd_c(nsteps_day),vwnd_c(nsteps_day),
     &             hum_c(nsteps_day),qs_c(nsteps_day),
     &             ql_c(nsteps_day),lw_c(nsteps_day))
          nsteps_allocated = nsteps_day
        else if (nsteps_allocated .ne. nsteps_day) then
          write(*,*) 'ERROR: dt changed during daily_ave run.'
          stop 1
        endif

C     Store current timestep values
        d = loop - (day_no*nsteps_day)
        if (d .lt. 1 .or. d .gt. nsteps_day) then
          write(*,*) 'ERROR: daily_ave index outside day:',d,
     &               nsteps_day,loop,day_no
          stop 1
        endif

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
                uvel_c(d,1:51)   =  uvel_var(1:51)
                vvel_c(d,1:51)   =  vvel_var(1:51)
                ph_c(d,1:51)     =  ph_var(1:51)
                co2ex_c(d,1:51)  =  co2ex_var(1:51)
                pCO2surf_c(d,1:51) =  pCO2surf_var(1:51)
                xco2_c(d,1:51)   =  xco2_var(1:51)

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

C     Calculate daily means after a complete day
        if(mod(loop,nsteps_day).eq.0)then

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
          ph_sum       =  0.0
          co2ex_sum    =  0.0
          pCO2surf_sum =  0.0
          xco2_sum     =  0.0

          uvel_sum     =  0.0
          vvel_sum     =  0.0
          ExcT0_sum    =  0.0
          ExcS0_sum    =  0.0
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

C     Accumulate values across the day
        do ii=1,nsteps_day

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

           SHout_total    = SHout_sum + SHout_c(ii,k)
           SHout_sum      = SHout_total

           SMout_total    = SMout_sum + SMout_c(ii,k)
           SMout_sum      = SMout_total

           uvel_total     = uvel_sum + uvel_c(ii,k)
           uvel_sum       = uvel_total

           vvel_total     = vvel_sum + vvel_c(ii,k)
           vvel_sum       = vvel_total

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

           ph_total       = ph_sum + ph_c(ii,k)
           ph_sum         = ph_total

           co2ex_total    = co2ex_sum + co2ex_c(ii,k)
           co2ex_sum         = co2ex_total

           pCO2surf_total = pCO2surf_sum + pCO2surf_c(ii,k)
           pCO2surf_sum   = pCO2surf_total

           xco2_total    = xco2_sum + xco2_c(ii,k)
           xco2_sum         = xco2_total

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

        end do

C     Convert daily totals to means
        temp_mean(k)       =  temp_total/real(nsteps_day)
        salt_mean(k)       =  salt_total/real(nsteps_day)
        ps_mean(k)         =  ps_total/real(nsteps_day)
        pl_mean(k)         =  pl_total/real(nsteps_day)
        zs_mean(k)         =  zs_total/real(nsteps_day)
        zl_mean(k)         =  zl_total/real(nsteps_day)
        zp_mean(k)         =  zp_total/real(nsteps_day)
        no3_mean(k)        =  no3_total/real(nsteps_day)
        nh4_mean(k)        =  nh4_total/real(nsteps_day)
        pom_mean(k)        =  pom_total/real(nsteps_day)
        dom_mean(k)        =  dom_total/real(nsteps_day)
        sioh4_mean(k)      =  sioh4_total/real(nsteps_day)
        opal_mean(k)       =  opal_total/real(nsteps_day)
        ca_mean(k)         =  ca_total/real(nsteps_day)
        caco3_mean(k)      =  caco3_total/real(nsteps_day)
        tco2_mean(k)       =  tco2_total/real(nsteps_day)
        talk_mean(k)       =  talk_total/real(nsteps_day)
        rlight_mean(k)    =  rlight_total/real(nsteps_day)
        SHout_mean(k)      =  SHout_total/real(nsteps_day)
        SMout_mean(k)      =  SMout_total/real(nsteps_day)
        qnet_mean          =  qnet_total/real(nsteps_day)
        sfx_in_mean        = sfx_in_total /real(nsteps_day)
        ppt_mean           = ppt_total/real(nsteps_day)
        evp_mean           = evp_total/real(nsteps_day)
        pme_mean           = pme_total/real(nsteps_day)
        Tair_mean          = Tair_total/real(nsteps_day)
        uwnd_mean          = uwnd_total/real(nsteps_day)
        vwnd_mean          = vwnd_total/real(nsteps_day)
        hum_mean           = hum_total/real(nsteps_day)
        qs_mean            = qs_total/real(nsteps_day)
        ql_mean            = ql_total/real(nsteps_day)
        lw_mean            = lw_total/real(nsteps_day)

        ph_mean(k)         =  ph_total/real(nsteps_day)
        co2ex_mean(k)      =  co2ex_total/real(nsteps_day)
        pCO2surf_mean(k)   =  pCO2surf_total/real(nsteps_day)
        xco2_mean(k)      =  xco2_total/real(nsteps_day)

        uvel_mean(k)       = uvel_total/real(nsteps_day)
        vvel_mean(k)       = vvel_total/real(nsteps_day)
        ExcT0_mean(k)      = ExcT0_total/real(nsteps_day)
        ExcS0_mean(k)      = ExcS0_total/real(nsteps_day)
        ExcNO3_mean(k)     = ExcNO3_total/real(nsteps_day)
        wvelocity_mean(k)       =  wvelocity_total/real(nsteps_day)
        sswin_mean(k)   = sswin_total/real(nsteps_day)

        RnewL_mean(k)       =  RnewL_total/real(nsteps_day)
        RnewS_mean(k)       = RnewS_total/real(nsteps_day)

        NO3_lim_mean(k)     = NO3_lim_total/real(nsteps_day)
        NH4_lim_mean(k)     = NH4_lim_total/real(nsteps_day)
        N_lim_mean(k)       = N_lim_total/real(nsteps_day)
        Si_lim_mean(k)     = Si_lim_total/real(nsteps_day)
        diff_N_and_Si_lim_mean(k) = diff_N_and_Si_lim_total/
     &      real(nsteps_day)

C     Calculate heat flux in the original output units
        qnet_mean_org = qnet_mean * 1000.0 * 4200.0 * dz

C     Write daily mean output
        write(*,*) temp_mean(k),salt_mean(k),ps_mean(k),pl_mean(k),
     & zs_mean(k),zl_mean(k),zp_mean(k),no3_mean(k),nh4_mean(k),
     & pom_mean(k),sioh4_mean(k), ca_mean(k) ,
     &tco2_mean(k),caco3_mean(k),talk_mean(k),dom_mean(k),
     &ph_mean(k),co2ex_mean(k),pCO2surf_mean(k),xco2_mean(k)

        end do

C     Advance the day counter and finish
        day_no=day_no+1
        end if

        return
        end
