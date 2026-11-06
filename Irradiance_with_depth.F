!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! treatement for surface light
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        subroutine light_treatment_by_band (sfx_in,km,sswin,psin,plin)
      
        !NOTE ****************************************
        !inputs ==== sfx_in,km,psin,plin
        !output ==== sswin
        !*********************************************

        real(4)rAlpha1,rAlpha2,rL1,rL2,rchla(51),rsum,sfx_in,sfx
        real(4)rinfraRed
        real(4)rvisible
        real(4)sswin(51)
        real(4)rvisible_extinct,dz
        real(4)psin(51),plin(51)
        real(4)rchla_pass
        integer kk,km
        
        dz = 5.0
!        psin(k)=psin(k)*1e-6
!        plin(k)=plin(k)*1e-6

        atm_wt_N_in_g_per_mol = 14.01 !g/mol

        rchla(1) = (((psin(1)+plin(1)))*atm_wt_N_in_g_per_mol)
!chla should in terms of mg/m3. with out /1e6 rchla will be in mg/m3.

! NOTE:   13-Dec-2022:  ps[d=2,k=@ave]*1e6*14.01,
! mole        gram     10^3 
! ----    x  ------ x  ----  = mg/m3
! m3          mole     10^3 
! NOTE:   13-Dec-2022:  ps[d=2,k=@ave]*1e6*14.01,


        do kk = 2,km
        rAlpha1    = 0.805
        rAlpha2    = 0.195
        rL1        = 0.043   ! 0.55 is taken from Figure
!2.3.3,
!/home/user/Desktop/PHd/my_work/penetrative_solar_radiation_treatment/Irradiance
!with depth.pdf
         rL2        = 0.42  !* sqrt(0.55)
        
        
!        rsum_ps = 0.0
!        rsum_pl = 0.0
!        do k=1,kk
!        rsum_ps = rsum_ps + psin(k) 
!        rsum_pl = rsum_pl + plin(k) 
!        enddo
!        rsum_ps = rsum_ps/k
!        rsum_pl = rsum_ps/k
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
        rchla(kk)   = (((psin(kk)+plin(kk)))  *
     &        (atm_wt_N_in_g_per_mol))  !0.02 !(ps+pl) in micro_mol/L



!        rchla(kk)   = (rsum_ps+rsum_pl)  *
!     &        (atm_wt_N_in_g_per_mol) !0.02 !(ps+pl) in micro_mol/L

converted to equvalent amt of N. in mg/L  : ref
!home/user/Desktop/PHd/my_work/unit_conversions/onvert_micromol_per_L_to_mg_per_L.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        
        sfx = sfx_in 

        rinfraRed = 0.0 !sfx *0.45 * ( rAlpha1 * exp(-1.0d0 * (kk-1) * 5.0 
!     &* rL1) + rAlpha2 * exp(-1.0d0 * (kk-1) * 5.0 * rL2) )                                                            
                                                           !check the
                                                           !brackets
                                                           !properly
!!!note: IR will absorbs with in the 5-10 m surface ocean, so dont allow it to penetrate deep in the ocean. which will warm the subsurface and producing instability, mixing increasing. so made it 0.0.
        
        !rchla_pass = max(0.02, rchla(kk)) !check it max function here.
        !rchla_pass =  max(0.02, rchla(kk)) !0.02 !check it max function here.
        rchla_pass =  max(0.02, rchla(1)) !0.02 !check it max function here.
        
        
        rvisible  = sfx *0.55 * rvisible_extinct (rchla_pass, kk-1)
        sswin(kk) = rinfraRed + rvisible
        enddo
        sswin(1)  = sfx
        
        
        end subroutine light_treatment_by_band
! treatement for surface light
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 	function rvisible_extinct(rchla, kk)


c
c table of coefficients from Morel and Amtoine (1994)
c J. Physical Oceanography, Vol. 24, 1652-1664
c to find the attenuation of PAR from chlorophyll-a
c concentrations.
c

     	 real(4) rVA (2,6), rZA (2,6)
         real(4) rvisible_extinct

        rVA(1,1) =   0.353
        rVA(1,2) =  -0.047
        rVA(1,3) =  0.083
        rVA(1,4) =  0.047
        rVA(1,5) =  -0.011
        rVA(1,6) = -0.009

        rVA(2,1) = 0.647
        rVA(2,2) = 0.047
        rVA(2,3) = -0.083
        rVA(2,4) = -0.047
        rVA(2,5) =  0.011
        rVA(2,6) =  0.009

        rZA(1,1) =  1.662
        rZA(1,2) = -0.605
        rZA(1,3) =  0.128
        rZA(1,4) = -0.033
        rZA(1,5) = -0.051
        rZA(1,6) = -0.004

        rZA(2,1) =  8.541
        rZA(2,2) = -8.924
        rZA(2,3) =  4.020
        rZA(2,4) = -0.077
        rZA(2,5) = -0.536
        rZA(2,6) =  0.055

!note: took the coeff. frm the 2nd table in Morel and Amtoine (1994)
!c J. Physical Oceanography, Vol. 24, 1652-1664, this values r for a
!profile of chla.

         chl_X = log10(rchla)
         chl_X2= chl_X**2
         chl_X3= chl_X**3
         chl_X4= chl_X**4
         chl_X5= chl_X**5

         rtempV1 = rVA(1,1) + rVA(1,2)*chl_X + rVA(1,3)*chl_X2
     &                      + rVA(1,4)*chl_X3+ rVA(1,5)*chl_X4
     &                      + rVA(1,6)*chl_X5

         rtempV2 = rVA(2,1) + rVA(2,2)*chl_X + rVA(2,3)*chl_X2
     &                      + rVA(2,4)*chl_X3+ rVA(2,5)*chl_X4
     &                      + rVA(2,6)*chl_X5
        

         rtempZ1 = rZA(1,1) + rZA(1,2)*chl_X + rZA(1,3)*chl_X2
     &                      + rZA(1,4)*chl_X3+ rZA(1,5)*chl_X4
     &                      + rZA(1,6)*chl_X5

         rtempZ2 = rZA(2,1) + rZA(2,2)*chl_X + rZA(2,3)*chl_X2
     &                      + rZA(2,4)*chl_X3+ rZA(2,5)*chl_X4
     &                      + rZA(2,6)*chl_X5
        

         rtemp1 = rtempV1*exp(-1.0d0*kk*5.0/rtempZ1)
         rtemp2 = rtempV2*exp(-1.0d0*kk*5.0/rtempZ2)
!	print *,rtemp1,rtemp2
	rvisible_extinct = rtemp1 + rtemp2 
                                           
        
        end function rvisible_extinct
