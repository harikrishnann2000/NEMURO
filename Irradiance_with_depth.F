        subroutine light_treatment_by_band (sfx_in,km,sswin,psin,
     &                                      plin,dz_input)

C--- Variable declarations ---------------------------------------------
        real(4)rAlpha1,rAlpha2,rL1,rL2,rchla(51),rsum,sfx_in,sfx
        real(4)rinfraRed
        real(4)rvisible
        real(4)sswin(51)
        real(4)rvisible_extinct,dz_input
        real(4)psin(51),plin(51)
        real(4)rchla_pass
        integer kk,km

C--- Chlorophyll-equivalent concentration ------------------------------
        atm_wt_N_in_g_per_mol = 14.01

        rchla(1) = (((psin(1)+plin(1)))*atm_wt_N_in_g_per_mol)

C--- Vertical attenuation by depth -------------------------------------
        do kk = 2,km
        rAlpha1    = 0.805
        rAlpha2    = 0.195
        rL1        = 0.043
         rL2        = 0.42

        rchla(kk)   = (((psin(kk)+plin(kk)))  *
     &        (atm_wt_N_in_g_per_mol))

        sfx = sfx_in

        rinfraRed = 0.0

        rchla_pass =  max(0.02, rchla(1))

        rvisible  = sfx *0.55 * rvisible_extinct (rchla_pass,
     &                                             kk-1,dz_input)
        sswin(kk) = rinfraRed + rvisible
        enddo
        sswin(1)  = sfx

        end subroutine light_treatment_by_band

C--- Visible-light attenuation function --------------------------------
        function rvisible_extinct(rchla, kk, dz_input)

         real(4) rVA (2,6), rZA (2,6)
         real(4) rvisible_extinct,rchla,dz_input
         integer kk

C--- Attenuation coefficients: Morel and Antoine (1994) ----------------
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

C--- Polynomial attenuation calculation --------------------------------
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

         rtemp1 = rtempV1*exp(-1.0d0*kk*dz_input/rtempZ1)
         rtemp2 = rtempV2*exp(-1.0d0*kk*dz_input/rtempZ2)
        rvisible_extinct = rtemp1 + rtemp2

        end function rvisible_extinct
