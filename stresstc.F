        subroutine stresstc(uwnd, vwnd, Tair,SST_in, 
     & tauxin, tauyin, vel)
!STRESSTC: computes the neutral wind stress following Smith (1988).
! tau = STRESSTC(vel,z,Tair,rho_a) computes the neutral wind stress given the 
! wind speed and air temperature at height z following Smith (1988),
! J. Geophys. Res., 93, 311-326. Air temperature and density are optional
! inputs. 
!   INPUT:  vel    - wind speed   [m/s]
!           z     - measurement height  [m]
!           Tair    - air temperature (optional) [C]
!           rho_a  - air density (optional)  [kg/m^3]
!   OUTPUT: tau   - wind stress  [N/m^2]
       
       real vel,z,Tair,uwnd,vwnd
       real tauxin,tauyin,cd!,cdntc
       real SST_in
#include "consts.h"
         
       z=10.0 ! 10m above the surface
       vel= sqrt (uwnd**2 + vwnd**2)
        
        !DRAG VALUES FROM SMITH (1988), JGR, TABLE-2. FOR T > 20 DEGREES
        ! AND T < 20
        !ref :/home/user/Desktop/PHd/my_work/coefficints to improve
        !SST.odt

       
!        if(SST_in.ge.20)then
!        if((vel.ge.0.0).and.(vel.lt.2.0))cd = 1.68e-3
!        if((vel.ge.2.0).and.(vel.lt.5.0))cd = 1.45e-3
!        if((vel.ge.5.0).and.(vel.lt.10.0))cd = 1.57e-3
!        if((vel.ge.10.0).and.(vel.lt.15.0))cd = 1.74e-3
!        if((vel.ge.15.0).and.(vel.lt.20.0))cd = 1.93e-3
!        if((vel.ge.20.0).and.(vel.lt.25.0))cd = 2.14e-3
!        end if
!
!        if(SST_in.lt.20.0)then
!        if((vel.ge.0.0).and.(vel.lt.2.0))cd = 1.61e-3
!        if((vel.ge.2.0).and.(vel.lt.5.0))cd = 1.40e-3
!        if((vel.ge.5.0).and.(vel.lt.10.0))cd = 1.52e-3
!        if((vel.ge.10.0).and.(vel.lt.15.0))cd = 1.71e-3
!        if((vel.ge.15.0).and.(vel.lt.20.0))cd = 1.91e-3
!        if((vel.ge.20.0).and.(vel.lt.25.0))cd = 2.12e-3
!        end if 
        
        
       cd = cdntc(vel,z,Tair)! will return value of cd

                
       tauxin = rho_a*(cd*uwnd*abs(uwnd))
       tauyin = rho_a*(cd*vwnd*abs(vwnd))
        
        
        
      
       return        
       end subroutine stresstc 

        

      function cdntc(vel,z,Tair)
        ! CTDTC: computes the neutral drag coefficient following Smith(1988).
        ! [cd,u10]=CDNTC(vel,z,Tair) computes the neutral drag coefficient and
        ! wind speed at 10m given the wind speed and air temperature at height z
        ! following Smith (1988), J. Geophys. Res., 93, 311-326.
        !
        ! INPUT:   vel - wind speed  [m/s]
        !          z - measurement height [m]
        !          Tair - air temperature (optional)  [C]
        !
        ! OUTPUT:  cd - neutral drag coefficient at 10m
        !          u10 - wind speed at 10m  [m/s]

        real kappa, cdntc
        integer mm

#include "consts.h"
        kappa=0.4
                ! iteration endpoint
      tol=.00001
      visc= 1.326e-5*(1 + 6.542e-3*Tair + 8.301e-6*Tair**2 - 4.84e-9*Tair**3)
        
        ! remove any vel==0 to prevent division by zero
      if (vel .eq. 0.0) then 
      vel=0.1 ! replacing al zero points with stress of 0.1
      end if
                ! initial guess
      ustaro=0.0
      ustarn=0.036*vel
                ! iterate to find z0 and ustar
10    if (ustarn-ustaro .gt. tol) then
      mm = 1
      else
      mm=0
      endif 
      if (mm .gt. 0) then
      ustaro=ustarn
      z0=Charnock_alpha*ustaro**2/g + (R_roughness*visc)/ustaro
      ustarn=vel*(kappa/log(z/z0))
      goto 10 
      endif
          
      sqrcd=kappa/log((10)/z0)
      cdntc=sqrcd**2
         
      end function

