!        program MTI

      SUBROUTINE MT1 (tempin, saltin, uvelin, vvelin, 
     &  psin, plin, zsin,zlin,zpin,no3in,nh4in,pomin,domin,sioh4in,
     &  opalin,cain,caco3in,tco2in,talkin,
     &  tauxin,tauyin, 
     &  qnet, pme, rLATin,dt_rg,SM,SH,sswin,salt,cp_w,sal_clim,nonlocal)

! ONE DIMENSIONAL MIXED LAYER MODEL
      DIMENSION UN(51),VN(51),TEN(51),SAN(51),
     1 U1(0:51),V1(0:51),TE1(0:51),SA1(0:51),Z(51),Z1(51),
     2 Q(51),SM(51),SH(51),R(51),YQ(51),SM1(51),SH1(51),
     3 PS1(51), PL1(51), ZS1(51), ZL1(51), ZP1(51),
     4 NO31(51), DOM1(51), POM1(51), NH41(51), SIOH41(51),
     5 Opal1(51), Ca1(51), CaCO31(51), TCO21(51), TALK1(51),
     6 PSN(51), PLN(51), ZSN(51), ZLN(51), ZPN(51),
     7 NO3N(51), DOMN(51), POMN(51), NH4N(51), SIOH4N(51),
     8 OpalN(51), CaN(51), CaCO3N(51), TCO2N(51), TALKN(51),
     9 TEN1(51),sswin(51),rlightk(51), nonlocal(51)  

        real tempin(51), saltin(51)
        real uvelin(51), vvelin(51)
	real psin(51), plin(51), zsin(51),zlin(51),zpin(51)
	real no3in(51),nh4in(51)
	real pomin(51),domin(51),sioh4in(51)
	real opalin(51),cain(51),caco3in(51),tco2in(51),talkin(51)
        real tauxin, tauyin, qnet, pme, rLATin, dt_rg
	real NO31, NO3N, NH41, NH4N,cp_w , nonlocal
        
        real add                                         
        real sal_clim
        real :: dq
        integer :: nloop
        logical :: NAN_check
 


! VARIABLES
! U,V=VELOCITY,TE=TEMP.,SA=SAL.,Z=DEPTH
! WITH N=NEW VALUE AT EACH TIME STEP, 1=OLD VALUE
! Q=TURBULENT KINETIC ENERGY
! SM,SH=MIXING COEFFICIENTS
! R=RICHARDSON NUMBER,YQ=MIXING LENGTH
! NUMBER OF GRIDS
         NAZ=50
      NAZ1=NAZ+1
      NAZ2=NAZ+2  
! GRID SIZE (THICKNESS IN METER)
      DZ=5.0
      DZZ=DZ*DZ
      DZ2=DZ*2
! ROUGHNESS
      Z0=0.1
! INITIAL CONDITION
      Y0=1.0
      DO 10 M=1,NAZ1
        UN(M)=uvelin(M)        !0.0
        VN(M)=vvelin(M)  
!    UN(M)=0.0
!      VN(M)=0.0
! STABLE STRATIFICATION
        TEN(M)   = tempin(M) !20.0 -0.1*(M-1)
        SAN(M)   = saltin(M)
	PSN(M)   = psin(M)   !micromolN/L
	PLN(M)   = plin(M)
	ZSN(M)   = zsin(M)
	ZLN(M)   = zlin(M)
	ZPN(M)   = zpin(M)
	NO3N(M)  = no3in(M)
	NH4N(M)  = nh4in(M)
	DOMN(M)  = domin(M)
	POMN(M)  = pomin(M)
	SIOH4N(M)= sioh4in(M)
	OpalN(M) = opalin(M)
	TCO2N(M) = tco2in(M)
	TALKN(M) = talkin(M)
	CaN(M)   = cain(M)
	CaCO3N(M)= caco3in(M)
 
!     TEN(M)=20.0-0.01*(M-1)
!     SAN(M)=33.0-0.000*(M-1)
      U1(M)=UN(M)
      V1(M)=VN(M)
      TE1(M)=TEN(M)
      SA1(M)=SAN(M)
	PS1(M)= PSN(M)
	PL1(M)= PLN(M)
	ZS1(M)=ZSN(M)
	ZL1(M)=ZLN(M)
	ZP1(M)= ZPN(M)
	NO31(M) = NO3N(M)
	NH41(M) =NH4N(M)
	DOM1(M) = DOMN(M)
	POM1(M) = POMN(M)
	SIOH41(M)= SIOH4N(M)
	Opal1(M)=OpalN(M)
	TCO21(M)=TCO2N(M)
	TALK1(M)=TALKN(M)
	Ca1(M)=CaN(M)
	CaCO31(M)=CaCO3N(M)
        
        TE1(0)  = TE1(1) !24_may
        SA1(0)  = SA1(1) !24_may
        U1(0)   = U1(1)  !24_may
        V1(0)   = V1(1)  !24_may
        TE1(51) = TE1(50)!24_may
        SA1(51) = SA1(50)!24_may
        U1(51)  = U1(50) !24_may
        V1(51)  = V1(50) !24_may
	
	
   10 CONTINUE
	
      DO 11 M=1,NAZ1
   11 Z(M)=DZ*(1-M)-Z0
      DO 12 M=1,NAZ1
   12 Z1(M)=DZ*(1-M)
! SECONDS IN ONE DAY
      DAY=86400.0
! TIME STEP
!      DT=0.002*DAY  !by shika
       DT=dt_rg/4
      nloop=dt_rg/DT
      DTZ=DT/DZZ
      !loop= DAY/DT  ! not using anywhere...
        
       

! WATER DENSITY
      RO=1025.0 
      R5=10.0**(-5) 
      Z0D=Z0/DZ

! GRAVITY ACCELERATION
      GA=9.8

! CORIOLIS PARAMETER
        deg2rad = 3.14/180.0
        F0= 2.0*(0.729e-4)*sin(rLATin*deg2rad) 
!	write(*,*) rLATin,F0
!      F0=1.2*10.0**(-4) 
       FT=F0*DT
! COEFFICIENTS IN THE TURBULENT MODEL
      A1=0.92
      A2=0.74
      B1=16.6
      B2=10.1
      C1=0.08
      G1=1.0/3.0-2.0*A1/B1  !Gamma 1 
      G2=(B2+6.0*A1)/B1     !Gamma2
      E1=3.0*A2             !Sh0
      E2=G1+G2 
      GC=G1/E2
      E3=A1/A2
      E4=B1*(G1-C1)
      E5=E4+6.0*A1+3.0*A2
      E6=B1*G1
      E7=B1*E2-3.0*A1
      F1=2.0*A1*E5
      F2=A1*E4
      F3=A2*E7
      F4=A2*E6
      G3=SQRT(B1*(3.0*A1*(G1-C1))**3) 
      G4=3.0*A2*G1
! SALINE EXPANSION COEFFICIENT
      BS=(-76.0)*10.0**(-5)
! VON KARMAN CONSTANT
      XK=0.4
      XK2=XK*XK
      G3=G3*XK2
      G4=G4*SQRT(B1*3.0*A1*(G1-C1))*XK2*Z0
      XK=(-1.0)*0.5*XK
      AL=0.04 !!0.04!AL can be change to 0.05,0.1
!ref:/home/user/Desktop/PHd/my_work/penetrative_solar_radiation_treatment/references/Mellor_and_Yamada_in_RAMS.pdf
! COEFFICIENTS FOR NEAR SURFACE INTERPOLATION
      H2=ALOG((1.0+Z0D)/(2.0+Z0D))
      H3=ALOG(Z0D/(1.0+Z0D))
      H1=H3-H2
      HH=H2/DZ+1.0/Z0
      HJ1=(H2+1.0/(Z0D+0.5))/(H2-H3)
      HJ2=1.0-HJ1
! TIME INTEGRATION STARTS
!      K=0
!    1 K=K+1
      K1=0
    2 K1=K1+1
! TIME INCREASING
      T=T+DT
! WIND STRESS
      TWX = tauxin
      TWY = tauyin
        
! TWX=1.0
! TWY=1.0
! VERTICAL MIXING
      TX=TWX/RO
      TY=TWY/RO
      TXY=SQRT(TX*TX+TY*TY)
      DO 22 M=1,NAZ 
      R(M)=0.0 
      SM(M)=0.0 
      SH(M)=0.0
   22 CONTINUE
      
      DO 23 M=1,NAZ
      SZ=0.5*(SA1(M)+SA1(M+1))-30.0 !why 0.5 here??
! THERMAL EXPANSION COEFFICIENT
      BT=((1.3-0.025*SZ)*(TE1(M)+TE1(M+1))*0.5 + (3.3+0.325*SZ))*R5
      IF(M.EQ.1) GO TO 24 
! STRATIFICATION
      BX=BT*(TE1(M-1)-TE1(M+1))+BS*(SA1(M-1)-SA1(M+1)) !central_diff
! VERTICAL SHEAR
      UV=(U1(M-1)-U1(M+1))**2+(V1(M-1)-V1(M+1))**2   !24_may
      GO TO 25
! NEAR SURFACE STRATIFICATION AND SHEAR
   24 BX=HJ1*(BT*(TE1(1)-TE1(2))+BS*(SA1(1)-SA1(2)))+ 
     & HJ2*(BT*(TE1(2)-TE1(3))+BS*(SA1(2)-SA1(3)))
      UU=HJ1*(U1(1)-U1(2))+HJ2*(U1(2)-U1(3))
      VV=HJ1*(V1(1)-V1(2))+HJ2*(V1(2)-V1(3))
      UV=UU**2+VV**2
   25 IF(UV.LT.10.0**(-8)) GO TO 23 
! R,SH AND SM CALCULATED
      BX=BX/UV
      BX=BX*GA*DZ*2 
      X1=F2+F3*BX
      X2=F4*BX
      X3=SQRT(X1**2-2.0*F1*X2)
      R(M)=(X1-X3)/F1    !gradient Richardson No.
      SH(M)=E1*(G1-E2*R(M))/(1.0-R(M))            !eq. 7a
      IF(R(M) .GT. GC) SH(M)=0.0                  ! GC -crti Richa no.
      SM(M)=SH(M)*E3*(E4-E5*R(M))/(E6-E7*R(M))    !eq. 7b
   23 CONTINUE



! ITERATIVE CALCULATION OF MIXING LENGTH
      J=0
   30 J=J+1 
      Q1=0.0
      Q2=0.0
      DO 31 M=1,NAZ
      ZK=XK*(Z(M)+Z(M+1)) 
      YQ(M)=ZK/(1.0+ZK/Y0) !YQ= mixing length.
      QX=(U1(M)-U1(M+1))**2+(V1(M)-V1(M+1))**2
      QX=QX*B1*SM(M)*(1.0-R(M))/DZZ                !eq.2 
      QX=SQRT(QX)
      Q(M)=QX*YQ(M) 
      Q1=Q1+Q(M)
      Q2=Q2-Q(M)*0.5*(Z(M)+Z(M+1)) ! 0.5 -> 0.1 then to 1.0
!it was 0.5 then changed it in to 0.1 again changed in to 1.0
   31 CONTINUE
      IF(Q1.EQ.0.0) GO TO 32
      YA=Y0
      Y0=AL*Q2/Q1                   !eq. 6
      YA=ABS(Y0-YA)
      IF(J.EQ.10) GO TO 32
      IF(YA.GT.0.01) GO TO 30 
   32 DO 33 M=1,NAZ
        
      SM(M)=SM(M)*YQ(M)*Q(M) !+ add   !DELETE
      SH(M)=SH(M)*YQ(M)*Q(M) !+ add   !DELETE

        

        
   33 CONTINUE
        
        SM(1) = max(1e-5,SM(1))!24_may
        SH(1) = max(1e-5,SH(1))!24_may
        SM(2) = max(1e-5,SM(2))!24_may
        SH(2) = max(1e-5,SH(2))!24_may
        
        
!************************ SMOOTHING*********************************************
      DO 34 M=1,NAZ1
      SM1(M)=SM(M)
      SH1(M)=SH(M)
   34 CONTINUE
      DO 35 M=2,NAZ
!!      SM(M)=(SM1(M-1)+2.0*SM1(M)+SM1(M+1))/4.0 !included march 15
!onwards 
!      SH(M)=(SH1(M-1)+2.0*SH1(M)+SH1(M+1))/4.0!""

   35 CONTINUE
        

!******************************************************************************
!! TIME INTEGRATION OF U,V,TE AND SA

		
      DO 40 M=1,NAZ
        UN(M) = UN(M) + FT*V1(M)- (U1(M)*1e-3) 
        VN(M) = VN(M) - FT*U1(M)- (V1(M)*1e-3) 
!        (U1(M)*1e-3),(V1(M)*1e-3) ! these 2 terms are linear damping, i
!        added later. it was not there in the code.
!        ie friction.
   40 CONTINUE
        


! CALCULATION OF THE TOP GRID VALUES
      UO=0.0
      VO=0.0
      UVS=0.01
      UV=SQRT(TXY/G3)
      IF(UV.LT.0.00001) GO TO 41 
      UO=TX/G3/UV 
      VO=TY/G3/UV 
      UVS=SQRT(UO**2+VO**2) 
   41 UC2=((UN(2)-UN(3))/DZ-UO/Z0)/HH 
      UN(1)=UN(2)+UC2*H1+UN(2)-UN(3)
      VC2=((VN(2)-VN(3))/DZ-VO/Z0)/HH 
      VN(1)=VN(2)+VC2*H1+VN(2)-VN(3)

!***********************************FROM SHIKA**************************
! SURFACE BOUNDARY CONDITION
! COOLING 0.1 DEG PER DAY
!        qf = qnet ! read value of heat flux here  
!        mf = pme  ! read value of freshwater flux in here in m/s
!        cp = 3958.0
!        h = DZ
!       dq= qf*DT/(RO*cp*h)
!       salt = mf*RO*35.0 ! Sref= 35.0   
!        dq=0.0
!        salt=0.0  
!***************************************************************************
        
!******************* PROVIDING COOLING***********************************

! qnet and rlight(M) will pass from main        
        
        
        TEN(1)=TEN(1) + qnet*DT  
                                 
        SAN(1)=SAN(1) + salt*DT  + 
     &       (sal_clim-SA1(1))*DT/(86400*3.0) 
                                 !ten(1)=te1(1)+dq ! heat flux, h- thickness of model layer, q-net flux
                                 !san(1)=sa1(1) + salt_change
!        write(*,*) sal_clim
       do M=2,NAZ 
                                 
        if (M .eq. 2) rpar_return = 1.0!0.45
        if (M .ne. 2) rpar_return = 1.0
        TEN(M) = TEN(M) + ((sswin(M-1) - sswin(M))*DT) ! for nemuro !light_treatment rlightk instead of sswin
     &                     / (RO * cp_w * DZ * rpar_return )! - for
                                              !nemuro light treatmnt 
                                               
        
        end do    
                                                    
        call implicit_vmix (SM+nonlocal, DZ, DT, NAZ1, UN, UN) !24_may
        call implicit_vmix (SM+nonlocal, DZ, DT, NAZ1, VN, VN) !24_may
        call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, TEN, TEN)!24_may
        call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, SAN, SAN)!24_may
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, PSN, PSN)   
        call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, PLN, PLN)   
        call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, ZSN, ZSN)  
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, ZLN, ZLN)  
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, ZPN, ZPN)   
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, NO3N, NO3N) 
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, NH4N, NH4N) 
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, POMN, POMN) 
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, DOMN, DOMN)  
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, SIOH4N, SIOH4N) 
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, OpalN, OpalN)   
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, CaN, CaN)        
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, CaCO3N, CaCO3N) 
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, TCO2N, TCO2N)   
	call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, TALKN, TALKN)   
!*******************************************************************************************        
!        TEN(1)=TEN(1) + qnet*DT
!        SAN(1)=SAN(1) + salt*DT
        
!        do M=1,NAZ
!         TEN(M) = TEN(M) + ((Q_zk(M) - Q_zk(M+1)  )*DT)
!     &                     / (RO * cp_w * DZ ) !here heat absorbed by
!!planktons not acconting .. dbt
                                                  
!       enddo
!*****************************************************************************8

        DO M=1,NAZ !50 or 51, why boundary avoided here. check it...
        i = 0
        if (UN(M) .ne. UN(M)) i = i + 1
        if (VN(M) .ne. VN(M)) i = i + 1
        if (TEN(M) .ne. TEN(M)) i = i + 1
        if (SAN(M) .ne. SAN(M)) i = i + 1
        if (PSN(M) .ne. PSN(M)) i = i + 1
       if (PLN(M) .ne. PLN(M)) i = i + 1
        if (ZSN(M) .ne. ZSN(M)) i = i + 1
        if (ZLN(M) .ne. ZLN(M)) i = i + 1
        if (ZPN(M) .ne. ZPN(M)) i = i + 1
        if (NO3N(M) .ne. NO3N(M)) i = i + 1
        if (NH4N(M) .ne. NH4N(M)) i = i + 1
        if (DOMN(M) .ne. DOMN(M)) i = i + 1
        if (POMN(M) .ne. POMN(M)) i = i + 1
        if (CaN(M) .ne. CaN(M)) i = i + 1
        if (CaCO3N(M) .ne. CaCO3N(M)) i = i + 1
        if (SIOH4N(M) .ne. SIOH4N(M)) i = i + 1
        if (OpalN(M) .ne. OpalN(M)) i = i + 1
        if (TCO2N(M) .ne. TCO2N(M)) i = i + 1
        if (TALKN(M) .ne. TALKN(M)) i = i + 1
        ENDDO

       DO 59 M=1,NAZ   
        if (i .ne. 0) then
      UN(M)=U1(M)
      VN(M)=V1(M)
      TEN(M)=TE1(M)  ! ASK SHIKHA
      SAN(M)=SA1(M)
	PSN(M) = PS1(M)
	PLN(M) = PL1(M)
	ZSN(M) = ZS1(M)
	ZLN(M) = ZL1(M)
	ZPN(M) = ZP1(M)
	NO3N(M) = NO31(M)
	NH4N(M) = NH41(M)
	DOMN(M) = DOM1(M)
	POMN(M) = POM1(M)
	CaN(M)  = Ca1(M)
	CaCO3N(M)=CaCO31(M)
	SIOH4N(M)=SIOH41(M)
	OpalN(M)=Opal1(M)
	TCO2N(M)=TCO21(M)
	TALKN(M)=TALK1(M)
        endif
   59 CONTINUE
        
! REPLACE THE OLD VALUES WITH THE NEW VALUES
        
       DO 52 M=1,NAZ   
        if (i .eq. 0) then
      U1(M)=UN(M)
      V1(M)=VN(M)
      TE1(M)=TEN(M)  
      SA1(M)=SAN(M)
	PS1(M) = PSN(M)
	PL1(M) = PLN(M)
	ZS1(M) = ZSN(M)
	ZL1(M) = ZLN(M)
	ZP1(M) = ZPN(M)
	NO31(M) = NO3N(M)
	NH41(M) = NH4N(M)
	DOM1(M) = DOMN(M)
	POM1(M) = POMN(M)
	Ca1(M)  = CaN(M)
	CaCO31(M)=CaCO3N(M)
	SIOH41(M)=SIOH4N(M)
	Opal1(M)=OpalN(M)
	TCO21(M)=TCO2N(M)
	TALK1(M)=TALKN(M)
        endif
   52 CONTINUE
! TIME STEPS FOR 0.1 DAY
! PRINT OUT THE RESULTS
!!! change by shikha to involve just one loop for all
!      IF(K1.LT.50) GO TO 2 
!         WRITE(6,101) T/DAY
  101 FORMAT(1H ,E15.7)
!          WRITE(6,102) (10.0*UN(M),M=1,51,1)
  102 FORMAT(1H ,51F8.2)
!           WRITE(6,102) (10.0*VN(M),M=1,51,1)
!          WRITE(6,102) (TEN(M),M=1,51,1)
!          WRITE(6,102) (SAN(M)-30.0,M=1,51,1)
!          WRITE(6,102) (100.0*SM(M),M=1,51,1)
!         WRITE(6,102) (100.0*SH(M),M=1,51,1)
! TIME STEPS FOR 0.1 DAY X 50
        
               
!      IF(K.LT. nloop) GO TO 1
      IF(K1.LT. nloop) GO TO 2

        DO M=1,NAZ
        tempin(M) = TEN(M)
        saltin(M) = SAN(M)
        uvelin(M) = UN(M)
        vvelin(M) = VN(M)
	psin(M)=PSN(M)
	plin(M)=PLN(M)
	zsin(M)=ZSN(M)
	zlin(M)=ZLN(M)
	zpin(M)=ZPN(M)
	no3in(M)=NO3N(M)
	nh4in(M)=NH4N(M)
	domin(M)=DOMN(M)
	pomin(M)=POMN(M)
	sioh4in(M)=SIOH4N(M)
	opalin(M)=OpalN(M)
	cain(M)=CaN(M)
	caco3in(M)=CaCO3N(M)
	tco2in(M)=TCO2N(M)
	talkin(M)=TALKN(M)
        ENDDO
        
!************************************
!        do M=1,NAZ  
!        print*,PSN(M),PLN(M),ZSN(M),ZLN(M),ZPN(M),
!     &  NO3N(M),NH4N(M),DOMN(M),POMN(M),SIOH4N(M),
!     &  opalN(M),CaN(M),CaCO3N(M),TCO2N(M),TALKN(M)     
!        end do
!**********************************
        RETURN
        STOP 
        END
        
! function for latent heat calculation , doubt declare everything in above also       
