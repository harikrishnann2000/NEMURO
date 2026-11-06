! Mixed-layer model routine
      SUBROUTINE MT1 (tempin, saltin, uvelin, vvelin,
     &  psin, plin, zsin,zlin,zpin,no3in,nh4in,pomin,domin,sioh4in,
     &  opalin,cain,caco3in,tco2in,talkin,
     &  tauxin,tauyin,
     &  qnet,pme,rLATin,dt_rg,SM,SH,sswin,salt,cp_w,sal_clim,
     &  nonlocal,dz_input)

! Declarations and model state
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
        real dz_input
	real NO31, NO3N, NH41, NH4N,cp_w , nonlocal
        real add
        real sal_clim
        real :: dq
        integer :: nloop
        logical :: NAN_check

! Grid configuration and initial conditions
         NAZ=50
      NAZ1=NAZ+1
      NAZ2=NAZ+2
      DZ=dz_input
      DZZ=DZ*DZ
      DZ2=DZ*2
      Z0=0.1
      Y0=1.0
      DO 10 M=1,NAZ1
        UN(M)=uvelin(M)
        VN(M)=vvelin(M)
        TEN(M)   = tempin(M)
        SAN(M)   = saltin(M)
	PSN(M)   = psin(M)
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
        TE1(0)  = TE1(1)
        SA1(0)  = SA1(1)
        U1(0)   = U1(1)
        V1(0)   = V1(1)
        TE1(51) = TE1(50)
        SA1(51) = SA1(50)
        U1(51)  = U1(50)
        V1(51)  = V1(50)
   10 CONTINUE
      DO 11 M=1,NAZ1
   11 Z(M)=DZ*(1-M)-Z0
      DO 12 M=1,NAZ1
   12 Z1(M)=DZ*(1-M)

! Time integration setup
      DAY=86400.0
       DT=dt_rg/4
      nloop=dt_rg/DT
      DTZ=DT/DZZ

! Physical constants and Coriolis parameter
      RO=1025.0
      R5=10.0**(-5)
      Z0D=Z0/DZ
      GA=9.8
        deg2rad = 3.14/180.0
        F0= 2.0*(0.729e-4)*sin(rLATin*deg2rad)
       FT=F0*DT

! Turbulence model coefficients
      A1=0.92
      A2=0.74
      B1=16.6
      B2=10.1
      C1=0.08
      G1=1.0/3.0-2.0*A1/B1
      G2=(B2+6.0*A1)/B1
      E1=3.0*A2
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
      BS=(-76.0)*10.0**(-5)
      XK=0.4
      XK2=XK*XK
      G3=G3*XK2
      G4=G4*SQRT(B1*3.0*A1*(G1-C1))*XK2*Z0
      XK=(-1.0)*0.5*XK
      AL=0.04
      H2=ALOG((1.0+Z0D)/(2.0+Z0D))
      H3=ALOG(Z0D/(1.0+Z0D))
      H1=H3-H2

! Surface interpolation coefficients
      HH=H2/DZ+1.0/Z0
      HJ1=(H2+1.0/(Z0D+0.5))/(H2-H3)
      HJ2=1.0-HJ1

! Main time-integration loop
      K1=0
    2 K1=K1+1
      T=T+DT
      TWX = tauxin
      TWY = tauyin

! Vertical mixing coefficients
      TX=TWX/RO
      TY=TWY/RO
      TXY=SQRT(TX*TX+TY*TY)
      DO 22 M=1,NAZ
      R(M)=0.0
      SM(M)=0.0
      SH(M)=0.0
   22 CONTINUE
      DO 23 M=1,NAZ
      SZ=0.5*(SA1(M)+SA1(M+1))-30.0
      BT=((1.3-0.025*SZ)*(TE1(M)+TE1(M+1))*0.5 + (3.3+0.325*SZ))*R5
      IF(M.EQ.1) GO TO 24
      BX=BT*(TE1(M-1)-TE1(M+1))+BS*(SA1(M-1)-SA1(M+1))
      UV=(U1(M-1)-U1(M+1))**2+(V1(M-1)-V1(M+1))**2
      GO TO 25
   24 BX=HJ1*(BT*(TE1(1)-TE1(2))+BS*(SA1(1)-SA1(2)))+
     & HJ2*(BT*(TE1(2)-TE1(3))+BS*(SA1(2)-SA1(3)))
      UU=HJ1*(U1(1)-U1(2))+HJ2*(U1(2)-U1(3))
      VV=HJ1*(V1(1)-V1(2))+HJ2*(V1(2)-V1(3))
      UV=UU**2+VV**2
   25 IF(UV.LT.10.0**(-8)) GO TO 23
      BX=BX/UV
      BX=BX*GA*DZ*2
      X1=F2+F3*BX
      X2=F4*BX
      X3=SQRT(X1**2-2.0*F1*X2)
      R(M)=(X1-X3)/F1
      SH(M)=E1*(G1-E2*R(M))/(1.0-R(M))
      IF(R(M) .GT. GC) SH(M)=0.0
      SM(M)=SH(M)*E3*(E4-E5*R(M))/(E6-E7*R(M))
   23 CONTINUE

! Iterative mixing-length calculation
      J=0
   30 J=J+1
      Q1=0.0
      Q2=0.0
      DO 31 M=1,NAZ
      ZK=XK*(Z(M)+Z(M+1))
      YQ(M)=ZK/(1.0+ZK/Y0)
      QX=(U1(M)-U1(M+1))**2+(V1(M)-V1(M+1))**2
      QX=QX*B1*SM(M)*(1.0-R(M))/DZZ
      QX=SQRT(QX)
      Q(M)=QX*YQ(M)
      Q1=Q1+Q(M)
      Q2=Q2-Q(M)*0.5*(Z(M)+Z(M+1))
   31 CONTINUE
      IF(Q1.EQ.0.0) GO TO 32
      YA=Y0
      Y0=AL*Q2/Q1
      YA=ABS(Y0-YA)
      IF(J.EQ.10) GO TO 32
      IF(YA.GT.0.01) GO TO 30
   32 DO 33 M=1,NAZ
      SM(M)=SM(M)*YQ(M)*Q(M)
      SH(M)=SH(M)*YQ(M)*Q(M)
   33 CONTINUE
        SM(1) = max(1e-5,SM(1))
        SH(1) = max(1e-5,SH(1))
        SM(2) = max(1e-5,SM(2))
        SH(2) = max(1e-5,SH(2))

! Smooth the mixing coefficients
      DO 34 M=1,NAZ1
      SM1(M)=SM(M)
      SH1(M)=SH(M)
   34 CONTINUE
      DO 35 M=2,NAZ
   35 CONTINUE

! Update momentum
      DO 40 M=1,NAZ
        UN(M) = UN(M) + FT*V1(M)- (U1(M)*1e-3)
        VN(M) = VN(M) - FT*U1(M)- (V1(M)*1e-3)
   40 CONTINUE

! Surface velocity boundary condition
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

! Apply surface heat, freshwater, and light forcing
        TEN(1)=TEN(1) + qnet*DT
        SAN(1)=SAN(1) + salt*DT  +
     &       (sal_clim-SA1(1))*DT/(86400*3.0)
       do M=2,NAZ
        if (M .eq. 2) rpar_return = 1.0
        if (M .ne. 2) rpar_return = 1.0
        TEN(M) = TEN(M) + ((sswin(M-1) - sswin(M))*DT)
     &                     / (RO * cp_w * DZ * rpar_return )
        end do

! Implicit vertical mixing of model tracers
        call implicit_vmix (SM+nonlocal, DZ, DT, NAZ1, UN, UN)
        call implicit_vmix (SM+nonlocal, DZ, DT, NAZ1, VN, VN)
        call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, TEN, TEN)
        call implicit_vmix (SH+nonlocal, DZ, DT, NAZ1, SAN, SAN)
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

! Check state variables for invalid values
        DO M=1,NAZ
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

! Restore previous values if invalid values are found
       DO 59 M=1,NAZ
        if (i .ne. 0) then
      UN(M)=U1(M)
      VN(M)=V1(M)
      TEN(M)=TE1(M)
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

! Save updated values for the next timestep
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
  101 FORMAT(1H ,E15.7)
  102 FORMAT(1H ,51F8.2)

! Return the updated profiles
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
        RETURN
        STOP
        END
