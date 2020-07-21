      integer, parameter :: nElements = 1
      integer, parameter :: nintpts = 8
      integer, parameter :: nsdv  = 93
      integer, parameter :: maxphasenum  = 5
      integer :: writeflag
      real*8 :: kFp, kgausscoords, kcurlFp
          
      COMMON/UMPS/kFp(nElements, nintpts, 9),
     +    kgausscoords(nElements,nintpts,3),
     +    kcurlFp(nElements, nintpts, 9),
     +    writeflag(maxphasenum)

      iphase = int(props(1))

      SELECT CASE(iphase)
      CASE(0) ! Example HCP System
          nSys = L0
          crystalsys = 0 ! 0 will load HCP slip systems
          
          TAUCSSD = 240
          TauPyra = 0.0011*(TEMP-273)**2-1.3801*(TEMP-273)+867.18
          crssratio = 3.01    
          
          burger1 = 3.20E-4
          caratio = 1.587

          ! Elastic Constants taken from (Hasija et al., 2003)
          E1 = 84.745E+3; E3 = 119.789E+3
          v12 = 0.46 
          v13 = 0.22        
          G12 = E1/(2.0*(1.0+v12)); G13 = 40000     

          gammast = 1.23 !increase if simulated hardening is lower than experimental

          rhossd = STATEV(54)
          rhoGND= STATEV(26)
          rho=rhoGND+rhossd
          XTAUC1 = TAUCSSD + (0.5 * G12 * burger1 * sqrt(rho)) !Taylor hardening equation
         
          XTAUC2 = crssratio*XTAUC1 !XTAUC1 is the CRSS of the weakest slip system (prismatic for Hf)
		                            !Then apply crssratio to find CRSS (XTAUC2) for slip system that is harder to activate
          burger2 = caratio*burger1

          alpha1 = 9.5D-6; alpha2 = alpha1; alpha3 = 0.5895*alpha1


          allocate(xNorm(L0,M),xDir(L0,M),tau(L0),gammadot(L0),
     +    gndold(L0),burgerv(L0),tauc(L0),ids(L0),tau2(L0),STAT=ialloc)
         
          burgerv(1:12) = burger1; burgerv(13:24) = burger2
          tauc(1:3) = XTAUC1*1.0D0;     !  <a>  Basal
          tauc(4:6) = XTAUC1*1.015D0;   !  <a>  Prismatic
          tauc(7:12) = XTAUC2           !  <a>  Pyramidal
          tauc(13:24) = XTAUC2          ! <c+a> Pyramidal

          rhossdm    = 5.0
          xhelmholtz = 7.64582499068784D-20*0.96  
          xfreq      = 1.0e+11
          gammazero  = 6.0e-4*0.22 !work conjugate of activation volume (appears in activation volume bit) 
          xboltzman  = 1.381e-23
          xtemp      = 293.0
          rhoinitial = 0.01 
            
      case(1) ! Example BCC System
          nSys = L1
          crystalsys = 1 ! 1 will load BCC slip systems
          XTAUC = 240.0E+9 
          burger = 3.20E-4
          
          !bcc properties Kim and Rokhlin (2009) J.Acoust.Soc.Am.
          E1 = 3.2024e+04; E3 = E1 
          v12 = 0.4556; v13 = v12
          G12 = 54900; G13 = G12
          
          alpha1 = 9.5e-6; alpha2 = alpha1; alpha3 = 0.5895*alpha1
          
          allocate(xNorm(L1,M),xDir(L1,M),tau(L1),gammaDot(L1),
     +    gndold(L1),burgerv(L1),tauc(L1),ids(L1),tau2(L1),STAT=ialloc)
         
          gammast = 1.23
          burgerv = burger
          tauc = xtauc

          rhossdm    = 0.01
          xhelmholtz = 2.6e-20
          xfreq      = 1.0e+11
          gammazero  = 1e-3
          xboltzman  = 1.381e-23
          xtemp      = 293.0
          rhoinitial = 0.01
     
      case(2) ! RR1000 FCC
          nSys = L2
          crystalsys = 2 ! 2 will load FCC slip systems
          XTAUCC=450 
          E1 = 210e+3; E3 = E1
          v12 = 0.28; v13 = v12
          G12 = 90e+3; G13 = G12  !
          
          alpha1 = 13.0e-6; alpha2 = alpha1; alpha3 = alpha1
          
          gammast = 150
          burger = 3.5072e-4
          rhossd = STATEV(54)
          rhognd = STATEV(26)
          rho = rhossd+rhognd
          xtauc = XTAUCC + (1.0 * G12 * burger * sqrt(rho))
         
          allocate(xNorm(L2,M),xDir(L2,M),tau(L2),gammaDot(L2),
     +    gndold(L2),burgerv(L2),tauc(L2),ids(L2),tau2(L2),STAT=ialloc)
         
          burgerv = burger
          tauc = xtauc

          rhossdm    = 0.05
          xhelmholtz = 3.456e-20
          xfreq      = 1.0e+11
          gammazero  = 8.33e-6
          xboltzman  = 1.381e-23
          xtemp      = 293.0
          rhoinitial = 0.05

      case(3) ! Oxide - Wang et al. 1991
          nSys = L2
          crystalsys = 2 ! 2 will load FCC slip systems
          XTAUCC=10e20
          E1 = 280e+3; E3 = E1
          v12 = 0.295; v13 = v12
          G12 = 90e3; G13 = G12  !
          
          alpha1 = 13.0e-6; alpha2 = alpha1; alpha3 = alpha1
          
          gammast = 150
          burger = 3.5072e-4
          rhossd = STATEV(54)
          rhognd = STATEV(26)
          rho = rhossd+rhognd
          xtauc = XTAUCC + (1.0 * G12 * burger * sqrt(rho))
         
          allocate(xNorm(L2,M),xDir(L2,M),tau(L2),gammaDot(L2),
     +    gndold(L2),burgerv(L2),tauc(L2),ids(L2),tau2(L2),STAT=ialloc)
         
          burgerv = burger
          tauc = xtauc

          rhossdm    = 0.05
          xhelmholtz = 3.456e-20
          xfreq      = 1.0e+11
          gammazero  = 8.33e-6
          xboltzman  = 1.381e-23
          xtemp      = 293.0
          rhoinitial = 0.05

      case default
      WRITE(6,*)
      WRITE(6,*)"Not sure what crystal type. Material constants."
      END SELECT

C     Write out important model parameters to .dat file for each phase
      if (kinc == 1 .and. npt == 1) THEN
        if (writeflag(iphase+1) .EQ. 0) THEN
          WRITE(6,*)
          WRITE(6,*) "***************************"
          WRITE(6,*) "**** PHASE ", iphase
          WRITE(6,*) "***************************"
          WRITE(6,*) "**** ELASTIC CONSTANTS ****"
          WRITE(6,*) "***************************"
          WRITE(6,*) "**** E1: ", E1
          WRITE(6,*) "**** E3: ", E3
          WRITE(6,*) "**** v12: ", v12
          WRITE(6,*) "**** v13: ", v13
          WRITE(6,*) "**** G12: ", G12
          WRITE(6,*) "**** G13: ", G13
          WRITE(6,*) "***************************"
          WRITE(6,*) "****  SLIP PROPERTIES  ****"
          WRITE(6,*) "***************************"
          WRITE(6,*) "**** Number of slip systems: ", nSys
          WRITE(6,*) "****    CRSS VALUES    ****"
          WRITE(6,*) TAUC
          WRITE(6,*) "***************************"
          WRITE(6,*) "****  Burgers vectors  ****"
          WRITE(6,*) burgerv
          WRITE(6,*) "***************************"
          WRITE(6,*) "**** Rho SSDm: ", rhossdm   
          WRITE(6,*) "**** Helmholtz: ", xhelmholtz
          WRITE(6,*) "**** Freq: ", xfreq     
          WRITE(6,*) "**** Gamma0: ", gammazero 
          WRITE(6,*) "**** Boltzman: ", xboltzman 
          WRITE(6,*) "**** Temp: ", xtemp     
          WRITE(6,*) "**** Rho ini: ", rhoinitial
          WRITE(6,*) "**** SSD Coefficient: ", gammast
          WRITE(6,*) "***************************"
          writeflag(iphase+1) = 1
        END if
      END if