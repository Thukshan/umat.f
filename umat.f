*****************************************************************************
**  UMAT FOR ABAQUS/STANDARD INCORPORATING ELASTIC BEHAVIOUR  FOR PLANE    **
**  STRAIN AND AXI-SYMMETRIC ELEMENTS.                                     **
*****************************************************************************
*****************************************************************************
**
**
**
*USER SUBROUTINE
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
C
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      INTEGER, PARAMETER:: M=3,N=3,L0=24,L1=24,L2=12,KM=6,KN=6
      REAL*8,parameter  :: zero=1.0e-8
C
      DIMENSION DSTRESS(6), gauss(8,3), gausscoords(3,8),xI(3,3),
     + xnat(20,3),xnat8(8,3)

C     *** USER DEFINED ARRAYS ***
      REAL*8 :: stressM(3,3),plasStrainRate(3,3), totalStrainInc(6),
     + prod(M),tempNorm(M), tempDir(M), 
     + xRot(M,N),plasStrainInc2(6),Lp(3,3),
     + totplasstran(6),Le(3,3),
     + tSigma(6,6),tStran(6,6),tSigmainv(6,6),prod6(6,6),
     + xsnt(M,N),xsnv(KM),xnsv(KM),xsnnst(KM,KN),xIden6(KM,KN),
     + xnst(M,N),tmat(KM,KN),
     + trialstress(6),xjfai(KM,KN),
     + xjfaiinv(KM,KN),stressV(6),fai(6),trialstressM(3,3),
     + dstressinc(6),xu(3,3),xuinver(3,3),
     + xstressmdef(M,N),xstressdef(6),
     + xStiff(6,6),xStiffdef(6,6),elasspin(3,3),plasspin(3,3),
     + tSig(6,6),tStr(6,6),tSiginv(6,6), T(6,6),Tinv(6,6),
     + tStraninv(6,6),gmatinv(3,3),devstress(3,3),
     + compliance(6,6), print1(3,3),print2(3,3),print3(3,3),
     + print4(3,3),xrote(3,3),dtotstran(6),totstran(6),
     + totStrainRate(3,3),spin(3,3), tempstrain(3,3),curlfp(3,3),
     + curlfe(3,3),fp(3,3),fe(3,3),feinv(3,3),fpinv(3,3),xfinv(3,3),
     + tempNormGND(3), rhogndsys(24),
     + gv(9), rhoOutput(36),
     + tauR(12), sslip(3),nnorm(3),tnorm(3),burgX(M),
     + screw(3,3),edgeN(3,3),edgeT(3,3),scol(9),encol(9),etcol(9),
     + dstressth(6),dstranth(6), SVARS(9*8*2),
     + Fdot(3,3), invF(3,3), F(3,3), L(3,3),thermat(3,3),
     + expanse33(3,3), sigthv(6),sigthm(3,3), gmatinvnew(3,3), rho,
     + TAUCSSD, LpFeinv(3,3), matrix(3,3), update(3,3), gammast,rhoGND,
     + rhossdm, xhelmholtz, xfreq, gammazero,xboltzman,xtemp, rhoinitial
 
      integer debug, MZ, NZ, slipsys, iphase, crystalsys
      character(len=*),parameter :: fmt20 = "(' ',20(F6.4,1X))",
     + fmt3="(3(' ',(ES11.3,1X)))",fmt6 = "(' ',6(F12.3,1X))",
     + fmt5 = "(' ',5(F12.8,1X))",fmt7 = "(' ',7(F12.8,1X))",
     + fmt2 = "(24(' ',(I2,1X)))", fmt24 = "(24(' ',(ES11.3,1X)))"
        
      REAL*8,dimension(:,:),allocatable :: xNorm,xDir
      REAL*8,dimension(:),allocatable :: tau, gammadot, gndold,
     +  burgerv, tauc, tau2 ,S,gndtotal
      INTEGER,dimension(:),allocatable :: ids
      DOUBLE PRECISION,dimension(:,:),allocatable :: V,U,Sinv,A,AM,
     + Ainv,diagm,tempA, UT

      PARAMETER ( xgauss = 0.577350269189626)

#include <SMAAspUserSubroutines.hdr>

      include 'kmodelparams.f'
 
      debug = 0
      do while (debug == 1 .and. npt == 1 .and. noel == 1 .and. kinc ==1)
          ! For debugging, set debug to 1 and attach debugger
      end do
      
C     *** ZERO ARRAYS ***
      tempstrain=0.;spin=0.
       
      knsdv = nsdv

C     WRITE PROPS INTO STATEV TO INITIALIZE (ONCE ONLY),
      if (kinc == 1 .and. kstep==1) then

        do i=1,NSTATV
            STATEV(i) = 0. 
        end do

        do i=1,3
            do j=1,3
                STATEV(j+(i-1)*3) = props(j+1+((i-1)*3)) !gmat
            end do
        end do
         

      !fp     
      STATEV(81) = 1.0
      STATEV(85) = 1.0
      STATEV(89) = 1.0

      STATEV(54) = 0.01 ! initial sessile SSD density -> IS THIS CORRECT?

      STATEV(90) = 0.0
      STATEV(91) = 0.0
      STATEV(93) = 1.0
      Gcycle = 0.0
      svars = 0

      call MutexLock( 1 )    
      do i=1,9      
        kFp(noel,npt,i) = 0.0
      end do
      call MutexUnlock( 1 )
              
      call MutexLock( 2 )

      do i =1,3
        kgausscoords(noel,npt,i) = coords(i)
      end do
      call MutexUnlock( 2 )
      
      call MutexLock( 3 )
      do i=1,9      
        kcurlFp(noel,npt,i) = 0.0
      end do
      call MutexUnlock( 3 )
      
      end if    ! if loop started in line 91 ends here

C     Create identity matrix array
      xI=0.            
      DO I=1,3
        xI(I,I)=1.
      END DO

100   FORMAT (4(2X, E20.5))
150   FORMAT (4(2X, E20.5))

      F=0.  
      invF = 0. 
      Fdot = 0.
      L=0.

C     DETERMINE DEFORMATION AND VELOCITY GRADIENTS - CORRECTED ET 20/05/15     
      F = DFGRD0   
      Fdot = (DFGRD1-DFGRD0)/DTIME 
      CALL lapinverse(F,3,info,invF) 
   
C     CALL KMAT FOR MATERIAL BEHAVIOUR - Global stiffness matrix C
C     and stress 
      L = matmul(Fdot,invF)

   !  call KMLT(Fdot,invF,L)
      
 
      DO i=1,9
        statev(37+i) = kcurlfp(noel,npt,i)
      END DO
      
C ======================================================================
C                              BEGIN KMAT
C ======================================================================

      h = 0.0
      
      gndon = 0 !1=on, 0=off



C     *** ZERO ARRAYS ***

      result = 0.0; totstran = 0.0; totplasstran = 0.0;devstress=0.
      plasStrainInc2=0.
      xStiff=0.0; xStiffdef=0.0; DDSDDE=0.0; trialstressM=0.0
      tmat=0.0;xjfai=0.
      xjfaiinv=0.;plasStrainRate=0.;trialstress=0.; stressV=0.
      fai=0.;dstressinc=0.;xIden6=0.; xu=0.
      xuinver=0.;xRot=0.;prod=0.
      T=0; Tinv=0 !NEW ADDITIONS
      tSig=0.;tStr=0.;tSiginv=0.;tStraninv=0.
      tSigmainv=0.;gmatinv=0.;Lp = 0.;compliance=0.;Le = 0.
      print1=0.;print2=0.; print3 = 0.;print4 = 0.
      thermat =0.
      dstranth=0.;expanse33=0.;sigthv=0.;sigthm=0.
      xrote=0.;tempstrain=0.;spin=0.;totStrainRate=0.
      curlfp=0.;curlfe=0.;fp=0.;fe=0.;feinv=0.;fpinv=0.;xfinv=0.
      gmatinvnew = 0.
      ! 0. -> '.' makes 0 of type REAL
      
      DO I=1,KM; xIden6(I,I)=1.; END DO      
      DO I=1,M; xRot(I,I) = 1.; END DO
      
      
      !Zero allocate arrays from above!
      xNorm=0.; xDir=0.; tau=0.; gammadot=0.
      gndold=0.; ids=0; tau2=0.
      gndcas=0.;gndcap=0.;gndapy=0.;gndapr=0.;gndab=0.;rhognd=0.
      
      thermat(1,1) = alpha1; thermat(2,2)=alpha2; thermat(3,3) = alpha3
      
C     *** SET UP ELASTIC STIFFNESS MATRIX IN LATTICE SYSTEM ***   
      compliance(1,1:3) = (/1./E1,-v12/E1,-v13/E1/)
      compliance(2,2:3) =         (/1./E1,-v13/E1/)
      compliance(3,3:3) =                 (/1./E3/)
      compliance(4,4:4) =                       (/1./G12/)
      compliance(5,5:5) =                       (/1./G13/)
      compliance(6,6:6) =                       (/1./G13/)
C
      DO i=2,6
         DO j=1,i-1
            compliance(i,j)=compliance(j,i)
         END DO
      END DO
      
      CALL lapinverse(compliance,6,info,xStiff)
     
    
C     *** INITIALIZE USER ARRAYS ***

      DO i=1,3
        DO j=1,3
         gmatinv(i,j) = STATEV(j+(i-1)*3)
        END DO
      END DO
C
      p = STATEV(10) ! effective plastic strain 
C
      DO i=1,6
        totplasstran(i) = STATEV(10+i)
      END DO
C
      DO i=1,6
        totstran(i) = STATEV(16+i)
      END DO
C        
      DO i=1,6
        xstressdef(i) = STATEV(47+i)
      END DO
C
      rhognd = STATEV(37)    
      rhossd = STATEV(54)
      
      r = STATEV(56)
                
      DO i=1,nSys
        gndold(i) = STATEV(56+i)
      END DO
C
      DO i=1,3
        DO j=1,3
          fp(i,j) = STATEV(80+j+((i-1)*3))
        END DO
      END DO              
C
      DO i=1,3
        DO j=1,3 
         curlfp(i,j) = STATEV(37+j+(i-1)*3)
        END DO
      END DO
       

      Nf = STATEV(91)

C     *** DIRECTIONS FROM LATTICE TO DEFORMED SYSTEM ***
        
      CALL kdirns(gmatinv,crystalsys,nSys,xDir,xNorm)
        

C     *** STIFFNESS FROM LATTICE TO DEFORMED SYSTEM ***

      CALL rotord4sig(gmatinv,tSig)
      CALL rotord4str(gmatinv,tStr)
      CALL lapinverse(tSig,6,info2,tSiginv)

      prod6 = matmul(tSiginv,xStiff)      
      xStiffdef = matmul(prod6,tStr)

      expanse33 = matmul(matmul(gmatinv,thermat),transpose(gmatinv))
      expanse33 = expanse33*DTEMP !dstrain = alpha*dT
      
      CALL kmatvec6(expanse33,dstranth)
      dstranth(4:6) = 2.0*dstranth(4:6)
            

    
C     *** DETERMINE INCREMENT IN TOTAL STRAIN (6X1 ***     

      tempstrain=(L+transpose(L))*0.5*dtime
      spin=(L-transpose(L))*0.5 

      CALL kmatvec6(tempstrain,dtotstran)
      dtotstran(4:6) = 2.0*dtotstran(4:6)


C     *** COMPUTE TRIAL STRESS ***

      stressV = xstressdef ! old stress
      trialstress = stressV+matmul(xStiffdef,dtotstran)-
     + matmul(xStiffdef,dstranth)            
      CALL kvecmat6(trialstress,trialstressM) 
            
      CALL kvecmat6(stressV,stressM) 
      trialstressM = trialstressM + (matmul(spin,stressM) - 
     + matmul(stressM,spin))*dtime 


 
C     *** CALCULATE RESOLVED SHEAR STRESS ON A SLIP SYSTEMS ***


      DO I=1,nSys
        tempNorm = xNorm(I,:); tempDir = xDir(I,:)
        prod = matmul(trialstressM,tempNorm)
        tau(I)= dot_product(prod,tempDir)
        IF(tau(I) .LT. 0.0) THEN
          tau(I) = -1.E0*tau(I) !ensures tau(I) is +ve
          DO K=1,3
            xDir(I,K)=-1.E0*xDir(I,K)
          END DO
        END IF
      END DO
        
          
      xtau = maxval(tau/tauc)  


C     *** PLASTIC DEFORMATION ***

      IF (xtau >= 1.0 ) THEN
      
      faivalue=1.
      xacc=1.e-8
      iter=0


C     *** USE NEWTON METHOD TO DETERMINE STRESS INCREMENT ***

      DO WHILE (faivalue .gt. xacc)      
      
      iter=iter+1
      ! For each iteration do everything from lines 346-424. 

      !============================================================================   
      !  Slip rule:
      !  Returns Lp and tmat required to define the material jacobian.  
      !============================================================================  

      CALL kslip6(xNorm,xDir,tau,tauc,burgerv,caratio,      
     + dtime,nSys,r,Lp,tmat,TEMP,rhossdm, xhelmholtz,
     + xfreq, gammazero, xboltzman, xtemp, rhoinitial)
      
      if (any(tmat /= tmat) .or. any(tmat-1 == tmat)) then                                                                                                                                                            
          ! The "sinh" has probably blown up  -- then try again with smaller dt                                                                                                                                      
          pnewdt = 0.5                                                                                                                                                                                               
          write(*,*) "W! The tmat is nan or inf (NOEL, NPT, time): ",                                                                                                                                              
     + NOEL, NPT, time                                                                                                                                                                                             
          return                                                                                                                                                                                                     
      end if
      
      !============================================================================  


C     *** DETERMINE PLASTIC STRAIN INCREMENETS FOR UEL

      plasStrainRate = (Lp+transpose(Lp))*0.5*dtime
      CALL kmatvec6(plasStrainRate,plasStrainInc2)
      plasStrainInc2(4:6) = 2.0*plasStrainInc2(4:6)            



C     *** CALCULATE THE STRESS INCREMENT ***

      xjfai =  xIden6 + matmul(xStiffdef,tmat)
      CALL lapinverse(xjfai,6,info3,xjfaiinv)
      fai = trialstress - stressV - matmul(xStiffdef,plasStrainInc2)
      dstressinc = matmul(xjfaiinv,fai)
      stressV = stressV + dstressinc
      CALL kvecmat6(stressV,stressM)      
      faivalue = sqrt(sum(fai*fai))        


C     *** UPDATE RESOLVED SHEAR STRESS ACCORDING TO NEW STRESS ***
       
      DO I=1,nSys
      
          tempNorm = xNorm(I,:); tempDir = xDir(I,:)    
          prod = matmul(stressM,tempNorm)
          tau(I)= dot_product(prod,tempDir)

          IF(tau(I) < 0.0) THEN
            tau(I) = -1.E0*tau(I)
            DO K=1,3
              xDir(I,K)=-1.E0*xDir(I,K)
            END DO
          END IF
          !=============================== ids
          IF(tau(i)/tauc(i) >= 1.0) THEN
            ids(i) = i; tau2(i) = tau(i)/tauc(i)
          END IF
          !===============================
          
      END DO
        
      xtau = maxval(tau/tauc) 
          
      
      IF (iter .gt. 50) THEN !i.e. if more than 50 iterations are required:
          WRITE(*,*) "WARNING NEWTON LOOP NOT CONVERGED: NOEL, NPT, 
     +     time:", NOEL, NPT, time
          pnewdt = 0.5
          return
          !CALL XIT 
      END IF
                 
      !*** THE END OF NEWTON ITERATION ***
      END DO


C     *** NOW CALCULATE THE JACOBIAN***

      DDSDDE = matmul(xjfaiinv,xStiffdef)
       
       
C     *** ROTATE STRESS BACK TO GLOBAL SYSTEM *** 

      xstressmdef = stressM


C     *** UPDATE OUTPUT VARIABLES ***

      plasStrainrate=(Lp+transpose(Lp))*0.5       
      pdot=sqrt(2./3.*sum(plasStrainrate*plasStrainrate))
      
      p = p + pdot*dtime
      r = r + h*pdot*dtime
      
C     *** UPDATE PLASTIC DEFORMATION GRADIENT
    
      print2 = 0.; print3 = 0.
      print2 = xI - Lp*dtime      
      CALL kdeter(print2,deter)      
      IF (deter /= 0.0) THEN
         CALL lapinverse(print2,3,info4,print3)
         fp = matmul(print3,fp)
      ELSE
         fp = fp
      END IF  
      
      

      
    !=========================================================================
    ! SSD Evolution 
     
      rhossd = rhossd + (gammast*pdot*dtime)
      
    !=========================================================================
               
C     *** ELASTIC DEFORMATION ***     
      ELSE
      xstressmdef = trialstressM
      DDSDDE = xStiffdef      
      END IF
      
      CALL kmatvec6(xstressmdef,xstressdef) !output stress
      devstress = xstressmdef - 1./3.*trace(xstressmdef)*xI
      vms = sqrt(3./2.*(sum(devstress*devstress))) !von mises stress 


    !=========================================================================
    ! *** DETERMINE DENSITY OF GNDs
    ! Update P.Ashton November 2015
    !=========================================================================
      
      IF(MAXVAL(ABS(curlfp)) <= 1.0e-8) THEN 
        rhoGND = 0.;gndab  = 0.;gndapr = 0.;gndapy = 0.
      
      ELSE

      IF (gndon == 0) THEN !Switching GND evolution on and off
        rhoGND = 0.;gndab  = 0.;gndapr = 0.;gndapy = 0.
      
      ELSE   
      
      MZ=9;NZ=36
      
      ALLOCATE(A(mz,nz),AM(mz,nz),V(nz,nz),S(mz),U(mz,mz),Sinv(nz,mz),
     + Ainv(nz,mz),tempA(nz,mz),UT(mz,mz),STAT=ialloc)
     
      A=0.;AM=0.;V=0.;S=0.;U=0.;Sinv=0.;Ainv=0.;tempA=0.;UT=0.
      
      gv = reshape(curlfp,(/9/))   
      

           
      tauR = tau2(1:12)
        
      DO i=1,12
        IF(tauR(i) > 1) THEN
           tauR(i) = tauR(i) - 1
        END IF 
      END DO 
      
      tausum = sum(tauR)
      tauBas = sum(tauR(1:3))
      tauPris= sum(tauR(4:6))
      tauPyr = sum(tauR(7:12))
       
  !    DO slipSys=1,2 ! Cycle through each system
      
  !    IF (ids(slipSys) > 0) THEN
   
        DO I=1,12 ! 1-12 for <a> type
        
            burgX = xDir(I,:)
            burgX = burgX*burger1 
         
            sslip = xDir(I,:); nnorm = xNorm(I,:) 
            CALL CrossProd(sslip,nnorm,tnorm) 
         
            CALL DyadicProd(sslip,burgX,screw)
            CALL DyadicProd(nnorm,burgX,edgeN)
            CALL DyadicProd(tnorm,burgX,edgeT)
         
            CALL Convert2Col(screw,scol)
            CALL Convert2Col(edgeN,encol)
            CALL Convert2Col(edgeT,etcol)
         
            DO J=1,9
                A(J,I)    = scol(J)
                A(J,I+12) = encol(J)
                A(J,I+24) = etcol(J)
            END DO
                
        END DO   ! Generate current A-matrix loop
      


  ! ***********************************************************************    
      ! Matrix inversion by singular value decomposition: A+ = [V][S+][U*]
  ! ***********************************************************************   
        CALL SVD(A,U,S,V,MZ,NZ) 
      
        DO i = 1, ubound(S,1)
            IF (S(i) > 1e-6) THEN
                Sinv(i,i)= 1.0 / S(i)
            END IF
        END DO 
        
            
        UT=transpose(U)
      
        NCOL=NZ;MX=NZ;NX=MZ 
        CALL KMLTM (V,Sinv,tempA,NCOL,MX,NX) 

      
        NCOL=MZ;MX=NZ;NX=MZ 
        CALL KMLTM (tempA,UT,Ainv,NCOL,MX,NX) 

        CALL KMLT36991(Ainv,gv,rhoOutput)
        
 
        rhos = sum(rhoOutput(1:12))
        rhoen = sum(rhoOutput(13:24))
        rhoet = sum(rhoOutput(25:36))
        rhofinal = sqrt((rhos*rhos) +(rhoen*rhoen) +(rhoet*rhoet))

     
      rhoGND = rhofinal      



      

      END IF ! End of GND on/off switch (line 488) 
      
      END IF ! if curl(Fp) < 1e-8 (line 483)


C     *** ORIENTATION UPDATE ***
      !Assuming that all rigid body rotatoin is lumped into Fe and that the elastic strians are small 
!     then the elastic spin is We = d(Fe)/dt inv(Fe)
      !L = We + Fe Lp inv(Fe) therefore 
      !We = L - Fe Lp inv(Fe)
      ! G(t+dt) = G(t) + We G(t)dt dt or an implicit update is G(t+dt)  = G(t)exp[We(t+dt)dt]  ~ inv[I - We(t+dt) dt] G(t) 
      
      ! We need Fe and inv(Fe) using F = Fe Fp gives Fe = F.inv(Fp)
      CALL kdeter(Fp,deter)      
      
      IF (deter /= 0.) THEN
         Fpinv = 0.
         CALL lapinverse(Fp,3,info5,Fpinv)
         Fe = matmul(F,Fpinv)          
      ELSE
         write(*,*) "Error in orientation update: finding inv(Fp)",noel,
     +    npt, kinc
         call XIT 
      
      END IF  
      
      
      CALL kdeter(Fe,deter)      
      
      IF (deter /= 0.) THEN
         Feinv = 0.
         CALL lapinverse(Fe,3,info5,Feinv)       
      ELSE
          write(*,*) "Error in orientation update: finding inv(Fe)",noel
     +     ,npt, kinc
         call XIT 
      
      END IF        
            
      LpFeinv = 0.; 
      LpFeinv = matmul(Lp, Feinv)
      Le = L - matmul(Fe,LpFeinv)        
      elasspin=(Le-transpose(Le))*0.5
      matrix = xI - elasspin*dtime      
      CALL kdeter(matrix,deter)      
      
           
          
      
      
      IF (deter /= 0.) THEN
         update = 0.
         CALL lapinverse(matrix,3,info5,update)
         IF(info5 /= 0) write(*,*) "inverse failure: print3 in kmat"
         gmatinvnew = matmul(update,gmatinv)
      ELSE         
         gmatinvnew = gmatinv
      write(*,*) "WARNING gmatinv not updated at noel,npt, kinc:", noel,
     + npt, kinc
      END IF      

      gmatinv =  gmatinvnew            
      
       if (maxval(gmatinv) > 1) then
          write(*,*) "something very wrong with gmatinv"
          call XIT
       end if

C     ********************* UPDATE STATE VARIABLES *********************
      DO i=1,3
        DO j=1,3
        STATEV(j+(i-1)*3) = gmatinv(i,j)
        END DO
      END DO

      STATEV(10) = p
 
      DO i=1,6
        STATEV(10+i) = totplasstran(i) + 
     +                                  plasStrainInc2(i) 
      END DO

      DO i=1,6
        STATEV(16+i) = totstran(i) + dtotstran(i)
      END DO

      !23-25 are rotations stored in UEL
      STATEV(26) = rhognd !all
      STATEV(27) = gndab  !a basal
      STATEV(28) = gndapr !a prismatic  
      STATEV(29) = gndapy !a pyramidal
      STATEV(30) = gndcap !c+a primary    
      STATEV(31) = gndcas !c+a secondary
      STATEV(32) = maxval(plasStrainrate)
      STATEV(33) = pdot
      STATEV(34) = xtau !XTAUC1 
      !35-37 are free
      !38-46 are curlfp terms calulated in gradient routine (kcurl)

      DO i=1,6
       STATEV(47+i) = xstressdef(i)
      END DO

      STATEV(54) = rhossd !ssd = mobile (glissile) dislocations; gnd = sessile dislocations
      STATEV(55) = vms
      STATEV(56) = r !isotropic hardening function (h = constant) 

      DO i=1,nSys !GNDs on indiviual systems - Only room for 24 atm
       STATEV(56+i) = gndold(i)
      END DO

      DO i=1,3
       DO j=1,3 
        STATEV(80+j+(i-1)*3) = fp(i,j)
       END DO
      END DO   

      ! Stored energy density and related calculations - A. Bergsmo 2019
      Gdot=0.0
      

      DO i = 1, 6
          Gdot = Gdot + ABS(0.05 * xstressdef(i) * plasStrainInc2(i))
      END DO

C     Combine SSDs and GNDs      
      DislDens = rhognd + rhossd
      DislDens = ABS(SQRT(DislDens))
      
      if (DislDens == 0) then
          write(*,*) "Dislocation Density Wrong!"
          DislDens = 0.01
      end if
      
      Gdot = Gdot / DislDens

      Gstored = STATEV(90) + Gdot

      Gcrit = 404.0 !J/m^-2 B. Chen 2019 RR1000
      cycletime = 22.0 !s TIME OF FATIGUE CYCLE

      Gcycle = STATEV(92) + Gdot
      cyc = STATEV(93)
       
      if(TIME(2) > cyc*cycletime) then
          Nf = Gcrit/Gcycle ! TIME TO INITIATION
          !write(*,*) "Fatigue life is:", Nf, TIME(2)
          Gcycle = 0
          cyc = cyc+1
      end if

      STATEV(90) = Gstored
      STATEV(91) = Nf
      STATEV(92) = Gcycle
      STATEV(93) = cyc
C ======================================================================
C                               END KMAT
C ======================================================================
   
C     RECOVER stress for calculating residual force
      DO K=1,6
        stress(K)=STATEV(47+K)
      END DO
      
      call MutexLock( 1 )
        DO i=1,9                                                      
            kFp(noel,npt,i)= statev(80+i)
        END DO
      call MutexUnlock( 1 )
      
      IF (npt == 8 ) THEN ! update curl Fp

C     SPECIFY GAUSS POINT LOCAL COORDS, AND WEIGHTING CONSTANTS 
      
        INCLUDE 'kgauss.f'     
        xnat8 = xnat(1:8,:)       
        nsvars = 9*8*2 ! 8 x 89 = 712       
         
          
        DO kint =1,8
            DO i=1,3         
                gausscoords(i,kint) = kgausscoords(noel,kint,i)                          
            END DO
    
            DO i=1,9          
                svars(72 + i + 9*(kint-1)) = kFp(noel,kint,i)         
            END DO
        END DO
   
         
C     A FULL INTEGRATION GRADIENT SCHEME

        CALL kcurl(nsvars,svars,9,xnat8,gauss,gausscoords)
      
        call MutexLock( 3 )
        DO kint =1, 8
            DO i=1, 9
                kcurlFp(noel,kint,i) = svars(i + 9*(kint-1))
            END DO
        END DO
        call MutexUnlock( 3 )
    
            
      END IF
          
      RETURN
      END

      include 'uexternaldb.f' ! Sets up Mutexes - Abaqus Subroutine
      include 'kdirns.f'      ! Loads the slip dirs and norms
      include 'kslip6.f'      ! Slip rule
      include 'ksvd2.f'       ! Lapack library used in GND calcs
      include 'kcurl.f'       ! Curl calculation for GNDs
      include 'kshapes.f'     ! Shape functions
      include 'utils.f'       ! Utility subroutines

            
            