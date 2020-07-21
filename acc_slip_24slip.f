    *** USER DEFINED ARRAYS ***
	 gammaslip(24), gammaslip2(24)

	*** INITIALIZE USER ARRAYS ***
	 !line 242:
	 STATEV(129)=STATE(12) 
	
	  DO i=1, nSys
	      gammaslip(i) = STATEV(102+i)
	  END DO
	 
	  STATEV(127)=STATEV(112)	!Abaqus doesn't recognize SDVs ending in 2
	  STATEV(128)=STATEV(122)
	  
	  DO i=1,nSys
	      gammaslip2(i) = STATEV(129+i)
	  END DO 

      STATEV(154)=STATEV(132)
	  STATEV(155)=STATEV(142)
	  STATEV(156)=STATEV(152) 
	  
	*** UPDATE STATE VARIABLES ***
	  DO i=1, nSys
	      STATEV(102+i) = gammaslip(i) 
	  END DO
	  
	  STATEV(127)=STATEV(112)	!Abaqus doesn't recognize SDVs ending in 2
	  STATEV(128)=STATEV(122)
	  
      ! accumulated slip on individual slip systems
	  DO i=1,nSys
          STATEV(129+i) = gammaslip2(i) !Change if using acc_slip_24slip code
      END DO
	  
	  STATEV(154)=STATEV(132)
	  STATEV(155)=STATEV(142)
	  STATEV(156)=STATEV(152) 