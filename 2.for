*
*     ***************************************************************************     *
*     ***************************************************************************     *
*
*                         HETVAL FOR THE CEMENT HYDRATION 
*
*     ***************************************************************************     *
*     ***************************************************************************     *
*
      SUBROUTINE HETVAL(CMNAME,TEMP,TIME,DTIME,STATEV,FLUX,
     1 PREDEF,DPRED)
*      
      INCLUDE 'ABA_PARAM.INC'
*
      CHARACTER*80 CMNAME
*
      DIMENSION TEMP(2),STATEV(20),PREDEF(*),TIME(2),FLUX(2),
     1 DPRED(*) 
*      
*     ---------------------------------------------------------------------------     *
*     Parameter input
*
      Cemass   = 437.0d0           !unit: kg/m^3,   cement content in one cubic meter conconcrete
      C_Q0     = 431000D0          !unit: J/kg,   maximum heat for one kg cement
*
      tau  = 14.0d0
      beta = 0.94d0
      
      tolerance = 0.0001d0


*     ---------------------------------------------------------------------------     *
      
      Time_eq       = STATEV(3)
      Time_eq_coeff = STATEV(4) 
      
      if(TIME(2) .le. tolerance) then
          hydro = 0.0d0
          
          dq_dt = 0.0d0
      else
          hydro = STATEV(2)
*           
          term_1  = Cemass * C_Q0
          term_2  = hydro * (tau / Time_eq) ** beta * (beta / Time_eq )
*          
          dq_dt   =  term_1 * term_2 * Time_eq_coeff
          
      end if
      
     
      
*
*     Q(t)     = Q0 * (1 - exp(-b*t^c))                                ! t is time in day
*     dQ(t)_dt = b*c*Q0 * t^(c-1) * exp(-b*t^c)                        ! t is time in day
*      
*     dQ(t)_dt = b*c*Q0 * (t/24)^(c-1) * exp(-b*(t/24)^c)   / 24       ! t is time in hour
* 
*
*      term_1  = C_Q0*b*c/24*(Time(2)/24)**(c-1)
*      term_2  = exp(-b*(Time(2)/24)**c) 
*
*      dq_dt   = Cemass* term_1 *  term_2
*
      FLUX(1) = dq_dt
      
      
*     WRITE(6,*) 'aa ', PREDEF(1)

      STATEV(1) = TEMP(1)        ! TEMP(1) - current temperature; TEMP(2) - temperature increment  
*     ---------------------------------------------------------------------------     *


      
      CALL AMB_TEMP( TIME(2), TEMP_ENV)

      STATEV(5) = TEMP_ENV 
      
      RETURN
      END


      
      
      
*
*     ***************************************************************************     *
*     ***************************************************************************     *
*
*                         subroutine: DFLUX for the heat transfer
*
*     ***************************************************************************     *
*     ***************************************************************************     *
*       
      
       SUBROUTINE DFLUX(FLUX,SOL,JSTEP,JINC,TIME,NOEL,NPT,COORDS,JLTYP,
     1 TEMP,PRESS,SNAME)
* 
      INCLUDE 'ABA_PARAM.INC'
*     ---------------------------------------------------------------------------     *

      DIMENSION COORDS(3),FLUX(2),TIME(2)
      CHARACTER*80 SNAME
*     
      WIND_SPEED = 1.0D0
      
      COEFF = 24.3D0 + 14.5D0 * WIND_SPEED      ! KJ/(m2 h C)
      
      
      
      COEFF = COEFF * 1000.0d0                  ! J/(m2 h C) 
      
      CALL AMB_TEMP( TIME(2), TEMP_ENV)
      
*      TEMP_ENV  =  30.0D0
      
      FLUX(1) = -COEFF * (SOL - TEMP_ENV)
      FLUX(2) = -COEFF
          
	  
*	write(6,*)  'TEMP_ENV = ', TEMP_ENV
*      write(6,*)  'NPT = ', NPT
*      write(6,*)  'SOL = ', SOL
*       write(6,*)  'FLUX(1) = ',FLUX(1)
*       write(6,*)  'FLUX(2) = ',FLUX(2)
      
      RETURN
      END
      
      
      
*
*     ***************************************************************************     *
*     ***************************************************************************     *
*
*                         subroutine: film for ambient condition
*
*     ***************************************************************************     *
*     ***************************************************************************     *
*      
      subroutine film(h,sink,temp,jstep,jinc,time,noel,npt,coords,
     1     jltyp,field,nfield,sname,jusernode,area)
*
      include 'aba_param.inc'
*
      dimension coords(3),time(2),field(nfield),h(2)
      character*80 sname
*     ---------------------------------------------------------------------------     *

      h(1) = 0.d0
      h(2) = 0.d0
      sink = 0.d0
      if (jltyp .eq. 0) then
         h0 = 0.91d0
         h1 = 0.0d0
         h(1) = h0 + h1 * temp
         h(2) = h1
*         sink = 100.0d0 * (1.d0 + time(1) / 3600.d0)
         sink = 0.0d0
      end if
*     ---------------------------------------------------------------------------     *
      RETURN
      END
      
      
      
*
*     ***************************************************************************     *
*     ***************************************************************************     *
*
*                         subroutine: USDFLD to define hydration degree
*
*     ***************************************************************************     *
*     ***************************************************************************     *
*          
      
      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
*
      INCLUDE 'ABA_PARAM.INC'
*
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
*     ---------------------------------------------------------------------------     *


*      WRITE(6,*) 'aa ', 1 
      

      w_c  = 0.42 d0
      hydro_u  = 1.031d0*w_c / (0.194d0 + w_c)
      tau  = 14.0d0
      beta = 0.94d0
      
      act_energy = 41841.0d0      !E_a,  (J mole-1)  
      gas_con    = 8.314d0        !universal gas constant,  (J K-1 mol-1) 
      
      tolerance = 0.0001d0
*     ---------------------------------------------------------------------------     *      
      Temp = STATEV(1)
      
      Temp_c = 273.0d0 + Temp  !Temp
      
      Temp_r = 273.0d0 + 20.0d0 
*     ---------------------------------------------------------------------------     *
      
      Time_eq = STATEV(3)
      
      Time_eq_coeff = exp( act_energy / gas_con * (1.0d0 / Temp_r - 1.0d0 / Temp_c ) )
      
      Time_eq = Time_eq + DTIME * Time_eq_coeff
      
      
      
      if ( Time_eq  .le. tolerance) then
          hydro = 0.0d0
      else
          hydro =  hydro_u * exp (-1.0d0 * (tau / Time_eq )**(beta)  )
      end if
      
      
      
      STATEV(6) = hydro - STATEV(2)
      
      
      FIELD(1) = hydro
      
      STATEV(2) = hydro
      
      STATEV(3) = Time_eq
      STATEV(4) = Time_eq_coeff
      
      
      

*     ---------------------------------------------------------------------------     *
      
      RETURN
      END




      SUBROUTINE UEXPAN(EXPAN,DEXPANDT,TEMP,TIME,DTIME,PREDEF,DPRED,
     $	              STATEV,CMNAME,NSTATV,NOEL)
      
      INCLUDE 'ABA_PARAM.INC'
*     TIME OF CURING
	PARAMETER(SHRINKAGE = 7.0D0)

	CHARACTER*80 CMNAME

	DIMENSION EXPAN(*),DEXPANDT(*),TEMP(2),TIME(2),PREDEF(*),
     $          DPRED(*),STATEV(NSTATV)
      
      DIMENSION CREEP_BASIC(6)

*
*     ----------------------------------------------------------------------------------      *      
*     1.0 Parameters
*     ----------------------------------------------------------------------------------      *      
*


*     ULTIMATE SHRINKAGE
*     ACI RECOMMENDATION EQ. A-2 
*
      SHU_ULTIMATE0 = -520.0D-6  
*
*     INPUT SLAB THICKNESS (TRUE THICKNSS = 2V/S), IN INCH. 
*
      
      w_c  = 0.42 d0 
      hydro_u  = 1.031d0*w_c / (0.194d0 + w_c)
      
      hydro_incre = STATEV(6)
      
      DELTA_SHRINKAGE =  SHU_ULTIMATE0 * (hydro_incre / hydro_u) /2.0d0 
      
         
      EXPAN(1) = TEMP(2) * 1.0d-5   + DELTA_SHRINKAGE
      EXPAN(2) = TEMP(2) * 1.0d-5   + DELTA_SHRINKAGE
      EXPAN(3) = TEMP(2) * 1.0d-5   + DELTA_SHRINKAGE
            
      STATEV(11) = STATEV(11) +    DELTA_SHRINKAGE 
          
      STATEV(12) = STATEV(12) +    EXPAN(1)
      
      STATEV(7)  = hydro_incre + STATEV(7)

*              
*
*     ----------------------------------------------------------------------------------      *


	RETURN
	END


*
*     ***************************************************************************     *
*     ***************************************************************************     *
*
*                         subroutine: USDFLD to define hydration degree
*
*     ***************************************************************************     *
*     ***************************************************************************     *
*  
      SUBROUTINE AMB_TEMP( TIME, TEMP_ENV)
      
      INCLUDE 'ABA_PARAM.INC'

      if ( (Time .ge.	0.0	d0)     .and. (Time .lt.0.03d0))	    TEMP_ENV  = 	33.00d0
      if ( (Time .ge.	0.03	d0) .and. (Time .lt.1.32667d0))	    TEMP_ENV  = 	32.80d0
      if ( (Time .ge.	1.32667	d0) .and. (Time .lt.2.41139d0))	    TEMP_ENV  = 	32.30d0
      if ( (Time .ge.	2.41139	d0) .and. (Time .lt.3.19528d0))	    TEMP_ENV  = 	30.30d0
      if ( (Time .ge.	3.19528	d0) .and. (Time .lt.3.95611d0))	    TEMP_ENV  = 	28.30d0
      if ( (Time .ge.	3.95611	d0) .and. (Time .lt.4.9125d0))      TEMP_ENV  = 	25.40d0
      if ( (Time .ge.	4.9125	d0) .and. (Time .lt.5.85778d0))	    TEMP_ENV  = 	25.30d0
      if ( (Time .ge.	5.85778	d0) .and. (Time .lt.6.94722d0))	    TEMP_ENV  = 	24.00d0
      if ( (Time .ge.	6.94722	d0) .and. (Time .lt.7.93222d0))	    TEMP_ENV  = 	23.70d0
      if ( (Time .ge.	7.93222	d0) .and. (Time .lt.8.49611d0))	    TEMP_ENV  = 	23.70d0
      if ( (Time .ge.	8.49611	d0) .and. (Time .lt.9.26389d0))	    TEMP_ENV  = 	23.10d0
      if ( (Time .ge.	9.26389	d0) .and. (Time .lt.10.32611d0))    TEMP_ENV  = 	23.40d0
      if ( (Time .ge.	10.32611d0) .and. (Time .lt.11.38639d0))	TEMP_ENV  = 	22.50d0
      if ( (Time .ge.	11.38639d0) .and. (Time .lt.12.39333d0))	TEMP_ENV  = 	22.30d0
      if ( (Time .ge.	12.39333d0) .and. (Time .lt.13.37278d0))	TEMP_ENV  = 	23.50d0
      if ( (Time .ge.	13.37278d0) .and. (Time .lt.14.27278d0))	TEMP_ENV  = 	23.40d0
      if ( (Time .ge.	14.27278d0) .and. (Time .lt.15.43833d0))	TEMP_ENV  = 	23.40d0
      if ( (Time .ge.	15.43833d0) .and. (Time .lt.16.25861d0))	TEMP_ENV  = 	23.70d0
      if ( (Time .ge.	16.25861d0) .and. (Time .lt.17.29611d0))	TEMP_ENV  = 	26.40d0
      if ( (Time .ge.	17.29611d0) .and. (Time .lt.18.25139d0))	TEMP_ENV  = 	29.30d0
      if ( (Time .ge.	18.25139d0) .and. (Time .lt.19.29167d0))	TEMP_ENV  = 	32.50d0
      if ( (Time .ge.	19.29167d0) .and. (Time .lt.20.48194d0))	TEMP_ENV  = 	30.20d0
      if ( (Time .ge.	20.48194d0) .and. (Time .lt.21.54167d0))	TEMP_ENV  = 	35.00d0
      if ( (Time .ge.	21.54167d0) .and. (Time .lt.22.67111d0))	TEMP_ENV  = 	38.00d0
      if ( (Time .ge.	22.67111d0) .and. (Time .lt.23.55972d0))	TEMP_ENV  = 	36.70d0
      if ( (Time .ge.	23.55972d0) .and. (Time .lt.24.38139d0))	TEMP_ENV  = 	35.05d0
      if ( (Time .ge.	24.38139d0) .and. (Time .lt.25.41444d0))	TEMP_ENV  = 	33.40d0
      if ( (Time .ge.	25.41444d0) .and. (Time .lt.26.59639d0))	TEMP_ENV  = 	30.40d0
      if ( (Time .ge.	26.59639d0) .and. (Time .lt.28.07306d0))	TEMP_ENV  = 	30.10d0
      if ( (Time .ge.	28.07306d0) .and. (Time .lt.30.30556d0))	TEMP_ENV  = 	27.00d0
      if ( (Time .ge.	30.30556d0) .and. (Time .lt.32.4425d0))	    TEMP_ENV  = 	24.90d0
      if ( (Time .ge.	32.4425	d0) .and. (Time .lt.33.31472d0))	TEMP_ENV  = 	24.40d0
      if ( (Time .ge.	33.31472d0) .and. (Time .lt.35.33306d0))	TEMP_ENV  = 	25.10d0
      if ( (Time .ge.	35.33306d0) .and. (Time .lt.37.43417d0))	TEMP_ENV  = 	24.20d0
      if ( (Time .ge.	37.43417d0) .and. (Time .lt.39.69944d0))	TEMP_ENV  = 	24.80d0
      if ( (Time .ge.	39.69944d0) .and. (Time .lt.42.23556d0))	TEMP_ENV  = 	31.20d0
      if (  Time .ge.	42.23556d0) 		                        TEMP_ENV  = 	31.20d0
         
         
      RETURN
      END