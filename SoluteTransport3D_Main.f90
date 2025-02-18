!***********************Program LPRWSTM***********************************
!
!	This program solves the advection-dispersion equations to
!   calculates the concentration in a 3D Aquifer due to a release 
!	of a instantaneous or continuous contaminant source using
!  Lagrangian Parcel Random Walk Method (LPRWSTM)
!   
!   Author: Uditha C. Bandara, PhD. P.E. 
!   Senior Water Resource Engineer, 
!   Sustainable Groundwater Management Office,
!   State of California, Department of Water Resources 
! *****************************************************************************

    
    PROGRAM AdvDiff3DLagrPar
    USE AdvDispLagrPar_Mod
    
    IMPLICIT NONE
    
   INTEGER ix, iy, np, m,  nnode, kx, ky, Iout, imon, iyr
   REAL (KIND=IKIND)  Time, Vtot,Term1, rr, WW !, BB
   REAL (KIND=IKIND) A1, A2, B2, ZZ, CN, C0, start_time, end_time
   LOGICAL        :: Current_STP

   CALL CPU_TIME(start_time)

   Iout=0
  
   !--READ TIME INVARIANT DATA------------------------------------------------------------------------------------------------------------------!
   
    CALL ReadInput ( )   
    CALL Initialize ( )   
  ! CALL ReadIbound ( )   ! To be completed
    IF (SteadyStateVel) CALL  ReadSSVelocity ( ) 
    CALL InitCond ( )
  

!--START THE TIME LOOP----------------------------------------------------------------------------------------------------------------------------!
  tot_mass_rel=0.00
  imon=1
  iyr=0
 Current_STP=.FALSE.
 
    DO  60 jt=1,num_time_steps
        
        Time=jt*del_t       ! in days (del_t should be in days)      
        CALL  Print_Time( jt, Time)     
        CALL Get_Year(imon,Time,Current_STP, iyr )

         !--Read Transient Velocity Data
        IF ( NOT (Current_STP)  .and.  NOT(SteadyStateVel) ) THEN
            CALL ReadTransVelocity(Current_STP ,Time) 
        END IF
        
       !--Release Parcels Corresponding to Discharges 
       IF ( ( Time<=t_spill) )  CALL Release_LP (iflg_vel, Time )

       CALL Advect_LP (iflg_vel) 

       !--Output Results
       IF ( nt_skip==1) THEN 
           CALL Output ( Iout, Time) 
       ELSE IF  (MOD(jt,nt_skip-1) == 0 ) THEN
           CALL Output ( Iout, Time)
       ENDIF 

60  CONTINUE	
!--END OF THE TIME LOOP--------------------------------------------------------------------------------------------------------------------------!
    
   CALL CPU_TIME(end_time)
   
   WRITE(*,*) "Elapsed Time:", end_time-start_time
   WRITE(fid_dbug,*) "Elapsed Time:", end_time-start_time
   
    STOP
    END PROGRAM AdvDiff3DLagrPar
    
    
    
    
    
    
    !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
      
  SUBROUTINE  ReadSSVelocity ( ) 
      USE AdvDispLagrPar_Mod

       !!--Read Steady State Velocity
            IF (iflg_vel==1) THEN
                    CALL Get_Cell_Velocity_SS ( )  ! Cell center velocity from MODFLOW
            ELSE 
                    CALL Get_Nodal_Velocity_SS ( ) 
            ENDIF
  END SUBROUTINE  ReadSSVelocity
  !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
  
  
  !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
  SUBROUTINE ReadTransVelocity (Cur_STP , aTime)
      USE AdvDispLagrPar_Mod

  REAL (KIND=IKIND), INTENT (IN)   :: aTime
  LOGICAL, INTENT (INOUT)                  :: Cur_STP
       !--Read Velocity for monthly intervals (IWFM outputs monthly average velocities)
                IF (iflg_vel==1) THEN
                        !CALL Get_Cell_Velocity (Time )  ! Cell center velocity from MODFLOW-not done yet
                ELSE 
                        CALL Get_Nodal_Velocity (aTime )  ! For IWFM    
                ENDIF      
                Cur_STP=.TRUE.
  END  SUBROUTINE ReadTransVelocity
  !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
 
   !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!      
  SUBROUTINE Release_LP (jflg_vel, aTime) 
       !--Particle Release for the discharg period
  !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
      USE AdvDispLagrPar_Mod

      INTEGER, INTENT (IN)  :: jflg_vel
      REAL (KIND=IKIND), INTENT (IN)          :: aTime
      
                IF (jflg_vel==1) THEN
                            !CALL Release_LP_Vcell ()                   ! For MODFLOW-needs to be updated
                ELSE 
                             CALL Release_LP_Vnode (aTime)    ! For IWFM    
                ENDIF

  END  SUBROUTINE Release_LP
  !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
  

!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
 SUBROUTINE Advect_LP (jflg_vel )
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
     USE AdvDispLagrPar_Mod

       INTEGER, INTENT (IN)  :: jflg_vel
          !--Advect particles 
             IF ( jflg_vel==1) THEN
                   !CALL  AdvDisp_Vcell ( )                  ! For MODFLOW-needs to be updated
             ELSE
                  CALL  AdvDisp_Vnode ( )               ! For IWFM    
            ENDIF
 
    END SUBROUTINE  Advect_LP
  !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     

!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
SUBROUTINE Output (Jout, aTime ) 
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     

     USE AdvDispLagrPar_Mod
     

    INTEGER, INTENT (INOUT)  ::  Jout
    REAL (KIND=IKIND), INTENT (IN)          :: aTime
     !--Output results
        Jout=Jout+1
        IF (iFlag_Output==1) THEN
            CALL Save_Output_Nodes_OneFile( Jout, aTime)
        ELSE
            CALL Save_Output_Nodes (Jout, aTime) 
        ENDIF 
        
        WRITE(fid_dbug,*)  'Total released parcel mass=', tot_mass_rel
        WRITE(fid_dbug,*)  'Total mass specified=', MassRate*t_spill

  END SUBROUTINE Output
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
 
    
    
    
 
    
    
 


  
    
    
    
    
    
    
    
    
    
    
   