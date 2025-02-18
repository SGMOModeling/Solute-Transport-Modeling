!---------------------------------------------------------------------------------------------------------------------------------------!
    SUBROUTINE Release_LP_Vnode ( time)
 !--Purpose: Release Lagrangian Parcels based on the release rate at each source location
!---------------------------------------------------------------------------------------------------------------------------------------!

    USE AdvDispLagrPar_Mod
   ! USE AdvDisp_BilinearInterp_Mod
    !USE Class_LLNode
    !USE GenericLinkedList
    USE LinkedList_Mod

    IMPLICIT NONE
    
    !INTEGER ,  INTENT(IN)               :: itype
    REAL (KIND=8),  INTENT(IN)  ::  time
    
    INTEGER                            ::  np , mp, ix, iy, iz, iele, lyr , num_par_dt !(np=new parcel being added)
    REAL (KIND=IKIND)     ::  u_ran, v_ran, w_ran,  dtsmall , xcord, ycord, zcord, distz
    REAL (KIND=IKIND)     ::  U_Vel_interp, V_Vel_interp, W_Vel_interp, U_Vel_tot, V_Vel_tot, W_Vel_tot
    REAL (KIND=4)              ::  rand, randx, randy, randz  ! Random number generator requires real 4
   INTEGER                           :: iStat, ielem
   TYPE (LPNodeType), POINTER :: current


    
!--Initialize the particle grid particle numbers
 	   Model%npar(1:num_xgrids, 1:num_ygrids, 1:num_zgrids)=0
       Element(:)%nLP=0
       
       current=> LPList%Tail
       !
       !-Find position to add to the list
       IF(.NOT.ASSOCIATED(current)) THEN
           np=1
       ELSE
           np=LPList%Tail%ID+1
        ENDIF
       
   i=0
  DO 80  k=1, num_sources   
     IF (time >= Source(k)%sp_start    .and.  time <=Source(k)%sp_end ) THEN
         DO 90 ip=1,Source(k)%npar_dt   !num_par_dt
             
             xcord=Source(k)%Xorg
             ycord=Source(k)%Yorg
             zcord=Source(k)%Zorg
             ielem=Source(k)%ElemID

            !CALL List_Insert(1, np)  !  REMOVE np=np+1 at the bottom 
            ! CALL LPList(np)%AddNode(np,iStat)    !  Add np=np+1 at the bottom 
           CALL  LPList%InsertLP (np, xcord, ycord, zcord, ielem , iStat )  !  Add np=np+1 at the bottom 
           
           Current=>LPList%Tail


             ! Check status
                IF (iStat /= 0) THEN
                    WRITE(*,*)
                    WRITE(*,*) "Error inserting node:", np
                    WRITE(*,*)    
                    WRITE(fid_dbug,*)
                    WRITE(fid_dbug,*) "Error inserting node:", np
                    WRITE(fid_dbug,*)    
                    PAUSE
                END IF

          

            !--For continous release figure out the small time step for each release
            dtsmall=del_t/float(Source(k)%npar_dt)*ip
    
                   CALL   Gauss (1.0,0.0,randx)
                   CALL   Gauss (1.0,0.0,randy)   
                   CALL   Gauss (1.0,0.0,randz)              
         
                        
                    !--Interpolate the velocity based on the nodal velocities 
                    !--Release Element is specified in the input
                   
                    ! iele=LP(np)%Elem
                     iele=Current%Elem

                     !--Calculate the layer number
                     lyr=Get_lyr(iele, zcord)           

                     
                    CALL  Element(iele)%Vel_Interp(xcord,ycord,zcord, lyr, U_Vel_interp, V_Vel_interp,W_Vel_interp)
                     
                    CALL Get_DispCoeff ( U_Vel_interp, V_Vel_interp, W_Vel_interp) 
              
                   u_ran=(SQRT(2.0*DispCoeff_xx/dtsmall))*randx      ! Random velocity component along X-dir
                   v_ran=(SQRT(2.0*DispCoeff_yy/dtsmall))*randy      ! Random velocity component along Y-dir
                   w_ran=(SQRT(2.0*DispCoeff_zz/dtsmall))*randz      ! Random velocity component along Y-dir      
                        
                   
				   U_Vel_tot=U_Vel_interp+u_ran                                         ! Total velocity component along X-dir
				   V_Vel_tot=V_Vel_interp+v_ran                                         ! Total velocity component along Y-dir
 				   W_Vel_tot=W_Vel_interp+W_ran                                     ! Total velocity component along Z-dir                  
                                    
				   distx=U_Vel_tot*dtsmall                                                    ! Displacement along X-dir
				   disty=V_Vel_tot*dtsmall                                                    ! Displacement along Y-dir
				   distz=W_Vel_tot*dtsmall                                                    ! Displacement along Z-dir                   
                   
                   !LP(np)%xcord=LP(np)%xcord+distx                   ! New parcle location along X-dir
                   !LP(np)%ycord=LP(np)%ycord+disty                  ! New parcle location along X-dir
                   !LP(np)%zcord=LP(np)%zcord+distz                 ! New parcle location along X-dir              
                   
                   Current%xcord=Current%xcord+distx                   ! New parcle location along X-dir
                   Current%ycord=Current%ycord+disty                  ! New parcle location along X-dir
                   Current%zcord=Current%zcord+distz                 ! New parcle location along X-dir      
                   
                                        
                     !--If the particle has moved to a new element update the new element number 
                    CALL Get_Elem_LP (Current%Elem, Current%xcord, Current%ycord)
                   
                   !--Update the element number in case if the parcel has moved out of the element
                   iele=Current%Elem
                   
                  !--Account for parcels leaving the domain 
                  !--Check if the new parcel has landed on a boundary and remove it and insert at the starting point
                  !-- Check if the LP has moved out of the domain
                   CALL Check_bounds ( Current%Elem, Current%xcord, Current%ycord,  Current%zcord) 

                 IF (Current%Ingrid)     THEN       
                       iy=CEILING((YUprRight-Current%ycord )/dy)                                   
                       ix=CEILING((Current%xcord-XLwrLeft)/dx)
                       
                        IF (idz_ReadFile) THEN
                                iz=Get_lyr(iele, Current%zcord)         
                       ELSE
                                iz=CEILING((Current%zcord-ZLwr)/dz)
                       ENDIF

                        Current%Igridx=FLOAT (ix)
                        Current%Igridy=FLOAT (iy) 
                        Current%Igridz=FLOAT (iz)                     
                        Model%npar(ix,iy,iz)=Model%npar(ix,iy,iz)+1
                       iele=Current%Elem
                       !Model%ElemID(ix,iy)=iele
                       Element(Current%Elem)%nLP=Element(Current%Elem)%nLP+1
                 ENDIF

            num_par_rel=num_par_rel+1
            np=np+1                                       !  REMOVE FOR THE F77 method
   90   CONTINUE 
END IF     
        WRITE (fid_dbug,*) 'Number of released particles from location: ',  k , '=', Source(k)%npar_dt 
80 CONTINUE 
   WRITE (fid_dbug,*) 'Total number of released particles: =', num_par_rel
     
     !--TO CHECK THE INITIAL CONCENTRATION TURN THIS LOOP ON and the parcel counter  Model%npar(ix,iy)=Model%npar(ix,iy)+1 above
     !WRITE (*,*) "Number of Parcels Released :", num_par_rel
          !DO  ix=1, num_xgrids	
          !              DO  iy=1, num_ygrids
          !                  Model%Con(ix,iy)= Model%npar(ix,iy)*mass_of_parcel/(dx*dy*dz)/Porosity*1000 ! g/m3=mg/l=ppm
          !              ENDDO
          !ENDDO
   
   
  !  CALL LPList%print_list()

    
    END SUBROUTINE Release_LP_Vnode
 !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   

 !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   

SUBROUTINE      Get_DispCoeff ( U_vel, V_vel, W_vel )  
!--Purpose: Calculate Dispersion Coefficients 
 !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   

USE AdvDispLagrPar_Mod

    REAL (KIND=IKIND), INTENT(IN)     :: U_vel, V_vel, W_vel
    REAL (KIND=IKIND)                               :: Usq,Vsq,Wsq


        Usq=U_vel*U_vel
        Vsq=V_vel*V_vel
        Wsq=W_vel*W_vel
        Vtot=SQRT(Usq+Vsq+Wsq)
            
        HorTrnsDisp=HorTrns2LongDispRat*LongDisp
        VerTrnsDisp=VerTrns2LongDispRat*LongDisp
        
        DispCoeff_xx=LongDisp*(Usq)/Vtot+HorTrnsDisp*(Vsq )/Vtot  + VerTrnsDisp*(Wsq)/Vtot ! (m2/day)
        DispCoeff_yy=LongDisp*(Vsq)/Vtot+HorTrnsDisp*(Usq)/Vtot+VerTrnsDisp*(Wsq)/Vtot   ! (m2/day)
        DispCoeff_zz=LongDisp*(Wsq)/Vtot+VerTrnsDisp*(Usq)/Vtot+VerTrnsDisp*(Vsq)/Vtot   ! (m2/day)
        
        DispCoeff_xy=(LongDisp-HorTrnsDisp)*U_vel*V_vel/Vtot
        DispCoeff_xz=(LongDisp-VerTrnsDisp)*U_vel*W_vel/Vtot
        DispCoeff_xz=(LongDisp-VerTrnsDisp)*V_vel*W_vel/Vtot
                    
                    
    END SUBROUTINE Get_DispCoeff
 !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   
    
    
 !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   
    
    SUBROUTINE AdvDisp_Vnode (  )
     USE AdvDispLagrPar_Mod
     USE LinkedList_Mod
     !USE AdvDisp_BilinearInterp_Mod
 !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   
   
     !REAL(8), INTENT (IN)   :: Time
     
    INTEGER                  :: np , mp, ix, iy, iz, irand, iele, lyr, npar
    REAL (KIND=IKIND)     ::  u_ran, v_ran, w_ran, dtsmall , xcrd, ycrd, zcrd, U_Vel_interp, V_Vel_interp, W_Vel_interp, U_vel, V_vel, W_vel, U_vel_tot, V_vel_tot,W_vel_tot, elem_thick
    REAL (KIND=4)     ::  rand, randx, randy, rand_num  ! Random number generator requires real 4
    
    TYPE (LPNodeType), POINTER :: current
    LOGICAL :: found
    
    found=.FALSE.

        !np=beg(1)
        current=> LPList%Head
        np= current%ID
        num_par=0    
       !-- Reset the particle grid particle numbers
 	   Model%npar(1:num_xgrids, 1:num_ygrids, 1:num_zgrids)=0
       Element(:)%nLP=0


        !DO 10, WHILE (np  /=   -1 ) 
        DO 10, WHILE (ASSOCIATED (current))    
            
                    np=current%ID     

                    CALL   Gauss (1.0,0.0,randx)
                    CALL   Gauss (1.0,0.0,randy)                    
                    CALL   Gauss (1.0,0.0,randz)
                   !
 
                    iele=Current%Elem
                       
                     !-- Check if the LP has moved out of the domain
                     ! CALL Check_bounds ( np, iele) 
                    CALL Check_bounds (Current%Elem, Current%xcord, Current%ycord,  Current%zcord) 


                     xcrd=Current%xcord
                     ycrd=Current%ycord
                     zcrd=Current%zcord
          
                     !--Calculate the layer number
                    lyr=Get_lyr(iele, Current%zcord)       
                     

                     
                    CALL  Element(Current%Elem)%Vel_Interp(Current%xcord, Current%ycord,Current%zcord, lyr, U_Vel_interp, V_Vel_interp,W_Vel_interp)
                     
                     CALL Get_DispCoeff ( U_Vel_interp, V_Vel_interp, W_Vel_interp) 
             
                   u_ran=(SQRT(2.0*DispCoeff_xx/del_t))*randx      ! Random velocity component along X-dir
                   v_ran=(SQRT(2.0*DispCoeff_yy/del_t))*randy      ! Random velocity component along Y-dir
                   w_ran=(SQRT(2.0*DispCoeff_zz/del_t))*randz      ! Random velocity component along Z-dir      
                                       
      
				   U_Vel_tot=U_Vel_interp+u_ran                                         ! Total velocity component along X-dir
				   V_Vel_tot=V_Vel_interp+v_ran                                         ! Total velocity component along Y-dir
				   W_Vel_tot=W_Vel_interp+W_ran                                         ! Total velocity component along Z-dir                   
                                    
				   distx=U_Vel_tot*del_t                                                    ! Displacement along X-dir
				   disty=V_Vel_tot*del_t                                                    ! Displacement along Y-dir
				   distz=W_Vel_tot*del_t                                                    ! Displacement along Z-dir                   
                   
                   Current%xcord=Current%xcord+distx                   ! New parcle location along X-dir
                   Current%ycord=Current%ycord+disty                  ! New parcle location along Y-dir
                   Current%zcord=Current%zcord+distz                  ! New parcle location along Y-dir                   
                   
         

                   !--If the particle has moved to a new element figure out the new element number 
                   !CALL Get_Elem_LP (np)
                   CALL Get_Elem_LP (Current%Elem, Current%xcord, Current%ycord)

                   
                   !--Update the element number in case if the parcel has moved out of the element
                  ! iele=LP(np)%Elem
                   

                    !-- Check if the LP has moved out of the domain
                  !   CALL Check_bounds ( np, iele) 
                   CALL Check_bounds (Current%Elem, Current%xcord, Current%ycord,  Current%zcord) 

 
       
                   !--Account for sinks-Delete parcels at the pumping wells 

                    DO I =1,  npump_wells          
                       iy=CEILING((YUprRight- Current%ycord )/dy)
                        ix=CEILING(( Current%xcord-XLwrleft)/dx)
                        IF (  (ix==PumpWell(i)%IR  .and.  iy==PumpWell(i)%IC)   )THEN
                             !LP(np)%Ingrid= .FALSE.
                             !mp=Link(np)
                            !CALL List_Delete (1, np)
                             !np=mp
                            !-Delete parcels
                             CALL  LPList%RemoveLP (np, found, iStat) 
                                          ! Check status
                              IF (iStat /= 0) THEN
                                   WRITE(*,*)
                                   WRITE(*,*) "Error removing node:", np
                                   WRITE(*,*)    
                                   WRITE(fid_dbug,*)
                                   WRITE(fid_dbug,*) "Error removing node:", np
                                   WRITE(fid_dbug,*)    
                                  PAUSE
                               END IF
                             IF (found) PumpWell(i)%nParDel=PumpWell(i)%nParDel+1
                            !goto 10
                               !-Reset release location
                             !CALL   RANDOM_NUMBER(rand_num)      ! Random number (rea1) between 0-1
                             !irand=CEILING(rand_num*(num_sources))             ! Random number (integer) between 1 and number 
                             !LP(np)%xcord=Source(irand)%Xorg
                             !LP(np)%ycord=Source(irand)%Yorg  
                             !LP(np)%zcord=Source(irand)%Zorg                               
                       ENDIF
                    ENDDO
                   
                   
                   !--CHECK IF THE PARCEL HAS LEFT THE BOUNDARY
             
                 !--Deal with particle in inactive domain
                 !--Remove the particle from inactive domain and place it at the origin
                 !IF (Model%Ibnd(ix,iy) /= 1 ) THEN
                 !        !-Reset release location
                 !        CALL   RANDOM_NUMBER(rand_num)      ! Random number (rea1) between 0-1
                 !        irand=CEILING(rand_num*(num_sources))             ! Random number (integer) between 1 and number 
                 !        LP(np)%xcord=Source(irand)%Xorg
                 !        LP(np)%ycord=Source(irand)%Yorg                  
                 !ENDIF
                 
                 
                 !--Calculate the number of parcels within each model cell
                 IF (Current%Ingrid)     THEN         
                       iy=CEILING((YUprRight-Current%ycord )/dy)
                       ix=CEILING((Current%xcord-XLwrLeft)/dx)
                       IF (idz_ReadFile) THEN
                                iz=Get_lyr(Current%Elem, Current%zcord)   
                       ELSE
                                 iz=CEILING((Current%zcord-ZLwr)/dz)
                       ENDIF
                     Current%Igridx=FLOAT (ix)
                     Current%Igridy=FLOAT (iy) 
                     Current%Igridz=FLOAT (iz)                      
                     Model%npar(ix,iy,iz)=Model%npar(ix,iy,iz)+1
                     iele=Current%Elem
                     !Model%ElemID(ix,iy)=iele
                     Element(iele)%nLP=Element(iele)%nLP+1
                 ELSE 
                     WRITE (*,*)  'LP is outside the grid'
                 ENDIF

                 

    
               ! np=Link(np)
                current=>current%Next
                num_par=num_par+1
10      CONTINUE
        
        WRITE (fid_dbug,*) 'Number of particles advected=', num_par
 
        !--Make sure the number of particles in the system are same as the numbers of particles released
        !npar=0
        !        DO iz=1, num_zgrids
			     !     DO  ix=1, num_xgrids	
        !                DO  iy=1, num_ygrids
        !                    npar=npar+ Model%npar(ix,iy,iz)
        !  
        !                    if (model%npar(ix,iy,iz)>0) then
        !                           iele= Model%ElemID(ix,iy)
        !                           elem_thick=(Element(iele)%AvgLyrTopElv(iz)-Element(iele)%AvgLyrTopElv(iz+1))   
        !
        !                           if (Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity>1.0 ) then
        !                                WRITE(*,  '( 3I4, F12.3)' )   ix, iy, iz, Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity*1000  ! g/m3=ppm
        !                            endif
        !                       !write(*,*) ix, iy, iz, model%npar(ix,iy,iz)
        !                    endif
        !                ENDDO
        !              ENDDO
        !        ENDDO
 
    
    
    END SUBROUTINE AdvDisp_Vnode  
  !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   
   
    
!  !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   
!   
!    SUBROUTINE Release_LP_Vcell( )
!    
!    USE AdvDispLagrPar_Mod
!    !USE AdvDisp_BilinearInterp_Mod
!    !USE Class_LLNode
!   !USE GenericLinkedList
!  !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   
!   
!   ! INTEGER ,  INTENT(IN)               :: itype
!    !REAL (KIND=8),  INTENT(IN)  ::  dt
!    
!    INTEGER                  :: np , mp, ix, iy,  iele,num_par_dt , npar_rel_prev, lyr !(np=new parcel being added)
!    REAL (KIND=IKIND)     ::  u_ran, v_ran, w_ran, dtsmall , xcrd, ycrd, zcrd, U_Vel_interp, V_Vel_interp, W_Vel_interp, U_vel, V_vel, W_vel, U_vel_tot, V_vel_tot,W_vel_tot
!    REAL (KIND=4)     ::  rand, randx, randy,randz  ! Random number generator requires real 4
!
!    
!!--Initialize the particle grid particle numbers
! 	   Model%npar(1:num_xgrids, 1:num_ygrids, 1:num_zgrids)=0
!       Element(:)%nLP=0
!
! 
!
!npar_rel_prev=0
!  DO 80  k=1, num_sources   
!     DO 70 ip=1,Source(k)%npar_dt     !,num_par_dt
!
!            CALL List_Insert(1, np)
!                 LP(np)%Igridx=0
!                 LP(np)%Igridy=0        
!                 LP(np)%xcord=Source(k)%Xorg
!                 LP(np)%ycord=Source(k)%Yorg
!                 LP(np)%zcord=Source(k)%Zorg
!                 LP(np)%Elem=Source(k)%ElemID
!                 LP(np)%Ingrid= .TRUE.         ! within domain Ingrid=T; Outside the domain Ingrid=F
!            
!            
!            !--For continous release figure out the small time step for each release
!            dtsmall=del_t/float(Source(k)%npar_dt)*ip
!    
!                    CALL   Gauss (1.0,0.0,randx)
!                    CALL   Gauss (1.0,0.0,randy)                    
!                    CALL   Gauss (1.0,0.0,randz)
!                   !  
!
!                   u_ran=(SQRT(2.0*DispCoeff_xx/del_t))*randx      ! Random velocity component along X-dir
!                   v_ran=(SQRT(2.0*DispCoeff_yy/del_t))*randy      ! Random velocity component along Y-dir
!                   w_ran=(SQRT(2.0*DispCoeff_zz/del_t))*randz      ! Random velocity component along Z-dir      
!                                       
!                   
!                   
!                       !!---USE CELL CENTERED VELOCITY-MODFLOW
!                       !!--Figure out the cell 
!                       iy=CEILING((YUprRight-LP(np)%ycord )/dy)  ! Row
!                       ix=CEILING((LP(np)%xcord-XLwrleft)/dx)    ! Col
!                       ! iz= CEILING((LP(np)%zcord-ZLwr)/dz)
!                      iz=Get_lyr(iele, zcrd)           
!
!                       
!                        !--Figure out the velocity
!                        U_vel=Model%Uvel(ix,iy)
!                        V_vel=Model%Vvel(ix,iy) 
!                        W_vel=Model%Wvel(ix,iy)                     
!                   
!                     iele=LP(np)%Elem
!                     xcrd=LP(np)%xcord
!                     ycrd=LP(np)%ycord
!                     zcrd=LP(np)%zcord
!                    
!                     
!                     !--Interpolate the velocity based on the MODFLOW cell face velocities -Need to be implemented
!         !            CALL  Element(iele)%Vel_Interp_Quad(xcrd, ycrd, U_Vel_interp, V_Vel_interp)
!				     !U_Vel_tot=U_Vel_interp+u_ran                                         ! Total velocity component along X-dir
!				     !V_Vel_tot=V_Vel_interp+v_ran                                         ! Total velocity component along Y-dir
!                     
!                    U_Vel_tot=U_vel+u_ran                                         ! Total velocity component along X-dir
!				    V_Vel_tot=V_vel+v_ran                                         ! Total velocity component along Y-dir
!                    W_vel_tot=W_vel+w_ran
!                                    
!				   distx=U_Vel_tot*dtsmall                                                    ! Displacement along X-dir
!				   disty=V_Vel_tot*dtsmall                                                    ! Displacement along Y-dir
!                   distz=W_Vel_tot*dtsmall       
!                   
!                   LP(np)%xcord=LP(np)%xcord+distx                   ! New parcle location along X-dir
!                   LP(np)%ycord=LP(np)%ycord+disty                  ! New parcle location along Y-dir
!                   LP(np)%zcord=LP(np)%zcord+distz                  ! New parcle location along Z-dir
!                   
!                   LPx=LP(np)%xcord
!                   LPy=LP(np)%ycord   
!                   LPy=LP(np)%zcord     
!                   
!                     !--If the particle has moved to a new element update the new element number 
!                  ! CALL Get_Elem_LP (np)
!                   
!                 !--CHECK IF THE PARCEL HAS LEFT THE BOUNDARY
!                 !  CALL Check_bounds ( np, iele)   -Needs to be tested for MODFLOW
!
!                 IF (LP(np)%Ingrid)     THEN       
!                       iy=CEILING((YUprRight-LP(np)%ycord )/dy)                                   
!                       ix=CEILING((LP(np)%xcord-XLwrLeft)/dx)
!                     ! iz=CEILING((LP(np)%zcord-ZLwr)/dz)
!
!                       iz=Get_lyr(iele, zcrd) 
!                       
!                     LP(np)%Igridx=FLOAT (ix)
!                    LP(np)%Igridy=FLOAT (iy)
!                    LP(np)%Igridz=FLOAT (iz)                    
!                    Model%npar(ix,iy, iz)=Model%npar(ix,iy, iz)+1
!                     iele=LP(np)%Elem
!                     Element(iele)%nLP=Element(iele)%nLP+1
!                 ENDIF
!
!            num_par_rel=num_par_rel+1
!70   CONTINUE 
!     
!        !WRITE (*,*) 'Number of released particles=', num_par_rel
!     
!        WRITE (fid_dbug,*) 'Number of released particles from location: ',  k , '=', Source(k)%npar_dt 
!80 CONTINUE 
!   WRITE (fid_dbug,*) 'Total number of released particles: =', num_par_rel
!     
!     !--TO CHECK THE INITIAL CONCENTRATION TURN THIS LOOP ON and the parcel counter  Model%npar(ix,iy)=Model%npar(ix,iy)+1 above
!     !WRITE (*,*) "Number of Parcels Released :", num_par_rel
!          !DO  ix=1, num_xgrids	
!          !              DO  iy=1, num_ygrids
!          !                  Model%Con(ix,iy)= Model%npar(ix,iy)*mass_of_parcel/(dx*dy*dz)/Porosity*1000 ! g/m3=mg/l=ppm
!          !              ENDDO
!          !ENDDO
!  
!    
!    END SUBROUTINE Release_LP_Vcell
! !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   
!    
    
    
!  !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   
!   
!    SUBROUTINE AdvDisp_Vcell (  )
!     USE AdvDispLagrPar_Mod
!     !USE AdvDisp_BilinearInterp_Mod
!   
!     !REAL(8), INTENT (IN)   :: Time
!     
!    INTEGER                  :: np , mp, ix, iy, irand, iele, lyr
!    REAL (KIND=IKIND)     ::u_ran, v_ran, w_ran, dtsmall , xcrd, ycrd, zcrd, U_Vel_interp, V_Vel_interp,  W_Vel_interp, U_vel, V_vel, W_vel, U_vel_tot, V_vel_tot,W_vel_tot
!    REAL (KIND=4)     ::  rand, randx, randy, randz, rand_num  ! Random number generator requires real 4
!
!        np=beg(1)
!        num_par=0    
!       !-- Reset the particle grid particle numbers
! 	   Model%npar(1:num_xgrids, 1:num_ygrids, 1:num_zgrids)=0
!       
!       Element(:)%nLP=0
!
!        DO 10, WHILE (np  /=   -1 ) 
!
!                    CALL   Gauss (1.0,0.0,randx)
!                    CALL   Gauss (1.0,0.0,randy)                    
!                    CALL   Gauss (1.0,0.0,randz)
!
!                   u_ran=(SQRT(2.0*DispCoeff_xx/del_t))*randx      ! Random velocity component along X-dir
!                   v_ran=(SQRT(2.0*DispCoeff_yy/del_t))*randy      ! Random velocity component along Y-dir
!                   w_ran=(SQRT(2.0*DispCoeff_zz/del_t))*randz      ! Random velocity component along Z-dir      
!                                       
!         
!                      !--Use Cell velocity (MODLFLOW)
!                       !!!--Figure out the cell 
!                       iy=CEILING((YUprRight-LP(np)%ycord )/dy)        !(29,15) element 1303
!                       ix=CEILING((LP(np)%xcord-XLwrleft)/dx)
!                       ! iz=CEILING((LP(np)%zcord-ZLwr)/dz)
!                       iz=Get_lyr(iele, zcrd) 
!                                   
!                        !--Figure out the velocity
!                        U_vel=Model%Uvel(ix,iy)
!                        V_vel=Model%Vvel(ix,iy) 
!                        W_vel=Model%Wvel(ix,iy)         
!                        
!                     !--Interpolate the velocity based on the MODFLOW cell face velocities -Need to be implemented
!                    !  iele=LP(np)%Elem
!                    ! xcrd=LP(np)%xcord
!                    ! ycrd=LP(np)%ycord
!                    !CALL  Element(iele)%Vel_Interp_Quad(xcrd, ycrd, U_Vel_interp, V_Vel_interp)
!                     !U_Vel_tot=U_Vel_interp+u_ran                                         ! Total velocity component along X-dir
!				    !V_Vel_tot=V_Vel_interp+v_ran                                         ! Total velocity component along Y-dir
!
!				   U_Vel_tot=U_vel+u_ran                                         ! Total velocity component along X-dir
!				   V_Vel_tot=V_vel+v_ran                                         ! Total velocity component along Y-dir
!				   W_Vel_tot=W_vel+v_ran                                         ! Total velocity component along Y-dir                   
!                                    
!				   distx=U_Vel_tot*del_t                                                    ! Displacement along X-dir
!				   disty=V_Vel_tot*del_t                                                    ! Displacement along Y-dir
!				   distz=W_Vel_tot*del_t                                                    ! Displacement along Y-dir                   
!                   
!                   LP(np)%xcord=LP(np)%xcord+distx                   ! New parcle location along X-dir
!                   LP(np)%ycord=LP(np)%ycord+disty                  ! New parcle location along Y-dir
!                   LP(np)%zcord=LP(np)%zcord+distz                  ! New parcle location along Z-dir                   
!                   
!                   LPx=LP(np)%xcord
!                   LPy=LP(np)%ycord         
!                   LPz=LP(np)%zcord                            
!
!                   !--If the particle has moved to a new element figure out the new element number 
!                   CALL Get_Elem_LP (np)
!                   
!       
!                   !--Account for sinks-Delete parcels at the pumping wells 
!                    iy=CEILING((YUprRight- LPy )/dy)
!                    ix=CEILING(( LPx-XLwrleft)/dx)
!                    !iz=CEILING((LPz-ZLwr)/dz)
!                    iz=Get_lyr(iele, zcrd) 
!                    DO I =1,  npump_wells                     
!                        IF (  (ix==PumpWell(i)%IR  .and.  iy==PumpWell(i)%IC)   )THEN
!                             !LP(np)%Ingrid= .FALSE.
!                             !mp=Link(np)
!                            !CALL List_Delete (1, np)
!                             !np=mp
!                            PumpWell(i)%nParDel=PumpWell(i)%nParDel+1
!                            !goto 10
!                               !-Reset release location
!                             CALL   RANDOM_NUMBER(rand_num)      ! Random number (rea1) between 0-1
!                             irand=CEILING(rand_num*(num_sources))             ! Random number (integer) between 1 and number 
!                             LP(np)%xcord=Source(irand)%Xorg
!                             LP(np)%ycord=Source(irand)%Yorg 
!                             LP(np)%zcord=Source(irand)%Zorg      
!                       ENDIF
!                   ENDDO
!                   
!                    !--CHECK IF THE PARCEL HAS LEFT THE BOUNDARY
!                 !  CALL Check_bounds ( np, iele) -Needs to be tested for MODFLOW
!
!                 !--Deal with particle in inactive domain
!                 !--Remove the particle from inactive domain and place it at the origin
!                 !IF (Model%Ibnd(ix,iy) /= 1 ) THEN
!                 !        !-Reset release location
!                 !        CALL   RANDOM_NUMBER(rand_num)      ! Random number (rea1) between 0-1
!                 !        irand=CEILING(rand_num*(num_sources))             ! Random number (integer) between 1 and number 
!                 !        LP(np)%xcord=Source(irand)%Xorg
!                 !        LP(np)%ycord=Source(irand)%Yorg                  
!                 !ENDIF
!                 
!                 
!                 !--Calculate the number of parcels within each model cell
!                 IF (LP(np)%Ingrid)     THEN         
!                       iy=CEILING((YUprRight-LP(np)%ycord )/dy)
!                       ix=CEILING((LP(np)%xcord-XLwrLeft)/dx)
!                       !iz=CEILING((LP(np)%zcord-ZLwr)/dz)
!                       iz=Get_lyr(iele, zcrd) 
!                     LP(np)%Igridx=FLOAT (ix)
!                     LP(np)%Igridy=FLOAT (iy) 
!                     LP(np)%Igridz=FLOAT (iz) 
!                     Model%npar(ix,iy, iz)=Model%npar(ix,iy,iz)+1
!                     iele=LP(np)%Elem
!                     Element(iele)%nLP=Element(iele)%nLP+1
!                 ENDIF
!
!                 
!
!    
!                np=Link(np)
!                num_par=num_par+1
!10      CONTINUE
!        
!        WRITE (fid_dbug,*) 'Total umber of particles advected=', num_par
!    
!    
!    END SUBROUTINE AdvDisp_Vcell 
! !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   
!    
!    
    
    
    

    
    
    
    