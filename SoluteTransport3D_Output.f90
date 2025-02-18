!***********************Program LPRWSTM****************************************!
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
! **********************************************************************************!
      
    
    SUBROUTINE Save_Output_Nodes_OneFile (Iout, Time )
        USE AdvDispLagrPar_Mod
        
        INTEGER, INTENT (IN)  :: Iout
        REAL(KIND=IKIND), INTENT (IN)   :: Time
        
        
        INTEGER                             :: kx, ky, m,  npar, npar_del_tot, nl, ielem, np, fid_ouput_Node,  fid_out_YZplane, fid_out_XYplane, fid_out_XZplane, ngrids
        REAL(KIND=IKIND)       :: CN, U_Vel, V_Vel, elem_thick
        REAL (KIND=IKIND)      :: x, y,z, xcrd, ycrd, zcord, Dxx, Dyy, Dzz, term1, term2, RR, UU
       CHARACTER (LEN=10)     :: IO       
        
        
            !--Output the calculated concentrations at each node
            avg_conc_sim=0.00
            avg_conc_ana=0.00
            
            !WRITE (fid_ouput_Node_All,  '(A18, A20, A19, 3A7,A10, 2A16, A15,  A10, A15 )')   'XCOORD', 'YCOORD', 'ZCOORD', 'ROW','COL','LYR','ELEMENT', 'LYR_TOP_ELEV', 'LYR_BOT_ELEV' , 'ELEM_THICK', 'NUM_LP',  'CONC_PPM'
            !--All in one file
            ngrids=num_xgrids* num_ygrids*num_zgrids
            WRITE( IO, '(F10.0)') Time
            WRITE (fid_ouput_Node_All, '(A, A,6x)', ADVANCE='NO') "Conc_", TRIM(ADJUSTL(IO))
            DO iz=1, num_zgrids
                DO  ix=1, num_xgrids	
                    DO  iy=1, num_ygrids
                            ielem= Model%ElemID(ix,iy)
                            IF (idz_ReadFile) THEN
                                    elem_thick=(Element(ielem)%AvgLyrTopElv(iz)-Element(ielem)%AvgLyrTopElv(iz+1))  
                            ELSE
                                     elem_thick=dz
                            ENDIF
                            WRITE(fid_ouput_Node_All,  '(<ngrids>(F12.3))'  , ADVANCE='NO'  )   Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity*1000  ! g/m3=ppm
                    ENDDO
                 ENDDO
            ENDDO
            WRITE(fid_ouput_Node_All,*)
            
          !  
          !!--XY plane
          !  DO  ix=1, num_xgrids	
          !      DO  iy=1, num_ygrids
          !              ielem= Model%ElemID(ix,iy)
          !              elem_thick=(Element(ielem)%AvgLyrTopElv(iz)-Element(ielem)%AvgLyrTopElv(iz+1))   
          !              WRITE(fid_ouput_Node_All,  '(<ngrids>(F12.4,2x))'  , ADVANCE='NO'  )   Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity*1000  ! g/m3=ppm
          !      ENDDO
          !  ENDDO   
          !  
          !!--XZ plane
          !      DO  iz=1, num_zgrids	
          !          DO  ix=1, num_xgrids
          !              ielem= Model%ElemID(ix,iy)
          !              elem_thick=(Element(ielem)%AvgLyrTopElv(iz)-Element(ielem)%AvgLyrTopElv(iz+1))   
          !              WRITE(fid_ouput_Node_All,  '(<ngrids>(F12.4,2x))'  , ADVANCE='NO'  )   Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity*1000  ! g/m3=ppm
          !      ENDDO
          !  ENDDO   
   
    
  END   SUBROUTINE Save_Output_Nodes_OneFile
    
    
    
    
    
    
    
    
    
    
    
    
    
    SUBROUTINE Save_Output_Nodes ( Iout , Time)
    
        USE AdvDispLagrPar_Mod
        
        INTEGER, INTENT (IN)  :: Iout
        REAL(KIND=IKIND), INTENT (IN)   :: Time
        
        
        INTEGER                             :: kx, ky, m,  npar, npar_del_tot, nl, ielem, np, fid_ouput_Node,  fid_out_YZplane, fid_out_XYplane, fid_out_XZplane
        REAL(KIND=IKIND)       :: CN, U_Vel, V_Vel, elem_thick
        REAL (KIND=IKIND)      :: x, y,z, xcrd, ycrd, zcord, Dxx, Dyy, Dzz, term1, term2, RR, UU
       CHARACTER (LEN=4)     :: IO       
        
        


            !--Output the calculated concentrations at each node
            fid_ouput_Node=fid_ouput_Node_start+IOut*10
            WRITE( IO, '(I3)') IOut
            avg_conc_sim=0.00
            avg_conc_ana=0.00

  	        OPEN (fid_ouput_Node, FILE=  'Output\ConLagrng_3D_'//trim(adjustl(IO))//trim(adjustl('.dat ')) , STATUS='UNKNOWN')
            WRITE (fid_ouput_Node,  '(A18, A20, A19, 3A7,A10, 2A16, A15,  A10, A15 )')   'XCOORD', 'YCOORD', 'ZCOORD', 'ROW','COL','LYR','ELEMENT', 'LYR_TOP_ELEV', 'LYR_BOT_ELEV' , 'ELEM_THICK', 'NUM_LP',  'CONC_PPM'
            npar=0
                DO iz=1, num_zgrids
			          DO  ix=1, num_xgrids	
                        DO  iy=1, num_ygrids
                            xcrd=ix*dx+XLwrLeft
                            ycrd=YUprRight-iy*dy
                              IF (idz_ReadFile) THEN
                                  ielem= Model%ElemID(ix,iy)
                                  elem_thick=(Element(ielem)%AvgLyrTopElv(iz)-Element(ielem)%AvgLyrTopElv(iz+1))   
                                  zcrd=0.5*(Element(ielem)%AvgLyrTopElv(iz)+Element(ielem)%AvgLyrTopElv(iz+1))
                                  !elem_thick=(Model%AvgLyrTElv(iz)-Model%AvgLyrTElv(iz+1))
                                  !zcrd=0.5*(Model%AvgLyrTElv(iz)+Model%AvgLyrTElv(iz+1))
                              ELSE
                                   zcrd=dz/2+(iz-1)*dz+ZLwr
                             ENDIF
                                WRITE (fid_ouput_Node, '( 2(f20.2,1x),f15.2,1x, 3I6, I10, 3(f15.2,1x), I10,  (f15.3,1x)  )')  xcrd, ycrd, zcrd, ix, iy, iz, Element(ielem)%ElemID, Element(ielem)%AvgLyrTopElv(iz), Element(ielem)%AvgLyrTopElv(iz+1),   elem_thick, Model%npar(ix,iy,iz), Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity*1000  ! g/m3=ppm
                            npar=npar+Model%npar(ix, iy ,iz)
                        ENDDO
                      ENDDO
                ENDDO
                
            !     !--LP locations output
            !    !DO i=1, num_par_rel
            !    !    WRITE (fid_ouput_Node, '( 3(f20.2,1x),  1(f15.3,1x)  )')  LP(i)%xcord, LP(i)%ycord, LP(i)%zcord, mass_of_parcel 
            !    !ENDDO
            !    
            !
            ! !    !--Output in XY-2D plane
            !    fid_out_XYplane=fid_out_XYplane_Start+IOut*10
            !    WRITE( IO, '(I3)') IOut
            !    OPEN (fid_out_XYplane, FILE=  '..\Output\ConLagrng_XY_'//trim(adjustl(IO))//trim(adjustl('.dat ')) , STATUS='UNKNOWN')
            !
            !        ! iz=Source(1)%Zorg/dz   ! layer 15
            !        ! iz=Source(1)%Zorg/dz+1 ! layer 16
            !         !iz=Source(1)%Zorg/dz+2   ! layer 17                
            !         
            !        iz=2   ! layer 2
            !        elem_thick=Model%AvgLyrTElv(iz)   
            !        dz=elem_thick             
			         ! DO  ix=1, num_xgrids	
            !            DO  iy=1, num_ygrids
            !                xcrd=ix*dx+XLwrLeft
            !                ycrd=iy*dy+YLwrLeft
            !                zcrd=dz/2+(iz-1)*dz+ZLwr
            !                np=Model%npar(ix,iy,iz) 
            !                !IF (np>0 ) THEN
            !                    ielem= Model%ElemID(ix,iy)
            !                    elem_thick=(Element(ielem)%AvgLyrTopElv(iz)-Element(ielem)%AvgLyrTopElv(iz+1))
            !                    zcrd=0.5*(Element(ielem)%AvgLyrTopElv(iz)+Element(ielem)%AvgLyrTopElv(iz+1))
            !                    WRITE (fid_out_XYplane, '( 3(f20.2,1x),  (f15.3,1x)  )') xcrd, ycrd, zcrd, Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity*1000  ! g/m3=ppm
            !                !ELSE 
            !                !    WRITE (fid_out_XYplane, '( 3(f20.2,1x),  2(f15.4,1x)  )') xcrd, ycrd, zcrd, Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*dz)/Porosity*1000 , Model%ConAnal(ix,iy,iz) ! g/m3=ppm
            !                !ENDIF
            !            ENDDO
            !          ENDDO
            !    
            !  !    !--Output in YZ-2D plane
            !      fid_out_YZplane=fid_out_YZplane_Start+IOut*10
            !      WRITE( IO, '(I3)') IOut
   	        !     OPEN (fid_out_YZplane, FILE=  '..\Output\ConLagrng_YZ_'//trim(adjustl(IO))//trim(adjustl('.dat ')) , STATUS='UNKNOWN')                 
            !     
            !   ! ix=Source(1)%Xorg/dx+1
            !     
            !    ! ix=33
            !     ix=33
            !     xcrd=ix*dx+XLwrLeft
            !     !elem_thick=100.0
            !     !dz=elem_thick
            !     DO  iy=1, num_ygrids
			         ! DO  iz=1, num_Lyrs                           
            !                ycrd=iy*dy+YLwrLeft
            !                np=Model%npar(ix,iy,iz)
            !                ielem= Model%ElemID(ix,iy)
            !                elem_thick=(Element(ielem)%AvgLyrTopElv(iz)-Element(ielem)%AvgLyrTopElv(iz+1))
            !                zcrd=0.5*(Element(ielem)%AvgLyrTopElv(iz)+Element(ielem)%AvgLyrTopElv(iz+1))       ! Variable zcord makes the transitions from figure to figure not smooth (in the movie)
            !                !zcrd=0.5*( Model%AvgLyrTElv(iz)+Model%AvgLyrTElv(iz+1) )
            !                WRITE (fid_out_YZplane, '( 3(f20.2,1x),  (f15.3,1x)  )') xcrd, ycrd, zcrd, Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity*1000  ! g/m3=ppm
            !            ENDDO
            !     ENDDO
            !
            !!    !--Output in XZ-2D plane
            !      fid_out_XZplane=fid_out_XZplane_Start+IOut*10
            !      WRITE( IO, '(I3)') IOut
   	        !     OPEN (fid_out_XZplane, FILE=  '..\Output\ConLagrng_XZ_'//trim(adjustl(IO))//trim(adjustl('.dat ')) , STATUS='UNKNOWN')      
            !     
            !
            !    !iy=Source(1)%Xorg/dy+1        
            !     
            !     iy=15 !14
            !    ! iy=15  !15
            !     !elem_thick=100.0
            !     !dz=elem_thick
            !     
            !    DO  ix=1, num_xgrids	
			         ! DO  iz=1, num_Lyrs                 
            !                xcrd=ix*dx+XLwrLeft
            !                np=Model%npar(ix,iy,iz)
            !                
            !               !IF (np>0 ) THEN
            !                    ielem= Model%ElemID(ix,iy)
            !                    elem_thick=(Element(ielem)%AvgLyrTopElv(iz)-Element(ielem)%AvgLyrTopElv(iz+1))
            !                    zcrd=0.5*(Element(ielem)%AvgLyrTopElv(iz)+Element(ielem)%AvgLyrTopElv(iz+1))       ! Variable zcord makes the transitions from figure to figure not smooth (in the movie)
            !                   !zcrd=0.5*( Model%AvgLyrTElv(iz)+Model%AvgLyrTElv(iz+1) )
            !                    WRITE (fid_out_XZplane, '( 3(f20.2,1x),  (f15.3,1x)  )') xcrd, ycrd, zcrd, Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity*1000  !! g/m3=ppm
            !                !ELSE                                              
            !                !     dz=Model%AvgLyrTElv(iz)   
            !                !     zcrd=dz/2+(iz-1)*dz+ZLwr
            !                !    WRITE (fid_out_XZplane, '( 3(f20.2,1x),  2(f15.4,1x)  )') xcrd, ycrd, zcrd, Model%npar(ix, iy ,iz)*mass_of_parcel/(dx*dy*dz)/Porosity*1000 , Model%ConAnal(ix,iy ,iz)  ! g/m3=ppm
            !                !ENDIF
            !            ENDDO
            !     ENDDO

                 
               

              
             !--COMMENTED FOR ANALYTICAL SOLUTION   
             !      DO iz=1, num_Lyrs-1    !num_zgrids
			          !DO  ix=1, num_xgrids	
             !           DO  iy= num_ygrids, 1-1
             !               xcrd=ix*dx+XLwrLeft
             !               ycrd=YUprRight-iy*dy
             !               np=MAXVAL(Model%npar(ix,iy,:)) 
             !               elem_thick=100.0
             !               IF (np>0 ) THEN
             !                   ielem= Model%ElemID(ix,iy)
             !                   elem_thick=(Element(ielem)%AvgLyrTopElv(iz)-Element(ielem)%AvgLyrTopElv(iz+1))
             !                   zcrd=0.5*(Element(ielem)%AvgLyrTopElv(iz)+Element(ielem)%AvgLyrTopElv(iz+1))
             !                   WRITE (fid_ouput_Node, '( 3(f20.2,1x),  (f15.3,1x)  )') xcrd, ycrd, zcrd, Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity  ! kg/m2
             !
             !               ENDIF
		           !      !WRITE (fid_ouput_Node, '( 3(f20.2,1x),  (f15.3,1x)  )') xcrd, ycrd, zcrd, Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity  ! kg/m2
             !           ENDDO
             !         ENDDO
             !      ENDDO
                   
                 
                 avg_conc_sim=avg_conc_sim/float(num_nodes)
                 avg_conc_ana=avg_conc_ana/float(num_nodes)
                 
                 WRITE(fid_dbug,*) 
                 WRITE(fid_dbug,*) 'Average simulated conc.=', avg_conc_sim
                 WRITE(fid_dbug,*) 'Average analytical conc.=', avg_conc_ana
                 WRITE(fid_dbug,*) 
                 
               !--Output the analytical concentrations
  	             !OPEN (fid_ouput_Anal, FILE= '..\Output\ConAnal.dat ', STATUS='UNKNOWN')
                ! DO  iy=1,num_ygrids
		              !  WRITE (fid_ouput_Anal,'(<num_xgrids>(f15.3,1x) )') (Model%ConAnal(ix, iy), ix= 1,num_xgrids)
                ! ENDDO
                npar_del_tot=0
                DO I =1,  npump_wells
                        WRITE (*,*) "Total parcel deleted at Pumping Wel",i ,  PumpWell(i)%nParDel
                        WRITE (fid_dbug,*) "Total parcel deleted at Pumping Wel",i ,  PumpWell(i)%nParDel   
                         npar_del_tot= npar_del_tot+PumpWell(i)%nParDel
                ENDDO 
                
                
               WRITE (*,*) 'Total number of particles released=',  npar+npar_del_tot
               WRITE (*,*)
               WRITE (fid_dbug,*) 'Total number of particles released=',  npar+npar_del_tot
               WRITE (fid_dbug,*)
               
             !-- Reset grid info for new time step
 		    Model%npar(1:num_xgrids, 1:num_ygrids, 1:num_zgrids)=0


        76   FORMAT ('X = ',f7.1,' (m),',1x\)    
        77		FORMAT (/'Y = ',f7.2,' (m)',1x \)
        78		FORMAT (f15.3,1x \)
      
  
           
        WRITE (fid_dbug,*) 'Mass of a parcel=', mass_of_parcel
           
        WRITE (*,*)
        WRITE (*,*) 'Total Number of Particles Finally Available=', npar
        WRITE (fid_dbug,*)
        WRITE (fid_dbug,*) 'Total Number of Particles Finally Available=', npar
        !
        !   
           
        !   CALL Output (Iout) 
          
           
    END SUBROUTINE Save_Output_Nodes
    
    
    