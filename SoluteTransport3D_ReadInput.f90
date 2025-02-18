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
    
    !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
    SUBROUTINE  ReadInput ( )
!--Purpose:  Read input data from files    
    !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
    
    USE AdvDispLagrPar_Mod
    

    idz_ReadFile=.FALSE.
    
    OPEN (fid_dbug,FILE='Output\dbug_Lgrng.dat',STATUS='REPLACE')
    OPEN (fid_Geninf,FILE='Input\GenInfo.dat',STATUS='OLD')
    
    READ (fid_Geninf,*) t_sim			            !Simulation duration (days)
	READ (fid_Geninf,*) del_t			            !Delta t (days)
	READ (fid_Geninf,*) nt_skip		            !No: of time steps skiped for output
    READ (fid_Geninf,*)  num_Lyrs           ! Number of layers
	READ (fid_Geninf,*) dx			            !dx (m)
	READ (fid_Geninf,*) dy			            !dy (m)  
    READ (fid_Geninf,*) dz                           ! Layer Thickness
    IF (dz ==0 ) THEN
        idz_ReadFile=.TRUE.                                             ! Read elevation data from another file
    ENDIF 
    !READ (fid_Geninf,*) num_nodes
    READ (fid_Geninf,*) num_sources
    ALLOCATE (Source(num_sources))
    DO I =1, num_sources
        READ (fid_Geninf,*)   Source(i)%Xorg, Source(i)%Yorg,  Source(i)%Zorg, Source(i)%ElemID, Source(i)%Conc, Source(i)%FlowRate, Source(i)%sp_start, Source(i)%sp_end  !, Source(i)%ist_day, Source(i)%iend_day
        Source(i)%t_spill=Source(i)%sp_end-Source(i)%sp_start
    ENDDO
    READ (fid_Geninf,*) XLwrLeft, YLwrLeft, ZLwr    ! Coordinate of the upper left corner
    READ (fid_Geninf,*) XUprRight, YUprRight, ZUpr    ! Coordinate of the upper right corner
   READ (fid_Geninf,*) ReleaseType
    READ (fid_Geninf,*) num_parcels              ! Number of Lagrangian Parcels
    
    
    t_spill=MAXVAL(Source(:)%t_Spill)
    
      !-Instantaneous release or continous release    
        IF (ReleaseType==1) THEN
            INST= .TRUE.
            WRITE(*,*)
            WRITE(*,*) 'Instantaneous Release'
            WRITE(*,*)
            
            WRITE(fid_dbug,*)
            WRITE(fid_dbug,*) 'Instantaneous Release'
            WRITE(fid_dbug,*)
        ELSE
            INST= .FALSE.
            WRITE(*,*)
            WRITE(*,*) 'Continous Release with ', num_sources, 'release locations'
            WRITE(*,*)
            
            WRITE(fid_dbug,*)
            WRITE(fid_dbug,*) 'Continous Release with ', num_sources, 'release locations'
            WRITE(fid_dbug,*)
        ENDIF

    !--Convert Source concentration into mass rate (kg/day/per release location/unit depth)
    TotalMass=0.0d00
    DO I=1, num_sources
         Source(I)%MassRate=Source(i)%Conc/1000* Source(i)%FlowRate ! (kg/m3 * m3/day= kg/day)
         TotalMass=TotalMass+Source(I)%MassRate*Source(i)%t_spill
        WRITE (*,55) 'Mass release rate at location:', i, "=", Source(i)%MassRate, 'kg/day'
        WRITE (fid_dbug, 55) 'Mass release rate at location:', i, "=", Source(i)%MassRate, 'kg/day'
    ENDDO
    
    !--NEEDS TO BE UPDATED FOR MULTIPLE SOURCES WITH DIFFERENT SPILL DURATIONS
   ! t_spill=Source(1)%t_spill
    
    !MassRate=Src_Conc/1000* Q_src*num_sources/dz  ! (kg/m3 * m3/day= kg/day/m) MassRate is in mass per unit lenght (lenght taken as dz)
    
     
     WRITE (*,60)
     WRITE (*,60) 'Total Mass released =', TotalMass/1000, 'kg/1000'
     WRITE (*,60)
     
     WRITE (fid_dbug, 60) 
     WRITE (fid_dbug, 60) 'Total Mass released =', TotalMass/1000, 'kg/1000'
     WRITE (fid_dbug, 60) 
     
      !--Report the input data 
      WRITE (fid_dbug, *) 'Source Data'
      WRITE (fid_dbug, '(A15,1x, A10,7x, A10, 1x, A10, A15, 1x, A10,2x, A10, 1x, A10 )')   "ElementID", "Xorigin", "Yorigin", "Zorigin", "Concentration", "FlowRate", "StartDate", "EndDate"
        DO I =1, num_sources
            WRITE (fid_dbug, '(I10, 2x, 2(F15.2, 2x), F8.2,  2(F12.2, 2x), 2(F10.1, 2x), 2(I5, 2x))')  Source(i)%ElemID,  Source(i)%Xorg, Source(i)%Yorg,  Source(i)%Zorg, Source(i)%Conc, Source(i)%FlowRate, Source(i)%sp_start, Source(i)%sp_end  !, Source(i)%ist_day, Source(i)%iend_day
       ENDDO
        

    CLOSE (fid_Geninf)

    !--READ INPUT DATA FOR SIMULATION
    OPEN (fid_Data, FILE='Input\Data.dat',STATUS='OLD')
    

	READ (fid_Data,*) LongDisp         !Longitudinal Dispersivity (m)   !DispCoeff_xx	    !DispCoeff_xx-Dispersion Coefficient (m2/day)
    READ (fid_Data,*) HorTrns2LongDispRat  ! Ratio of  Horizontal Transverse to Longitudinal Dispersivity    
    READ (fid_Data,*) VerTrns2LongDispRat  ! Ratio of Vertical  Transverse to Longitudinal Dispersivity
    READ (fid_Data,*) Porosity           ! Effective Porosity 
    READ (fid_Data,*) FileIbnd           ! Ibound-active boundary 
    IF (idz_ReadFile) THEN
        READ (fid_Data,*) FileLayElev     ! Layer elevation data from C2VSimFG hydrostrat file
    ENDIF
    READ (fid_Data,*) FileVelInfo     ! Velocity information file
    READ (fid_Data,*)  FileElem           ! Element data file
    READ (fid_Data,*)  FileNode           ! Node data file  
    READ (fid_Data,*)  iFlag_Output
    READ (fid_Data,*)  ID_pw  ! Pumping well availalbe ID_pw=1; not available ID_pw=0
    IF (ID_pw==1) THEN
        READ (fid_Data,*) FilePumpWell
        CALL ReadPumpWell ( )
    ENDIF


    
    
    !!--OPEN VELOCITY INFORMATION FILE
       OPEN (fid_vel_info, FILE='Input\' //FileVelInfo//'' ,STATUS='OLD')
       
         READ (fid_vel_info,*) iflg_vel             !  iflg=0 => one velocity for the entire domain; iflg_vel=1 => velocity from cell/element (MODFLOW), iflg_vel=2 =>velocity from FEM node (IWFM)
         IF (iflg_vel>0) THEN
             READ (fid_vel_info,*)  iflg_vel_SS      ! 1-steady state velocity   2-transient velocity from IWFM
             READ (fid_vel_info,*) FileVel              ! Velocity file
             READ (fid_vel_info,*) VelConvFact   ! Conversion factor for veloity to convert to ft/day
         
             IF (iflg_vel_SS==1) THEN
                 SteadyStateVel=.TRUE.
             ELSE
                 SteadyStateVel=.FALSE.
             ENDIF
        ENDIF

55   FORMAT (A35, I2, 2x, A5, F12.3, A10)           
60   FORMAT (A35, F12.3, A10)
     
    CLOSE (fid_Data)
    
        
    
    END  SUBROUTINE  ReadInput 
    !---------------------------------------------------------------------------------------------------------------------------------------------------------!   

   
   
    
    
    !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
   SUBROUTINE Read_Mesh ( )
    !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
        Use AdvDispLagrPar_Mod
          INTEGER    :: mmax,  nnodes, num_nodes_all, mnode,m, ii
          INTEGER, ALLOCATABLE :: NodeList(:)
          
          
   !--Read Node input file-------------------------------------------------------------------------------------------!    
       OPEN (fid_Node, FILE='Input\' //FileNode//'' ,STATUS='OLD')
     
         READ (fid_Node, *) num_nodes
         READ (fid_Node, *)
     
         !--Allocate and Initialize Node array
         ALLOCATE (Node(num_nodes))
         ALLOCATE (NodeList(num_nodes*4))
        NodeList(:)=0
          DO i=1, num_nodes
             Node(i)%ID=0
             Node(i)%Elem(:)=0
             Node(i)%ElemCount=0
             Node(i)%xcord=0.0d00
             Node(i)%ycord=0.0d00   
             Node(i)%Conc=0.0d00  
             Node(i)%ConcAna=0.0d00 
             Node(i)%Vel=0.0d00 
            ALLOCATE( Node(i)%UVel(num_Lyrs))
            ALLOCATE( Node(i)%VVel(num_Lyrs))
            ALLOCATE( Node(i)%WVel(num_Lyrs))
             Node(i)%UVel(1:num_Lyrs)=0.0d00 
             Node(i)%VVel(1:num_Lyrs)=0.0d00 
             Node(i)%WVel(1:num_Lyrs)=0.0d00 
             ALLOCATE(Node(i)%Elev(1:num_Lyrs+1))
             Node(i)%Elev(:)  =0.0d00 
          ENDDO
 
        !--Read coordinates of nodes
         DO i=1,num_nodes
            READ (fid_Node, *) j, Node(i)%xcord , Node(i)%ycord  ! Coordinates needs to be fixed 
            Node(i)%ID=j
            NodeList(j)=i
         ENDDO         
         
      CLOSE (fid_Node) 
 !--End Read Node input file-------------------------------------------------------------------------------------------!        
                   
      
 
 !--Read ELEMENT input file-------------------------------------------------------------------------------------------!         
    OPEN (fid_Elem, FILE='Input\'  //FileElem//'' ,STATUS='OLD')
 
         mmax=0    
         READ (fid_Elem, *) num_elem
     
         !--Allocate and Initialize Element Array
        ALLOCATE (Element(num_elem))
         DO i=1, num_elem
             Element(i)%NodeID(:)=0 
             Element(i)%NeghbrElemID(:)=0
              Element(i)%NegbrElems=0   
              Element(i)%nLP=0   
              Element(i)%Area=0.d00   
              Element(i)%UVel=0.d00    
              Element(i)%Vvel=0.d00 
              ALLOCATE(Element(i)%AvgLyrTopElv(num_Lyrs+1))
             Element(i)%AvgLyrTopElv(1:num_Lyrs+1) =0.d00               
         ENDDO
 
     
          !--Read  node numbers associated with elements
         DO i=1,num_elem
            READ (fid_Elem, *) j,  (Element(j)%NodeID(k), k=1,max_nodes_elem)
            Element(j)%ElemID=j
            DO k=1,max_nodes_elem
                mnode=Element(j)%NodeID(k)
                IF (mnode>0) THEN
                    Node(mnode)%ElemCount=Node(mnode)%ElemCount+1
                    Node(mnode)%Elem(Node(mnode)%ElemCount )=j
               ENDIF
           ENDDO
         ENDDO 
 
         !--Sort the NodeList array
        ! CALL SORT(num_nodes, NodeList)
 
         CLOSE (fid_Elem)
 !--END Read ELEMENT input file-------------------------------------------------------------------------------------------!         
      
         !--Assign coordinates to element nodes
         DO i=1, num_elem
             IF (Element(i)%NodeID(4)==0) THEN
                   nnodes=3
             ELSE
                   nnodes=4
             ENDIF
            Element(i)%Nnodes=nnodes 
             DO j=1,nnodes
                 k=Element(i)%NodeID(j)                 
                 Element(i)%xcord(j)=Node(NodeList(k))%xcord
                 Element(i)%ycord(j)=Node(NodeList(k))%ycord
                 Element(i)%zcord(j)=Node(NodeList(k))%zcord                 
             ENDDO         
             
                 
             !--Calculate the area of elements
             Element(i)%Area=GetElemArea(Element(i)%ElemID)
             
         ENDDO
         
   
        !--Find the neighbouring elements of each element
         DO I=1, num_elem
             m=0
             DO J=1, max_nodes_elem
                 k=Element(I)%NodeID(J)
                 IF (k>0) THEN
                     DO iele=1, Node(k)%ElemCount
                          IF ( Node(k)%Elem(iele)  /=  I ) THEN
                              ii=1
                              !--Make sure the neighbouring element is not already counted
                              DO, WHILE  (Element(I)%NeghbrElemID(ii) /= Node(k)%Elem(iele)   .and.  ii < max_neghbr_elems) 
                                    ii=ii+1
                              ENDDO
                             IF (ii==max_neghbr_elems)  THEN
                                   m=m+1
                                  Element(I)%NeghbrElemID(m)=Node(k)%Elem(iele) 
                            ENDIF
                           ENDIF
                     ENDDO    
                  ENDIF
             ENDDO
             Element(I)%NegbrElems=m
         ENDDO
         
     
         


        DEALLOCATE (NodeList)
         
 
   END SUBROUTINE  Read_Mesh   
    !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
  
    
    
    
    !---------------------------------------------------------------------------------------------------------------------------------------------------------!   

    SUBROUTINE ReadLyrElev (  )
    !--Purpose: Read layer elevation data at each node using C2VSimFG hydrostart input file
    !--In C2VSimFG frist column is the ground elevation and then the model cell areas and  thicknessess are provided 
    !--in subsequent columns
    !--Also, calculates the avearage elevation of layers in each element
    !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
  
    USE AdvDispLagrPar_Mod
    REAL (KIND=IKIND), ALLOCATABLE   :: vec(:)
     
     INTEGER                          :: nnds, iele, fid_lyr
     REAL (KIND=IKIND)   :: xcord, ycord
     
     CHARACTER(len=10) :: filenum
     
     fid_lyr=5050
    
     OPEN (fid_LElev, FILE='Input\'  //FileLayElev//'' ,STATUS='OLD')
          
     ALLOCATE( vec( num_Lyrs*2+1)); vec=0.0d00
 
          DO I =1, nskip_Lelev_C2VSim
                READ( fid_LElev, *) 
          ENDDO
          
          DO I=1, num_nodes
             READ (fid_LElev , *) J, (vec(K),K=1,num_Lyrs*2+1)
             Node(J)%Elev(1)=vec(1)
             DO K=2,num_Lyrs+1 
                  Node(j)%Elev(K)=Node(J)%Elev(K-1)-vec(2*K-1)
             ENDDO
          ENDDO
          
        DO K=1, num_elem
            Elem=>Element(K)
            nnds=Elem%nnodes
            DO  I=1, num_Lyrs+1
                DO J=1, nnds
                      Element(K)%AvgLyrTopElv(I)=Element(K)%AvgLyrTopElv(I)+Node(Elem%NodeID(J))%Elev(I)
                ENDDO
            ENDDO
              Element(K)%AvgLyrTopElv(:)= Element(K)%AvgLyrTopElv(:)/float(nnds)
             DO  I=1, num_Lyrs+1
                    Model%AvgLyrTElv(I)=Model%AvgLyrTElv(I)+Element(K)%AvgLyrTopElv(I)
             ENDDO
        ENDDO
        Model%AvgLyrTElv(:)=Model%AvgLyrTElv(:)/FLOAT(num_elem)
        
        !--Calculate the total volume of each layer for calculating volume weighted concentrations
        DO I=1, num_Lyrs
            DO J=1, num_elem
                Model%VolLyr(I)=Model%VolLyr(I)+ (Element(J)%AvgLyrTopElv(I)- Element(J)%AvgLyrTopElv(I+1))* Element(J)%Area
            ENDDO                
        ENDDO        
        
        !!--Map the IWFM-FEM grid on to concentration grid 
        DO  ix=1, num_xgrids	
            DO  iy=1, num_ygrids
                       xcord=ix*dx+XLwrLeft     ! xcoordinate of the cell center
                       ycord=YUprRight-iy*dy   ! ycoordinate of the cell center
                       !--Do the nearest neighbor search to find out the grid elevations
                       CALL Nearest_Elem(xcord, ycord, iele)               
                        Model%ElemID(ix,iy)=iele              
            ENDDO
        ENDDO
        !
        
        !--Output layer elevations for plotting
        do i=1, 5
            WRITE(filenum, '(I1)') i   
            OPEN (fid_lyr, FILE='Output\ElevationOfLayer_' //trim(adjustl(filenum))//trim(adjustl('.csv ')) , STATUS='UNKNOWN')
            write(fid_lyr,'(A)') 'XCORD, YCORD,  LYRELV'
            do j=1,num_nodes
                 write (fid_lyr,  '( 2(f15.2,A,2x), F8.3)'   )  Node(j)%Xcord, ",", Node(j)%Ycord, ",", Node(j)%Elev(i)
            enddo
             fid_lyr=fid_lyr+10
        enddo
        

        
          CLOSE (fid_LElev)

    
    END SUBROUTINE ReadLyrElev
    !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
    
    !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
  SUBROUTINE Nearest_Elem(xcrd, ycrd, ielem)    
  !--Purpose: Search nearest element given x and y coordinates
    !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
      USE AdvDispLagrPar_Mod
  
  INTEGER, INTENT (OUT)                     :: ielem
  REAL(KIND=IKIND) , INTENT (IN)  :: xcrd, ycrd   
  
  INTEGER                          :: enods
  REAL(KIND=IKIND)     :: elem_center_x, elem_center_y, dist_min, dist
    
 dist_min=1e10
  DO I =1, num_elem   
      enods=Element(I)%nnodes
      elem_center_x=0.0d00
      elem_center_y=0.0d00      
      DO J=1, enods 
          elem_center_x=elem_center_x+Element(I)%Xcord(J) 
          elem_center_y=elem_center_y+Element(I)%Ycord(J) 
      ENDDO
      elem_center_x=elem_center_x/float(enods)
      elem_center_y=elem_center_y/float(enods)
      
      !--Check the distance between the cell center and element center
      dist=sqrt(( elem_center_x-xcrd)**2.+ (elem_center_y-ycrd)**2.) 
      IF (dist < dist_min) THEN
          dist_min=dist
          ielem=I
       ENDIF 
      
  ENDDO
  
  
  END SUBROUTINE Nearest_Elem
      !---------------------------------------------------------------------------------------------------------------------------------------------------------!   

  
 !---------------------------------------------------------------------------------------------------------------------------------------------------------!         

   SUBROUTINE ReadPumpWell ( )
   
!---------------------------------------------------------------------------------------------------------------------------------------------------------!      
    USE AdvDispLagrPar_Mod
    
    
    
    
      OPEN (fid_PW, FILE='Input\'  //FilePumpWell//'' ,STATUS='OLD')
    
      READ (fid_PW,*)   npump_wells
      
      ALLOCATE (PumpWell(npump_wells))
      PumpWell(:)%IR=0
      PumpWell(:)%IC=0
      PumpWell(:)%nParDel=0
      
      DO I=1, npump_wells
          READ (fid_PW,*)  PumpWell(i)%Xcrd, PumpWell(i)%Ycrd
             PumpWell(i)%ID=I
           !--Figure out the cell 
            PumpWell(i)%IR=CEILING((YUprLeft-PumpWell(i)%Ycrd )/dy)
            PumpWell(i)%IC=CEILING((PumpWell(i)%Xcrd-XLwrleft)/dx)
      ENDDO
      
    END  SUBROUTINE ReadPumpWell
    
    
    
 !---------------------------------------------------------------------------------------------------------------------------------------------------------!         

   SUBROUTINE ReadIbound ( )
   
!---------------------------------------------------------------------------------------------------------------------------------------------------------!      
    USE AdvDispLagrPar_Mod
   
    INTEGER                                  :: ID, IO,nElem,ID_active
    CHARACTER (LEN=150)    :: Header
    
    OPEN (fid_Ibnd, FILE='Input\'  //FileIbnd//'' ,STATUS='OLD')
    
    ID_active=1
    READ (fid_Ibnd,*) Header 
     nElem=0
       DO
               READ (fid_Ibnd,*,IOSTAT=IO) ID, k, j,i,  ID_active
                   Model%Ibnd(i,j)=ID_active
                IF (IO/=0) EXIT
                nElem = nElem + 1
       END DO
    
       WRITE (*,*) "Number of Elements =", nElem
    
    CLOSE(fid_Ibnd)
    
    END SUBROUTINE ReadIbound 
 !---------------------------------------------------------------------------------------------------------------------------------------------------------!       
    
    
!---------------------------------------------------------------------------------------------------------------------------------------------------------!         

   SUBROUTINE Get_Nodal_Velocity ( time)
   !--Purpose: Read nodal velocities from IWFM. Velocity is written for each stress period (month) in the by IWFM. So the
   !--velocity is read in monthly intervals.
   !--IWFM calculates the Darcy velocity, so it needs to be converted Seepage Velocity (=Darcy Velocity/Porosity)-NEED TO VERYFY ?
!---------------------------------------------------------------------------------------------------------------------------------------------------------!      
     Use AdvDispLagrPar_Mod
     
   REAL (KIND=8), INTENT(IN)  :: time
   INTEGER, PARAMETER               ::  nlines=3621, nchnk=100, nvel=362148      ! Parameters from IWFM velocity output file 
   INTEGER                                          :: ID_vel, nElem, IO, ID, ID_active, kx, ky, nnodes, knode, iend
   REAL (KIND=IKIND)                   :: Xcrd, Ycrc, Zcrc, Conv_Fact, U_Vel, V_Vel, W_vel, U_max, V_max, W_max, U_min, V_min, W_min,  U_avg, V_avg, W_avg, U_abs_max, V_abs_max, W_abs_max, V_Vel_tot,U_Vel_tot
   REAL (KIND=IKIND)                   :: tmpArray(nchnk), tmpVel(nvel)
   CHARACTER (LEN=150)           :: Header, Line1, Line2
   CHARACTER (LEN=20)              ::  text3*5, text4*5, text5*5, text6*5
   
   
   

  WRITE(fid_dbug,45) 'Reading Velocity at day =', INT(time)
  WRITE(fid_dbug,*) 
  
  
   READ(fid_Vel, '(A)', end=99) Line1
   IF (Line1(1:4) /= "TEXT") THEN
       WRITE (fid_dbug, *) "Velocity read error at", time
  ENDIF
98 READ (fid_Vel, '(A)', end=99) Line2
  WRITE (*,'(100A)') Line2
  WRITE (fid_dbug , '(100A)') Line2
  
99   j=1
    DO i=1,nlines
          READ (fid_Vel, *)   tmpArray        
          tmpVel(j:j+nchnk-1)=tmpArray
          j=j+nchnk
    ENDDO
    iend=(nvel-j)
    READ (fid_Vel, *)   (tmpArray(k), k=1,iend)  
     tmpVel(j:nvel)=tmpArray(1:iend)
     
       Node(1:num_nodes)%U1=tmpVel(num_nodes*0+1: num_nodes*1)*VelConvFact/Porosity    ! Velocity converted from  ft/mon  to  m/day Pore velocity =Dacry Velocity/porosity
       Node(1:num_nodes)%V1=tmpVel(num_nodes*1+1:num_nodes*2)*VelConvFact/Porosity
       Node(1:num_nodes)%W1=tmpVel(num_nodes*2+1:num_nodes*3)*VelConvFact/Porosity
       
       Node(1:num_nodes)%U2=tmpVel(num_nodes*3+1: num_nodes*4)*VelConvFact/Porosity
       Node(1:num_nodes)%V2=tmpVel(num_nodes*4+1:num_nodes*5)*VelConvFact/Porosity
       Node(1:num_nodes)%W2=tmpVel(num_nodes*5+1:num_nodes*6)*VelConvFact/Porosity      
       Node(1:num_nodes) %U3=tmpVel(num_nodes*6+1:num_nodes*7)*VelConvFact/Porosity  
       Node(1:num_nodes)%V3=tmpVel(num_nodes*7+1:num_nodes*8)*VelConvFact/Porosity
       Node(1:num_nodes)%W3=tmpVel(num_nodes*8+1:num_nodes*9)*VelConvFact/Porosity      
       Node(1:num_nodes)%U4=tmpVel(num_nodes*9+1:num_nodes*10)*VelConvFact /Porosity 
       Node(1:num_nodes)%V4=tmpVel(num_nodes*10+1:num_nodes*11)*VelConvFact /Porosity
       Node(1:num_nodes)%W4=tmpVel(num_nodes*11+1:num_nodes*12)*VelConvFact /Porosity
       
       do I=1,num_nodes   
            Node(I)%Uvel(1)=Node(I)%U1
            Node(I)%Vvel(1)=Node(I)%V1          
            Node(I)%Wvel(1)=Node(I)%W1  
            Node(I)%Uvel(2)=Node(I)%U2
            Node(I)%Vvel(2)=Node(I)%V2          
            Node(I)%Wvel(2)=Node(I)%W2  
            Node(I)%Uvel(3)=Node(I)%U3
            Node(I)%Vvel(3)=Node(I)%V3          
            Node(I)%Wvel(3)=Node(I)%W3   
            Node(I)%Uvel(4)=Node(I)%U4
            Node(I)%Vvel(4)=Node(I)%V4          
            Node(I)%Wvel(4)=Node(I)%W4              
       end do
       


   !
   !Node(1:num_nodes)%VX_L2=tmpVel(num_nodes*3+1: num_nodes*4)
   !Node(1:num_nodes)%VY_L2=tmpVel(num_nodes*4+1:num_nodes*5)
   !Node(1:num_nodes)%VZ_L2=tmpVel(num_nodes*5+1:num_nodes*6)        
   !Node(1:num_nodes) %VX_L3=tmpVel(num_nodes*6+1:num_nodes*7)  
   !Node(1:num_nodes)%VY_L3=tmpVel(num_nodes*7+1:num_nodes*8)
   !Node(1:num_nodes)%VZ_L3=tmpVel(num_nodes*8+1:num_nodes*9)      
   !Node(1:num_nodes)%VX_L4=tmpVel(num_nodes*9+1:num_nodes*10)  
   !Node(1:num_nodes)%VY_L4=tmpVel(num_nodes*10+1:num_nodes*11)  
   !Node(1:num_nodes)%VZ_L4=tmpVel(num_nodes*11+1:num_nodes*12) 

   
    U_max=0.d00; V_max=0.d00; W_max=0.d00
    U_min=10;V_min=10; W_min=10    
    
  
    DO J=1, num_Lyrs
            U_avg=0.d00;V_avg=0.d00;W_avg=0.d00
              DO I =1, num_nodes               
                    IF (Node(i)%Uvel(J)> GW_Vel_Limit) Node(i)%Uvel(J)=GW_Vel_Limit
                    IF (Node(i)%Uvel(J)< -1*(GW_Vel_Limit))  Node(i)%Uvel(J)=-1*(GW_Vel_Limit)
                    IF (Node(i)%Vvel(J)> GW_Vel_Limit) Node(i)%Vvel(J)=GW_Vel_Limit
                    IF (Node(i)%Vvel(J)< -1*(GW_Vel_Limit))  Node(i)%Vvel(J)=-1*(GW_Vel_Limit)   
                    IF (Node(i)%Wvel(J)> GW_Vel_Limit)  Node(i)%Wvel(J)=GW_Vel_Limit
                    IF (Node(i)%Wvel(J)< -1*(GW_Vel_Limit))  Node(i)%Wvel(J)=-1*(GW_Vel_Limit)  
                    
                    !-Average velocity
                    U_avg=U_avg+Abs(Node(i)%Uvel(J))
                    V_avg=V_avg+Abs(Node(i)%Vvel(J))
                    W_avg=W_avg+Abs(Node(i)%Wvel(J))
                    
                    !-Maximum Velocity
                    IF (Node(i)%Uvel(J)>U_max) U_max=Node(i)%Uvel(J)
                    IF (Node(i)%Vvel(J)>V_max) V_max=Node(i)%Vvel(J)                    
                    IF (Node(i)%Wvel(J)>W_max) W_max=Node(i)%Wvel(J)
                    
                    !-Minimum Velocity
                    IF (Node(i)%Uvel(J)<U_min) U_min=Node(i)%Uvel(J)
                    IF (Node(i)%Vvel(J)<V_min) V_min=Node(i)%Vvel(J)                    
                    IF (Node(i)%Wvel(J)<W_min) W_min=Node(i)%Wvel(J)                    
                    
              ENDDO
              U_avg=U_avg/float(num_nodes)
              V_avg=V_avg/float(num_nodes) 
              W_avg=W_avg/float(num_nodes) 

              !U_max=MAXVAL(Node(:)%Uvel(J))
              !V_max=MAXVAL(Node(:)%Vvel(J))
              !W_max=MAXVAL(Node(:)%Wvel(J))
              !
              !U_min=MINVAL(Node(:)%Uvel(J))
              !V_min=MINVAL(Node(:)%Vvel(J))
              !W_min=MINVAL(Node(:)%Wvel(J))
              !
       
            WRITE(fid_dbug,52) 'Average U Velocity Layer :', J ,' =', U_avg, 'm/day'
            WRITE(fid_dbug,52) 'Average V Velocity Layer : ', J ,'  =', V_avg, 'm/day'
            WRITE(fid_dbug,52) 'Average W Velocity Layer : ', J ,'  =', W_avg, 'm/day'
   
            WRITE(fid_dbug,52) 'Max U Velocity Layer : ', J ,' =', U_max, 'm/day'
            WRITE(fid_dbug,52)  'Max V Velocity Layer : ', J ,'  =', V_max, 'm/day'
            WRITE(fid_dbug,52)  'Max W Velocity Layer : ', J ,'  =', W_max, 'm/day'
   
            WRITE(fid_dbug,52) 'Min U Velocity Layer : ', J ,'  =', U_min, 'm/day'
            WRITE(fid_dbug,52)  'Min V Velocity Layer : ', J ,' =', V_min, 'm/day'  
            WRITE(fid_dbug,52)  'Min W Velocity Layer : ', J ,' =', W_min, 'm/day'        
            WRITE(fid_dbug,*) 

    ENDDO  
    
    
      CALL SetTrnspTimeStep (Abs(U_max), Abs(V_max))

    !U_vel=U_avg
    !V_vel=V_avg
        
        
  ! WRITE(*,55) 'Average U Velocity =', U_avg, 'm/day'
  ! WRITE(*,55) 'Average V Velocity =', V_avg, 'm/day'
   
   !WRITE(*,50) 'Max U Velocity =', U_max, 'm/day'
  ! WRITE(*,50)  'Max V Velocity =', V_max, 'm/day'
   
  ! WRITE(*,50) 'Min U Velocity =', U_min, 'm/day'
  ! WRITE(*,50)  'Min V Velocity =', V_min, 'm/day'

    
45 FORMAT (A25,  I8 )     
50 FORMAT (A25,  F12.6,A8 )    
52 FORMAT (A25, I4, A5,  F12.6,A8 )
55 FORMAT (A25, I4, F12.6, A8 )
   
    
!    Vtot=sqrt(U_avg*U_avg+V_avg*V_avg)
!    DispCoeff_xx=LongDisp*(U_avg*U_avg)/Vtot+HorTrnsDisp*(V_avg*V_avg)/Vtot   ! (m2/day)
 !   DispCoeff_yy=HorTrnsDisp*(V_avg*V_avg)/Vtot+HorTrnsDisp*(U_avg*U_avg)/Vtot   ! (m2/day)
    
    !--Assign Dispersion Coefficients
    !Dxx=DispCoeff_xx
    !Dyy=DispCoeff_yy
    !Dyx=0.0d00
    !Dxy=0.0d00   
 
    !--Calculates Grid Peclet Number
!    Pex=abs(U_avg)*dx/DispCoeff_xx
!    Pey=abs(V_avg)*dy/DispCoeff_yy   
    
    ! WRITE(*,*)
    ! WRITE(*,'(A30, F8.3)') 'Average Grid Peclet Number=', max(Pex,Pey)  
    !WRITE(*,*)
    !
     !WRITE(fid_dbug,*)
     !WRITE(fid_dbug,'(A30, F8.3)') 'Average Grid Peclet Number=', max(Pex,Pey)  
     !WRITE(fid_dbug,*)
     
     
     !!--Calculate Element Velocity
     !DO I =1, num_elem
     !     Element(i)%Uvel=0.0d00
     !     Element(i)%Vvel=0.0d00
     !      nnodes=Element(I)%Nnodes
     !       DO J=1, nnodes
     !             knode=Element(i)%NodeID(J)
     !             Element(i)%Uvel=Element(i)%Uvel+Node( knode)%Uvel
     !             Element(i)%Vvel=Element(i)%Vvel+Node( knode)%Vvel
     !       ENDDO
     !       !-Element velocity is calculate as the average nodal velocity
     !       Element(i)%Uvel=Element(i)%Uvel/float(nnodes)
     !       Element(i)%Vvel=Element(i)%Vvel/float(nnodes)            
     !ENDDO

   END SUBROUTINE Get_Nodal_Velocity
!---------------------------------------------------------------------------------------------------------------------------------------------------------!         
    
    
    
    
    
!---------------------------------------------------------------------------------------------------------------------------------------------------------!         

   SUBROUTINE Get_Cell_Velocity_SS ( )
   
!---------------------------------------------------------------------------------------------------------------------------------------------------------!      
    USE AdvDispLagrPar_Mod
   !REAL (KIND=IKIND), INTENT (IN) :: time
   INTEGER                                  :: ID_vel, nElem, IO, ID, ID_active
   REAL (KIND=8)                    :: W_vel, Conv_Fact, U_Vel, V_Vel,U_max, V_max, U_min, V_min, U_avg, V_avg,U_abs_max, V_abs_max, V_Vel_tot,U_Vel_tot
   CHARACTER (LEN=150)   :: Header
   
   OPEN (fid_Vel, FILE='Input\'  //FileVel//'' ,STATUS='OLD')
   
     !--Read the FEM mesh (nodes and elements) first
      CALL Read_Mesh ( )

   
   READ (fid_Vel, *) ID_vel
   
   IF (ID_vel==1) THEN 
          READ (fid_Vel,*) U_Vel	                        !Stream velocity (seepage v/n) in  x direction (u) in m/day 
          READ (fid_Vel,*) V_Vel	                         !Stream velocity (seepage v/n) in  y direction (u) in m/day
          READ (fid_Vel,*) W_Vel	                         !Stream velocity (seepage v/n) in z direction (u) in m/day
         U_max=U_Vel
         V_max=V_Vel
         U_min=U_Vel
         V_min=V_Vel
         U_avg=U_Vel
         V_avg=V_Vel
         
        !Model%Uvel(:,:)=U_Vel
       ! Model%Vvel(:,:)=V_Vel         

         
   ELSE
         READ (fid_Vel,*) Conv_Fact
         READ (fid_Vel,*) Header                         ! Read Header

         U_max=0.d00
         V_max=0.d00
         U_min=100000
         V_min=100000
         U_avg=0.d00
         V_avg=0.d00
         
         nElem=0
         DO
               READ (fid_Vel,*,IOSTAT=IO) ID, k, j, i, U_vel, V_vel, W_vel, ID_active
               
               U_vel=U_vel*Conv_Fact      ! m/day
               V_vel=V_vel*Conv_Fact      ! m/day
               W_vel=W_vel*Conv_Fact    ! m/day
               
               !-Cap the velocity at predefined value ( 0.5 m/day)
              IF (U_vel> GW_Vel_Limit) U_vel=GW_Vel_Limit
               IF (U_vel< -1*(GW_Vel_Limit))  U_vel=-1*(GW_Vel_Limit)
               IF (V_vel> GW_Vel_Limit) V_vel=GW_Vel_Limit
               IF (V_vel< -1*(GW_Vel_Limit))  V_vel=-1*(GW_Vel_Limit)               
               
               IF (U_vel>U_max) U_max=U_vel
               IF (V_vel>V_max) V_max=V_vel
               IF (U_vel<U_min) U_min=U_vel
               IF (V_vel<V_min)V_min=V_vel
               U_avg=U_avg+U_vel
               V_avg=V_avg+V_vel
               
                Model%Uvel(i,j)=U_vel
               Model%Vvel(i,j)=V_vel

               
                IF (IO/=0) EXIT
                nElem = nElem + 1
         END DO
           
        U_avg=U_avg/float(nElem)
        V_avg=V_avg/float(nElem) 
        
        U_vel=U_avg
        V_vel=V_avg
         
   ENDIF

   
   WRITE(*,50) 'Average U Velocity =', U_avg, 'm/day'
   WRITE(*,50) 'Average V Velocity =', V_avg, 'm/day'
   
   WRITE(*,50) 'Max U Velocity =', U_max, 'm/day'
   WRITE(*,50)  'Max V Velocity =', V_max, 'm/day'
   
   WRITE(*,50) 'Min U Velocity =', U_min, 'm/day'
   WRITE(*,50)  'Min V Velocity =', V_min, 'm/day'
   
    WRITE(fid_dbug,50) 'Average U Velocity =', U_avg, 'm/day'
    WRITE(fid_dbug,50) 'Average V Velocity =', V_avg, 'm/day'
   
    WRITE(fid_dbug,50) 'Max U Velocity =', U_max, 'm/day'
    WRITE(fid_dbug,50)  'Max V Velocity =', V_max, 'm/day'
   
    WRITE(fid_dbug,50) 'Min U Velocity =', U_min, 'm/day'
    WRITE(fid_dbug,50)  'Min V Velocity =', V_min, 'm/day'  
    
50 FORMAT (A25, F6.3,A8 )


   
   
    HorTrnsDisp=VerTrns2LongDispRat*LongDisp
    
    
    Vtot=sqrt(U_avg*U_avg+V_avg*V_avg)
    DispCoeff_xx=LongDisp*(U_avg*U_avg)/Vtot+HorTrnsDisp*(V_avg*V_avg)/Vtot   ! (m2/day)
    DispCoeff_yy=HorTrnsDisp*(V_avg*V_avg)/Vtot+HorTrnsDisp*(U_avg*U_avg)/Vtot   ! (m2/day)
    
    !--Assign Dispersion Coefficients
    Dxx=DispCoeff_xx
    Dyy=DispCoeff_yy
    Dyx=0.0d00
    Dxy=0.0d00   
    
           
        
    !--Calculates Grid Peclet Number
    Pex=abs(U_avg)*dx/DispCoeff_xx
    Pey=abs(V_avg)*dy/DispCoeff_yy   
    
     WRITE(*,*)
     WRITE(*,'(A30, F8.3)') 'Average Grid Peclet Number=', max(Pex,Pey)  
    WRITE(*,*)
    
     WRITE(fid_dbug,*)
     WRITE(fid_dbug,'(A30, F8.3)') 'Average Grid Peclet Number=', max(Pex,Pey)  
     WRITE(fid_dbug,*)

    
    !  !--Calculates Grid Peclet Number
    !Pe=sqrt(U_Vel*U_Vel+V_Vel*V_Vel)*sqrt(dx*dx+dy*dy)/sqrt(DispCoeff_xx*DispCoeff_xx+DispCoeff_yy*DispCoeff_yy)
    !Pex=U_vel*dx/DispCoeff_xx
    !Pey=V_vel*dy/DispCoeff_yy
    !
    !write (*,*) 'Grid Peclet Number=', max(Pex,Pey)
    
    BB=2*DispCoeff_xx/U_Vel
    
   !--Chek Courant Condition and adjust the time step if necessary
    U_abs_max= max(U_max, abs(U_min))
    V_abs_max= max(V_max, abs(V_min))    
    CALL SetTrnspTimeStep (U_abs_max, V_abs_max )
    
    !if (del_t>abs(dx/U_max)  .or.   del_t>abs(dy/V_max)) then
    !    write(*,*) 'Cournt condition is not met'
    !    !del_t=del_t/2
    !        
    !   CALL SetTrnspTimeStep ( )
    !endif


    CLOSE (fid_Vel)

    END SUBROUTINE Get_Cell_Velocity_SS
 !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
    
    
 !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
    SUBROUTINE     SetTrnspTimeStep (Umax, Vmax )
 !---------------------------------------------------------------------------------------------------------------------------------------------------------!   
     USE AdvDispLagrPar_Mod
    
    REAL (KIND=8), INTENT (IN)      ::  Umax, Vmax
    REAL (KIND=8)      :: dt_sigma, dt_nue, dt_crit
    LOGICAL         :: IDT  

            dt=del_t
            !--Check the Advection Time Step Criteria (CFL)
            IDT=.FALSE.
            if (dt .LE. dx/abs(Umax)   .and.   dt .LE. dy/abs(Vmax)) then
                WRITE(*,*) '!----CFL  condition is satisfied. Transport time step is not adjusted.'
                WRITE(fid_dbug,*)  '!----CFL  condition is satisfied. Transport time step is not adjusted.'
            else 
                write(*,*) '!----CFL condition is not satisfied. Transport time step is being adjusted.'
               dt =min(dx/abs(Umax),dy/abs(Vmax))
               write(*, 160)' New transport time step is =', dt
               WRITE(fid_dbug,160) ' New transport time step is =', dt
               IDT=.TRUE.
            endif
    
            !if ( IDT) then
            !    num_time_steps=nint(t_sim/dt)				                        !No of time steps
            !    nt_skip=nint(t_sim/(dt*num_out_time_steps))             !Keep the same output steps but change the skiping step	  
            !      write(*,*) ' !----Adjusted timestep skip for output is =', nt_skip
            !    ! write(*,'(A45, I6 )') ' !----Adjusted timestep skip for output is =', nt_skip
            !endif
    
          !--Check for the Dispersion Time Step Criteria 
          if ( dt .LE. 0.5*dx*dx/DispCoeff_xx .and. dt .LE.0.5*dy*dy/DispCoeff_yy) then
               write(*,*) '!----Dispersion time criteria is satisfied. Transport time step is not adjusted.'
               WRITE(fid_dbug,*) '!----Dispersion time criteria is satisfied. Transport time step is not adjusted.'
            else 
               WRITE(*, *) '!----Dispersion time criteria is not satisfied. Transport time step is being adjusted.'
               dt =min(0.5*dx*dx/DispCoeff_xx  , 0.5*dy*dy/DispCoeff_yy)
               WRITE(*, 160)'Adjusted transport time step is =', dt
               WRITE(fid_dbug,160) ' Adjusted transport time step is =', dt
               WRITE(fid_dbug,*)
               WRITE(*, *)
            endif
    
            if ( IDT) then
                del_t=dt
                num_time_steps=nint(t_sim/del_t)				                        !No of time steps
                nt_skip=nint(t_sim/(del_t*num_out_time_steps))             !Keep the same output steps but change the skiping step	  
                !num_parcels_dt=
                  WRITE(*, 150)                ' Adjusted timestep skip for output is =', nt_skip
                  WRITE(fid_dbug,150) ' Adjusted timestep skip for output is =', nt_skip
                  WRITE(fid_dbug,*)
                  WRITE(*, *)
                ! write(*,'(A45, I6 )') ' !----Adjusted timestep skip for output is =', nt_skip   
                  
                      
                WRITE (*,160) 'Adjusted time step is =', del_t
                WRITE (*,150) 'Adjusted number of time steps=', num_time_steps
               ! WRITE (*,150) 'Adjusted of parcels released per day=', num_parcels_dt
                WRITE (*,150) 'Adjustedtime step skip for output is=', nt_skip
                WRITE(*, 150)
            
                WRITE(fid_dbug,160) 'Adjusted time step is =', del_t
                WRITE(fid_dbug,150) 'Adjusted number of time steps=', num_time_steps
                !WRITE (fid_dbug,150)'Adjusted of parcels released per day=', num_parcels_dt
                WRITE(fid_dbug,150) 'Adjusted time step skip for output is=', nt_skip
            WRITE(fid_dbug,150)
            
            endif
            
            !--Time step for stable solution
            !--Zienkiewicz and Taylor, The Finite Element Method, Fith Edition, Volume 3: Fluid Dynamics

            !dt_sigma=min( dx/abs(Umax) , dy/abs(Vmax) )
            !dt_nue=dx**2/(2*DispCoeff_xx)
            !
            !dt_crit=dt_sigma*dt_nue/(dt_sigma+dt_nue)
            !
            !IF (dt >dt_crit) THEN
            !     dt=dt_crit
            !     write (*,*) ' Time step is greater than critical value for stability. Switiching to critical time step'
            !ENDIF
        
      !!--Update the time step      
      !      del_t=dt    
      !      num_time_steps=nint(t_sim/del_t)	
      !      nt_skip=nint(t_sim/(del_t*num_out_time_steps))
      !      !num_parcels_dt=num_parcels/(t_spill)  
        
150         FORMAT (A50, I12 )  
160         FORMAT (A50, F8.3 )  
            
    
    END  SUBROUTINE   SetTrnspTimeStep      
     !---------------------------------------------------------------------------------------------------------------------------------------------------------!   

    
    
!  !---------------------------------------------------------------------------------------------------------------------------------------------------------!         
!
   SUBROUTINE Get_Nodal_Velocity_SS ( )
   !--Purpose: Read nodal velocities from IWFM. Velocity is written for each stress period (month) in the by IWFM. So the
   !--velocity is read in monthly intervals.
!---------------------------------------------------------------------------------------------------------------------------------------------------------!      
     Use AdvDispLagrPar_Mod
     
   !REAL (KIND=8), INTENT(IN)  :: time
   INTEGER                                          :: ID_vel, nElem, IO, ID, ID_active, kx, ky, nnodes, knode
   REAL (KIND=8)                            :: Xcrd, Ycrc, W_vel, Conv_Fact, U_Vel, V_Vel,U_max, V_max, U_min, V_min, U_avg, V_avg,U_abs_max, V_abs_max, V_Vel_tot,U_Vel_tot
   CHARACTER (LEN=150)           :: Header
   
   OPEN (fid_Vel, FILE='Input\'  //FileVel//'' ,STATUS='OLD')
   !
   !
   !!--Read the FEM mesh (nodes and elements) first
  !    CALL Read_Mesh ( )

    READ (fid_Vel,*) Conv_Fact
    READ (fid_Vel,*) Header                         ! Read Header
   
    U_max=0.d00
    V_max=0.d00
    U_min=100000
    V_min=100000
    U_avg=0.d00
    V_avg=0.d00
         
    nElem=0
    DO
            READ (fid_Vel,*,IOSTAT=IO) I, Xcrd, Ycrc, U_vel, V_vel, W_vel
        
            
            U_vel=U_vel*Conv_Fact/Porosity
            V_vel=V_vel*Conv_Fact/Porosity
            
               !-Cap the velocity at predefined value ( 2 m/day)
              IF (U_vel> GW_Vel_Limit) U_vel=GW_Vel_Limit
               IF (U_vel< -1*(GW_Vel_Limit))  U_vel=-1*(GW_Vel_Limit)
               IF (V_vel> GW_Vel_Limit) V_vel=GW_Vel_Limit
               IF (V_vel< -1*(GW_Vel_Limit))  V_vel=-1*(GW_Vel_Limit)               
              ! 
               IF (U_vel>U_max) U_max=U_vel
               IF (V_vel>V_max) V_max=V_vel
               IF (U_vel<U_min) U_min=U_vel
               IF (V_vel<V_min)V_min=V_vel      
            
            Node(I)%UVel=U_vel
            Node(I)%VVel=V_vel
            Node(I)%WVel=W_vel
        
            IF (U_vel>U_max) U_max=U_vel
            IF (V_vel>V_max) V_max=V_vel
            IF (U_vel<U_min) U_min=U_vel
            IF (V_vel<V_min)V_min=V_vel
            U_avg=U_avg+U_vel
            V_avg=V_avg+V_vel
        
            IF (IO/=0) EXIT
            nElem = nElem + 1
    END DO
           
    U_avg=U_avg/float(nElem)
    V_avg=V_avg/float(nElem) 
        
    U_vel=U_avg
    V_vel=V_avg
        
        
   WRITE(*,55) 'Average U Velocity =', U_avg, 'm/day'
   WRITE(*,55) 'Average V Velocity =', V_avg, 'm/day'
   
   WRITE(*,50) 'Max U Velocity =', U_max, 'm/day'
   WRITE(*,50)  'Max V Velocity =', V_max, 'm/day'
   
   WRITE(*,50) 'Min U Velocity =', U_min, 'm/day'
   WRITE(*,50)  'Min V Velocity =', V_min, 'm/day'

    WRITE(fid_dbug,50) 'Average U Velocity =', U_avg, 'm/day'
    WRITE(fid_dbug,50) 'Average V Velocity =', V_avg, 'm/day'
   
    WRITE(fid_dbug,50) 'Max U Velocity =', U_max, 'm/day'
    WRITE(fid_dbug,50)  'Max V Velocity =', V_max, 'm/day'
   
    WRITE(fid_dbug,50) 'Min U Velocity =', U_min, 'm/day'
    WRITE(fid_dbug,50)  'Min V Velocity =', V_min, 'm/day'  
    
50 FORMAT (A25, F12.6,A8 )
55 FORMAT (A25, F12.6,A8 )
   
    HorTrnsDisp=VerTrns2LongDispRat*LongDisp
    
    
    Vtot=sqrt(U_avg*U_avg+V_avg*V_avg)
    DispCoeff_xx=LongDisp*(U_avg*U_avg)/Vtot+HorTrnsDisp*(V_avg*V_avg)/Vtot   ! (m2/day)
    DispCoeff_yy=HorTrnsDisp*(V_avg*V_avg)/Vtot+HorTrnsDisp*(U_avg*U_avg)/Vtot   ! (m2/day)
    
    !--Assign Dispersion Coefficients
    Dxx=DispCoeff_xx
    Dyy=DispCoeff_yy
    Dyx=0.0d00
    Dxy=0.0d00   
 
    !--Calculates Grid Peclet Number
    Pex=abs(U_avg)*dx/DispCoeff_xx
    Pey=abs(V_avg)*dy/DispCoeff_yy   
    
    ! WRITE(*,*)
    ! WRITE(*,'(A30, F8.3)') 'Average Grid Peclet Number=', max(Pex,Pey)  
    !WRITE(*,*)
    !
     WRITE(fid_dbug,*)
     WRITE(fid_dbug,'(A30, F8.3)') 'Average Grid Peclet Number=', max(Pex,Pey)  
     WRITE(fid_dbug,*)
     
     
     !!--Calculate Element Velocity
     !DO I =1, num_elem
     !     Element(i)%Uvel=0.0d00
     !     Element(i)%Vvel=0.0d00
     !      nnodes=Element(I)%Nnodes
     !       DO J=1, nnodes
     !             knode=Element(i)%NodeID(J)
     !             Element(i)%Uvel=Element(i)%Uvel+Node( knode)%Uvel
     !             Element(i)%Vvel=Element(i)%Vvel+Node( knode)%Vvel
     !       ENDDO
     !       !-Element velocity is calculate as the average nodal velocity
     !       Element(i)%Uvel=Element(i)%Uvel/float(nnodes)
     !       Element(i)%Vvel=Element(i)%Vvel/float(nnodes)            
     !ENDDO
     
   REWIND(fid_Vel)
   END SUBROUTINE Get_Nodal_Velocity_SS
!---------------------------------------------------------------------------------------------------------------------------------------------------------!          
!    
!    
    
   !!--------------------------------------------------------------------------------------------------------!       
   !SUBROUTINE Read_Mesh ( )
   !!--------------------------------------------------------------------------------------------------------!    
   !
   !USE AdvDispLagrPar_Mod
   !
   !
   !       INTEGER    :: mmax, m, nnodes,num_nodes_all, n,mnode
   !       INTEGER, ALLOCATABLE :: NodeList(:)
   !       
   !       
   !!--Read Node input file-------------------------------------------------------------------------------------------!    
   !    OPEN (fid_Node,FILE='..\Input\' //FileNode//'' ,STATUS='OLD')
   !  
   !      READ (fid_Node, *) num_nodes
   !      READ (fid_Node, *)
   !  
   !      !--Allocate and Initialize Node array
   !      ALLOCATE (Node(num_nodes))
   !     ALLOCATE (NodeList(num_nodes*4))
   !     NodeList(:)=0
   !       DO i=1, num_nodes
   !          Node(i)%ID=0
   !          Node(i)%Elem(:)=0
   !          Node(i)%ElemCount=0
   !          Node(i)%xcord=0.0d00
   !          Node(i)%ycord=0.0d00   
   !          Node(i)%Conc=0.0d00  
   !          Node(i)%ConcAna=0.0d00 
   !          Node(i)%Vel=0.0d00 
   !          Node(i)%Area=0.0d00 
   !       ENDDO
   !
   !     !--Read coordinates of nodes
   !      DO i=1,num_nodes
   !         READ (fid_Node, *) j, Node(i)%xcord , Node(i)%ycord  ! Coordinates needs to be fixed 
   !         Node(i)%ID=j
   !         NodeList(j)=i
   !      ENDDO         
   !      
   !   CLOSE (fid_Node) 
 !--End Read Node input file-------------------------------------------------------------------------------------------!        
                   
      
 
 !!--Read ELEMENT input file-------------------------------------------------------------------------------------------!         
 !   OPEN (fid_Elem, FILE='..\Input\'  //FileElem//'' ,STATUS='OLD')
 !
 !        mmax=0    
 !        READ (fid_Elem, *) num_elem
 !    
 !        !--Allocate and Initialize Element Array
 !       ALLOCATE (Element(num_elem))
 !        DO i=1, num_elem
 !            Element(i)%NodeID(:)=0 
 !        ENDDO
 !
 !    
 !         !--Read  node numbers associated with elements
 !        DO i=1,num_elem
 !           READ (fid_Elem, *) j,  (Element(j)%NodeID(k), k=1,max_nodes_elem)
 !           DO k=1,max_nodes_elem
 !               mnode=Element(j)%NodeID(k)
 !               IF (mnode>0) THEN
 !                   Node(mnode)%ElemCount=Node(mnode)%ElemCount+1
 !                   Node(mnode)%Elem(Node(mnode)%ElemCount )=j
 !              ENDIF
 !          ENDDO
 !        ENDDO 
 !
 !        !--Sort the NodeList array
 !       ! CALL SORT(num_nodes, NodeList)
 !
 !        CLOSE (fid_Elem)
 !!--END Read ELEMENT input file-------------------------------------------------------------------------------------------!         
 !     
 !        !--Assign coordinates to element nodes
 !        DO i=1, num_elem
 !            IF (Element(i)%NodeID(4)==0) THEN
 !                  nnodes=3
 !            ELSE
 !                  nnodes=4
 !           ENDIF
 !            DO j=1,nnodes
 !                k=Element(i)%NodeID(j)                 
 !                Element(i)%xcord(j)=Node(NodeList(k))%xcord
 !                Element(i)%ycord(j)=Node(NodeList(k))%ycord
 !           ENDDO
 !        ENDDO
         
    
!        DEALLOCATE (NodeList)
         
 
!  END SUBROUTINE  Read_Mesh   
   !--------------------------------------------------------------------------------------------------------!         
