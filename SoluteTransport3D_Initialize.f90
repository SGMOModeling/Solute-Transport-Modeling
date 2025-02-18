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

! *****************************************************************************
    SUBROUTINE Initialize ( ) 
!--Purpose: Initialize variables, allocates arrays and calculates some intial 
!   parameter values  
! *****************************************************************************
    
    
    USE AdvDispLagrPar_Mod
     USE LinkedList_Mod

IMPLICIT NONE

    
    INTEGER ::jk, ielem, ix, iy, iz, igrid, ngrid
    REAL (KIND=IKIND) :: xcrd, ycrd, zcrd
    CHARACTER (LEN=150) :: Header, text1, text2, OutHeader*20
    

    
 !-- FOR MULTIPLE RELEASES WITH DIFFERENT RATES/DURATIONS THIS NEEDS TO BE MODIFIED
    
    num_time_steps=int((t_sim/del_t))		                                     !No of time steps (will get updated if the time step chaged for stability criteria)
	num_out_time_steps=nint(t_sim/(del_t*nt_skip))      	     !No of output time steps (will get updated if the time step chaged for stability criteria)

    
	num_ygrids=nint((YUprRight-YLwrLeft)/dx)				                     !No of x-grids
	num_xgrids=nint((XUprRight-XLwrLeft ) /dy)				                     !No of y-grids
    num_zgrids=num_Lyrs                                                                         !--No of z-grids
    num_nodes_conc=(num_xgrids+1)*(num_ygrids+1)     ! Number of nodes for calculating concentrations
    

    mass_of_parcel=TotalMass/float(num_parcels)                  !Mass of a parcel (kg/day per location/per depth)  -MassRate is in mass per unit lenght (lenght taken as dz)    


    WRITE(fid_dbug,*) 
    WRITE(fid_dbug,*) "Time related input data provided "
    WRITE(fid_dbug,*) "These maybe adjusted later to satify the stability criteria "
    WRITE(fid_dbug,160) "Time step specified=", del_t
    WRITE(fid_dbug,150) " Number of time steps specified=",  num_time_steps
    WRITE(fid_dbug,150) " Number of output time steps specified=",num_out_time_steps
    WRITE(fid_dbug,160) " Mass of a parcel specified=",mass_of_parcel

    WRITE(fid_dbug,*) 
    
    
    DO I=1, num_sources
        Source(I)%npar_rel=int(Source(I)%MassRate*Source(i)%t_spill/TotalMass*num_parcels)
        Source(I)%npar_dt=int( Source(I)%npar_rel/Source(i)%t_spill*del_t)
        !-Update the number of parcels to conserve mass
        !Source(I)%npar_rel=Source(I)%npar_dt*del_t*Source(I)%MassRate/TotalMassRate*num_time_steps 
        Source(I)%tot_mass=Source(I)%npar_dt*del_t*int((Source(I)%t_spill/del_t))
     
        WRITE(fid_dbug,150) " Number parcels released per dt at location:",I, "=",Source(I)%npar_dt

    ENDDO
    WRITE(fid_dbug,*)

150 FORMAT (A50, I2, A2, I12 )   
160 FORMAT (A50, F12.6 )  
    
!--ALLOCATE AND INITIALIZE ARRAYS     

 	    ALLOCATE (Model%ConAnal(1:num_xgrids,1:num_ygrids, 1:num_zgrids));Model%ConAnal(:,:,:) =0.00      ! (col, row)              
        ALLOCATE (Model%ElemID(1:num_xgrids, 1:num_ygrids)); Model%ElemID(:,:) =0 
        ALLOCATE (Model%Ibnd(1:num_xgrids+2, 1:num_ygrids+2)); Model%Ibnd(:,:) =1
        ALLOCATE (Model%Npar(1:num_xgrids, 1:num_ygrids+2, 1:num_zgrids)); Model%Npar(:,:, :) =0.00
        ALLOCATE (Model%AvgLyrTElv(1:num_Lyrs+1)); Model%AvgLyrTElv(:) =0.00
        ALLOCATE (Model%VolLyr(1:num_Lyrs)); Model%VolLyr(:) =0.00       
     
!--Prepare node concentrations for output

    !--Allocate and Initialize Node array
    !ALLOCATE (NodeC(num_nodes_conc))
    !DO i=1, num_nodes_conc
    !    NodeC(i)%ID=0
    !    NodeC(i)%xcord=0.0d00
    !    NodeC(i)%ycord=0.0d00    
    !    NodeC(i)%Conc=0.0d00  
    !    NodeC(i)%ConcAna=0.0d00  
    !ENDDO
   
   ALLOCATE(Link(0:num_parcels+1));   Link(:)=0

  !--Create a Linked List for LPs
  CALL LPList%Initialize( )

  CALL Make_List(num_parcels)

     !-- Initialize the particle grid particle numbers
 	Model%npar(1:num_xgrids, 1:num_ygrids, 1:num_zgrids)=0
   
    
    !--Read the FEM mesh (nodes and elements) first
      CALL Read_Mesh ( )   
    !  --READ LAYER ELEVATION DATA 
    IF (idz_ReadFile)  CALL ReadLyrElev ( ) 
  
  
  !--Open file for reading velocity data   and read initial data
     OPEN (fid_Vel, FILE='Input\'  //FileVel//'' ,STATUS='OLD')
        DO i=1, iskip_beg_Fvel
              READ (fid_Vel, *) Header
        ENDDO
    
        DO i=1, num_nodes
            READ (fid_Vel, *) text1     ! Node(i)%xcord, Node(i)%ycord
        ENDDO          
    
    !read (infile,*) ajunk
      DO i=1, nnode_con
            READ (fid_Vel, *)  text2     ! Nconn(i)%ID, Nconn(i)%CN1,Nconn(i)%CN2,Nconn(i)%CN3
      ENDDO
    
   !--Prepare for output
 ngrid=num_xgrids*num_ygrids*num_zgrids

   ALLOCATE (Grid(ngrid)) 
      !--Open File for output
     OPEN (fid_ouput_Node_All, FILE=  'Output\ConLagrng_3D_AllInOne.dat ' , STATUS='REPLACE')
     
        igrid=1
        DO iz=1, num_zgrids
            DO  ix=1, num_xgrids	
                DO  iy=1, num_ygrids
                    xcrd=ix*dx+XLwrLeft
                    ycrd=YUprRight-iy*dy
                    !zcrd=dz/2+(iz-1)*dz+ZLwr
                        ielem= Model%ElemID(ix,iy)
                        zcrd=0.5*(Element(ielem)%AvgLyrTopElv(iz)+Element(ielem)%AvgLyrTopElv(iz+1))
                        Grid(igrid)%Xcrd =xcrd
                        Grid(igrid)%Ycrd= ycrd
                        Grid(igrid)%Zcrd= zcrd
                        Grid(igrid)%Ix=ix
                        Grid(igrid)%Iy=iy
                        Grid(igrid)%Iz=iz                        
                        igrid=igrid+1
                        !WRITE (fid_ouput_Node, '( 2(f20.2,1x),f15.2,1x, 3I6, I10, 3(f15.2,1x), I10,  (f15.3,1x)  )')  xcrd, ycrd, zcrd, ix, iy, iz, Element(ielem)%ElemID, Element(ielem)%AvgLyrTopElv(iz), Element(ielem)%AvgLyrTopElv(iz+1),   elem_thick, Model%npar(ix,iy,iz), Model%npar(ix,iy,iz)*mass_of_parcel/(dx*dy*elem_thick)/Porosity*1000  ! g/m3=ppm
                ENDDO
            ENDDO
    ENDDO
    
    WRITE (fid_ouput_Node_All, '(A15)', ADVANCE='NO') "XCOORD"
     WRITE (fid_ouput_Node_All, '(<ngrid>(F12.2))', ADVANCE='NO') (Grid(igrid)%Xcrd, igrid=1,ngrid)
     WRITE(fid_ouput_Node_All,*)
    WRITE (fid_ouput_Node_All, '(A15)', ADVANCE='NO') "YCOORD"
     WRITE (fid_ouput_Node_All, '(<ngrid>(F12.2))', ADVANCE='NO') (Grid(igrid)%Ycrd, igrid=1,ngrid)
     WRITE(fid_ouput_Node_All,*)
    WRITE (fid_ouput_Node_All, '(A15)', ADVANCE='NO') "ZCOORD"
     WRITE (fid_ouput_Node_All, '(<ngrid>(F12.2))', ADVANCE='NO') (Grid(igrid)%Zcrd, igrid=1,ngrid)
     WRITE(fid_ouput_Node_All,*)
    WRITE (fid_ouput_Node_All, '(A15)', ADVANCE='NO') "ROW"
     WRITE (fid_ouput_Node_All, '(<ngrid>(I12))', ADVANCE='NO') (Grid(igrid)%Ix, igrid=1,ngrid)
     WRITE(fid_ouput_Node_All,*)
    WRITE (fid_ouput_Node_All, '(A15)', ADVANCE='NO') "COLUMN"
     WRITE (fid_ouput_Node_All, '(<ngrid>(I12))', ADVANCE='NO') (Grid(igrid)%Iy, igrid=1,ngrid)
     WRITE(fid_ouput_Node_All,*)
     WRITE (fid_ouput_Node_All, '(A15)', ADVANCE='NO') "LAYER"
     WRITE (fid_ouput_Node_All, '(<ngrid>(I12))', ADVANCE='NO') (Grid(igrid)%Iz, igrid=1,ngrid)
     WRITE(fid_ouput_Node_All,*)    
     
     DEALLOCATE (Grid)
     
    END SUBROUTINE Initialize 
    ! *****************************************************************************
