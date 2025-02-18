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
     
    
    MODULE AdvDispLagrPar_Mod
    
    INTEGER, PARAMETER :: PI=3.14159265358979D0 , IKIND=8

    INTEGER, PARAMETER ::fid_dbug=50, fid_Data=110, fid_vel=105, fid_vel_info=102,fid_Ibnd=115, fid_Geninf=120,  fid_ouput=130, fid_ouput_Anal=150, fid_Elem=160,  fid_Node=170
     INTEGER, PARAMETER :: fid_ouput_Node_start=1000, fid_out_XYplane_Start=2000, fid_out_YZplane_Start=3000,fid_out_XZplane_Start=4000,  fid_PW=170, fid_LElev=180
     INTEGER, PARAMETER ::  fid_ouput_Node_XY=1010, fid_ouput_Node_XZ=1020,fid_ouput_Node_YZ=1030, fid_ouput_Node_All=1500
    INTEGER, PARAMETER :: NMAX=256  !NUMBER OF TERMS USED IN GAUSS-LEGENDRE NUMERICAL INTEGRATION TECHNIQUE (MUST EQAUL 4, 20, 60,104, or 256)
     INTEGER, PARAMETER ::num_par_max=5000000, n_partype=2, nodes_per_elem=4, max_nodes_elem=4, max_neghbr_elems=14, iskip_beg_Fvel=3, nnode_con=37171 !, num_Lyrs=40
    INTEGER, PARAMETER  :: nskip_Lelev_C2VSim=105
    REAL, PARAMETER        ::  GW_Vel_Limit =10.0,    fact_bnd_offset=0.2 !(m/mon)  !ConvFactFtMtrs-converstion factor feet to meters


     
     
    INTEGER                                 :: i, j, k,  ih, jt, iv, ip, ind, TT, num_parcels, num_parcels_dt,  num_par_t, num_nodes, num_elem, num_par_rel, num_par, num_sources, ReleaseType, IConcFlg, iflg_vel, num_Lyrs, iFlag_Output
    INTEGER                                 :: num_out_time_steps,num_out_grids, num_xgrids, num_ygrids, num_zgrids, num_time_steps, nx_skip,nt_skip, ind_outgridx, ind_outgridy, npump_wells, num_nodes_conc
    INTEGER, ALLOCATABLE  :: Link(:)
    INTEGER                                 :: beg(n_partype),ID_pw, NDaysMonCum(12)
    REAL (KIND=IKIND)         :: Mass, MassRate, Porosity, DispCoeff_xx, DispCoeff_yy, DispCoeff_zz, DispCoeff_xy, DispCoeff_xz, DispCoeff_yz, LongDisp, HorTrnsDisp, VerTrnsDisp,VerTrns2LongDispRat, HorTrns2LongDispRat, VelConvFact
    REAL (KIND=IKIND)         :: mass_of_parcel,t_spill, Src_Conc, Q_src, tot_mass_rel, TotalMass
    REAL (KIND=IKIND)         :: t_sim, del_t, aquifer_length, aquifer_width, aquifer_height,dx, dy, dz, LPx,  LPy, Pe, Pex, Pey, Xorg, Yorg, distx, disty, BB, XX, YY, XUprRight, YUprRight, XLwrLeft, YLwrLeft, ZLwr, ZUpr
    LOGICAL                                 :: INST, SteadyStateVel,  idz_ReadFile

    CHARACTER (LEN=60)    :: dir, FileVel, FileVelInfo, FileIbnd, FilePumpWell, FileElem, FileNode, FileLayElev
    
       DATA  NDaysMonCum/31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/
        
    TYPE ModelType
         INTEGER, ALLOCATABLE      :: Npar(:,:,:), Ibnd(:,:)
         REAL(KIND=IKIND)                :: Xcrd, Ycrd, Zcrd
         REAL(KIND=IKIND), ALLOCATABLE  :: Con(:,:,:), ConAnal(:,:, :), Uvel(:,:), Vvel(:,:), Wvel(:,:), ElemID(:,:), Elev(:,:,:), AvgLyrTElv(:), VolLyr(:) 
    END TYPE ModelType
    
     TYPE(ModelType) :: Model
     

  !! -------------------------------------------------------------
  !! --- LAGRANGIAN PARCEL TYPE
  !! -------------------------------------------------------------
  !   TYPE TYPELagrParcel
  !      INTEGER                      ::   Igridx, Igridy, Igridz,Elem
  !      REAL(KIND=IKIND) ::   Xcord, Ycord, Zcord 
  !      LOGICAL                       ::   Ingrid          
  !  END TYPE TYPELagrParcel
  !  
  ! TYPE ( TYPELagrParcel), ALLOCATABLE   :: LP (:)
   
  TYPE GridTYPE
           INTEGER                      ::   Ix, Iy, Iz,Elem
           REAL(KIND=IKIND) ::   Xcrd, Ycrd, Zcrd 
   END TYPE GridTYPE
       
   TYPE (GridTYPE), ALLOCATABLE   :: Grid (:)    
       
       

      
   
      TYPE  ELEMENT_TYPE
             INTEGER                :: nnodes,NodeID(max_nodes_elem), ElemID, NeghbrElemID(max_neghbr_elems),NegbrElems, nLP
             REAL(KIND=8)   :: Xcord (max_nodes_elem), Ycord(max_nodes_elem), Zcord(max_nodes_elem), Area, UVel, VVel, WVel
             REAL(KIND=8), ALLOCATABLE   :: AvgLyrTopElv(:)
      CONTAINS
        !PROCEDURE, PASS  :: test_print=> Elem_print
         PROCEDURE, PASS  :: Vel_Interp=> Velocity_Interpolation         ! For quadlilateral and  triangular elements
         PROCEDURE, PASS  :: Vel_Interp_Quad=>  TrilinearInterpolation3D    ! For quadlilateral element
         PROCEDURE, PASS  :: Vel_Interp_Triang=>  TrilinearInterpolation3D_Triang    ! For triangular element
   END  TYPE ELEMENT_TYPE
   
   TYPE (ELEMENT_TYPE), POINTER :: Element (:), Elem
   !TYPE (ELEMENT_TYPE), POINTER :: ThisElement

     TYPE NODE_TYPE
         INTEGER               :: ID, Elem(8), ElemCount
        REAL(KIND=8)   :: Xcord , Ycord, Zcord, Conc,ConcAna, Area, Vel, U1, U2, U3, U4, V1, V2, V3, V4, W1, W2, W3, W4
        REAL(KIND=8), POINTER   ::    Elev(:),Uvel(:), Vvel(:), Wvel(:)
        LOGICAL                :: Used
   END TYPE NODE_TYPE
   
   TYPE (NODE_TYPE), POINTER :: Node (:), NodeC(:)
   
   TYPE RELEASE_TYPE
         INTEGER                                        :: ID, ElemID, ist_day, iend_day, tspill, npar_rel, npar_dt
         REAL (KIND=IKIND)                 :: Xorg, Yorg, Zorg, Conc, FlowRate, MassRate, t_spill, sp_start, sp_end,  tot_mass
   END TYPE RELEASE_TYPE
   
   TYPE(RELEASE_TYPE), ALLOCATABLE  :: Source(:)
   
   TYPE PumpWellType
       INTEGER                                       :: ID, CellID, IR, IC, IL, nParDel
       REAL(KIND=IKIND)                 :: Xcrd, Ycrd, Zcrd
   END TYPE PumpWellType
   
   TYPE (PumpWellType), ALLOCATABLE    :: PumpWell (:) 
   
   
  type quadrilateral
    real(kind=8), dimension(4) :: x, y
    real(kind=8), dimension(4) :: u, v
  end type quadrilateral
   
    CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
    SUBROUTINE Get_Year (jmon, atime, Curnt_STP , jyr)
    ! INTEGER, INTENT (IN) :: jmon
     INTEGER, INTENT (INOUT):: jyr
     INTEGER :: jmon
     LOGICAL, INTENT (INOUT) :: Curnt_STP

     REAL (KIND =IKIND):: atime
     

         IF (atime-365*jyr >NDaysMonCum(jmon) ) THEN
              Curnt_STP=.FALSE.
              jmon=jmon+1
              IF (jmon>12) THEN
                  jmon=jmon-12
                  jyr=jyr+1
              ENDIF
         ENDIF   
    END SUBROUTINE
     
  !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
   SUBROUTINE Print_Time (it, T ) 
    
    INTEGER, INTENT (IN) :: it
    REAL (KIND=IKIND), INTENT(IN) :: T
    
         WRITE(*,'(A10,I5, A2, F6.3, A6 )' ) 'Time step=', it , ":",  T/365, 'years'   
          WRITE(fid_dbug, *) 
          WRITE(fid_dbug,'(A10,I5, A2, F6.3, A6 )' ) 'Time step=', it  , ":",  T/365, 'years'   
          
    END SUBROUTINE Print_Time
    
    
    

  !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
    SUBROUTINE Velocity_Interpolation (ThisElement, xx, yy, zz, jlyr, interpolated_u, interpolated_v,  interpolated_w)
                                                                                         
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
    !--Purpose: Use bilinear/linear interpolation in 2D qudalilateral/triangle horizontal plane and
    !--use linear interpolation in vertical direction
    
      CLASS(ELEMENT_TYPE),  INTENT(IN) :: ThisElement
          INTEGER                                                            ::  nID, jlyr
          REAL (KIND=IKIND), INTENT (IN)          :: xx, yy, zz
          REAL (KIND=IKIND), INTENT (OUT)      :: interpolated_u, interpolated_v, interpolated_w
          REAL (KIND=IKIND)                                     :: intrp_u1, intrp_v1,intrp_w1, intrp_u2, intrp_v2,intrp_w2, zcrd1,zcrd2
   
    !--2D Bilinear/Linear  interpolation of velocity of the top layer
    CALL Velocity_Interpolation_2D (ThisElement, xx, yy, jlyr, intrp_u1, intrp_v1,  intrp_w1)
    
    !--2D Bilinear/Linear interpolation of velocity of the bottom layer
    CALL Velocity_Interpolation_2D (ThisElement, xx, yy, jlyr+1, intrp_u2, intrp_v2,  intrp_w2)
    
   
    !--1D Linear interpolation in the vertical direction
        zcrd1=ThisElement%AvgLyrTopElv(jlyr)
        zcrd2=ThisElement%AvgLyrTopElv(jlyr+1) 
        CALL Velocity_Interpolation_1D ( zz, zcrd1, zcrd2, intrp_u1, intrp_u2, ilyr, interpolated_u)
        CALL Velocity_Interpolation_1D ( zz, zcrd1, zcrd2, intrp_v1, intrp_v2, ilyr, interpolated_v)
        CALL Velocity_Interpolation_1D ( zz, zcrd1, zcrd2, intrp_w1, intrp_w2,  ilyr, interpolated_w)
        
        ! IF (ABS(interpolated_u) >GW_Vel_Limit .OR. ABS(interpolated_v) >GW_Vel_Limit .OR. ABS(interpolated_w) >GW_Vel_Limit ) THEN
        !    write (*,*)  
        !    write (*,*) ' Warning: Interpolated velocity is greater than specified limit'
        !    write (*,*)
        !ENDIF
    

    
   
    END SUBROUTINE Velocity_Interpolation
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     

       
  !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
       SUBROUTINE Velocity_Interpolation_1D ( z, z1, z2, u1, u2, ilyr, interp_u)
    !--Purpose: Use linear  interpolation to caculate velocity at given particle 
    !--used for interpolating only in the vertical direction
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
       INTEGER, INTENT(IN)                                :: ilyr
       REAL(KIND=IKIND), INTENT(IN)          :: z, z1, z2, u1, u2
       REAL(KIND=IKIND), INTENT(OUT)      ::  interp_u
       
       interp_u=u1+(z-z1)/(z2-z1)*(u2-u1)
      

       END SUBROUTINE  Velocity_Interpolation_1D    
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
    
    

    
    
    
 !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
       SUBROUTINE Velocity_Interpolation_2D (ThisElem, x, y, ilyr, interp_u, interp_v, interp_w)
    !--Purpose: Use bilinear or linear interpolation to caculate velocity at given particle 
    !--location. Bilinear interpolation of a qualilateral requires coordinate transformation 
    !--to a unit square
    !--Quadlilateral -Bilinear interpolation ; 
    !--Triangular-Linear Interpolation-using barycentric coordinates.
    !-- Barycentric coordinates are weights that sum up to 1
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
       
        CLASS(ELEMENT_TYPE),  INTENT(IN) :: ThisElem
          INTEGER                                                            ::  nID
          INTEGER                                                            :: max_iter, iter
          INTEGER , INTENT (IN)                                ::  ilyr
          REAL (KIND=IKIND), INTENT (IN)          :: x, y
          REAL (KIND=IKIND), INTENT (OUT)      :: interp_u, interp_v, interp_w

          REAL, PARAMETER :: eps = 1.0e-6 ! To avoid division by zero
          
          REAL (KIND=IKIND) :: x1, x2, x3, x4, y1, y2, y3, y4, aa, bb, cc, dd, dx1, dx2, dx3, dx4, dy1, dy2, dy3, dy4
          REAL (KIND=IKIND) :: u1, u2, u3, u4, v1, v2, v3, v4, denom_x, denom_y, w1,w2,w3,w4
          REAL (KIND=IKIND) ::  lambda1, lambda2, lambda3, alpha, beta
          REAL (KIND=IKIND) ::  detT, N1, N2, N3, N4
          REAL (KIND=IKIND)  ::  xi, eta, tol
          REAL (KIND=IKIND), DIMENSION(4)  :: x_vertices, y_vertices



           
              nID=ThisElem%NodeID(1)
              x1=Node(nID)%Xcord
              y1=Node(nID)%Ycord
              u1=Node(nID)%Uvel(ilyr)          
              v1=Node(nID)%Vvel(ilyr)   
              w1=Node(nID)%Wvel(ilyr)                
             
             nID=ThisElem%NodeID(2)
              x2=Node(nID)%Xcord
              y2=Node(nID)%Ycord
              u2=Node(nID)%Uvel(ilyr)
              v2=Node(nID)%Vvel(ilyr)     
              w2=Node(nID)%Wvel(ilyr)                   
             
             nID=ThisElem%NodeID(3)
              x3=Node(nID)%Xcord
              y3=Node(nID)%Ycord
              u3=Node(nID)%Uvel(ilyr)
              v3=Node(nID)%Vvel(ilyr) 
              w3=Node(nID)%Wvel(ilyr)                   
             
             nID=ThisElem%NodeID(4)
    

              IF (nID >0 ) THEN   ! Quadlilateral-Use bilinear interpolation 
              !-
                          x4=Node(nID)%Xcord
                          y4=Node(nID)%Ycord
                          u4=Node(nID)%Uvel(ilyr)
                          v4=Node(nID)%Vvel(ilyr)  
                         w4=Node(nID)%Wvel(ilyr)                               
                          
                          x_vertices=(/ x1, x2, x3, x4 /)
                          y_vertices=(/y1, y2, y3, y4 /)
                      
                          ! Initial guess for (xi, eta)
                          xi = 0.5d0
                          eta = 0.5d0
                          tol = 1.0d-6
                          max_iter = 100
                      
                        ! Call the Newton-Raphson method
                        !CALL newton_raphson_method (x, y, x_vertices, y_vertices, xi, eta, tol, max_iter, iter)
                          
                        !Perform bilinear interpolation
                         !interp_u=(1-eta)*(1-xi)*u1+xi*(1-eta)*u2+(1-xi)*eta*u3+xi*eta*u4
                         !interp_v=(1-eta)*(1-xi)*v1+xi*(1-eta)*v2+(1-xi)*eta*v3+xi*eta*v4
                         !interp_w=(1-eta)*(1-xi)*w1+xi*(1-eta)*w2+(1-xi)*eta*w3+xi*eta*w4
                        
                       ! Find out the (xi, eta) coordinates for the given (x,y) coordinates by solving two simultaneous equations using Newton Raphson iteration   
                        CALL NewtonRaphson(x, y, x1, y1, x2, y2, x3, y3, x4, y4, xi, eta, max_iter, tol)
                        
                            N1 = 0.25 * (1.0 - xi) * (1.0 - eta)
                            N2 = 0.25 * (1.0 + xi) * (1.0 - eta)
                            N3 = 0.25 * (1.0 + xi) * (1.0 + eta)
                            N4 = 0.25 * (1.0 - xi) * (1.0 + eta)

                             ! Perform bilinear interpolation to compute the velocity at (xi, eta)
                           interp_u = N1 * u1 + N2 * u2 + N3 * u3 + N4 * u4
                           interp_v = N1 * v1 + N2 * v2 + N3 * v3 + N4 * v4
                           interp_w = N1 * w1 + N2 * w2 + N3 * w3 + N4 * w4
                           
                        ! Determine minimum and maximum nodal velocities
                            u_min = MIN(u1, u2, u3, u4)
                            u_max = MAX(u1, u2, u3, u4)
                            v_min = MIN(v1, v2, v3, v4)
                            v_max = MAX(v1, v2, v3, v4)
                            w_min = MIN(w1, w2, w3, w4)
                            w_max = MAX(w1, w2, w3, w4)

              ELSEIF (nID==0) THEN ! Triangle-Use linear interpolation using Barycentric coordinates
                  detT = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
                  lambda1 = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / detT
                  lambda2 = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / detT
                  lambda3 = 1.0 - lambda1 - lambda2
  
                  interp_u = lambda1 * u1 + lambda2 * u2 + lambda3 * u3
                  interp_v = lambda1 * v1 + lambda2 * v2 + lambda3 * v3
                  interp_w = lambda1 * w1 + lambda2 * w2 + lambda3 * w3          
                  
            ! Determine minimum and maximum nodal velocities
                            u_min = MIN(u1, u2, u3)
                            u_max = MAX(u1, u2, u3)
                            v_min = MIN(v1, v2, v3)
                            v_max = MAX(v1, v2, v3)
                            w_min = MIN(w1, w2, w3)
                            w_max = MAX(w1, w2, w3)

              ENDIF
              
                ! Apply slope limiting to prevent overshooting
                IF (interp_u < u_min) THEN
                    interp_u = u_min
                ELSE IF (interp_u > u_max) THEN
                    interp_u = u_min
                END IF
                            
                IF (interp_v < v_min) THEN
                    interp_v = v_min
                ELSE IF (interp_v > v_max) THEN
                    interp_v = v_min
                END IF
                            
                IF (interp_w < w_min) THEN
                    interp_w = w_min
                ELSE IF (interp_w > w_max) THEN
                    interp_w = w_min
                END IF  
              
              
    
       END SUBROUTINE Velocity_Interpolation_2D    
   !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    

 !  !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
 !! Subroutine for Newton-Raphson method
 ! subroutine newton_raphson_method(x, y, x_vertices, y_vertices, xi, eta, tol, max_iter, iter)
 !
 !!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
 !   implicit none
 !   real(8), intent(in) :: x, y, tol
 !   real(8), intent(in), dimension(4) :: x_vertices, y_vertices
 !   integer, intent(in) :: max_iter
 !   real(8), intent(inout) :: xi, eta
 !   integer, intent(out) :: iter
 !   real(8) :: fx, fy, dfdxi_x, dfdxi_y, dfdeta_x, dfdeta_y
 !   real(8) :: jacobian(2, 2), det_jacobian
 !   real(8) :: delta_xi, delta_eta
 !
 !   do iter = 1, max_iter
 !     ! Evaluate the bilinear shape functions and their partial derivatives
 !     call evaluate_functions(x, y, x_vertices, y_vertices, xi, eta, fx, fy, dfdxi_x, dfdxi_y, dfdeta_x, dfdeta_y)
 !
 !     ! Form the Jacobian matrix
 !     jacobian(1, 1) = dfdxi_x
 !     jacobian(1, 2) = dfdeta_x
 !     jacobian(2, 1) = dfdxi_y
 !     jacobian(2, 2) = dfdeta_y
 !
 !     ! Calculate the determinant of the Jacobian
 !     det_jacobian = jacobian(1, 1) * jacobian(2, 2) - jacobian(1, 2) * jacobian(2, 1)
 !
 !     ! Check for singularity
 !     if (abs(det_jacobian) < 1.0d-12) then
 !       print*, 'Jacobian is singular!'
 !       exit
 !     end if
 !
 !     ! Calculate the increments using the inverse of the Jacobian
 !     delta_xi = (jacobian(2, 2) * fx - jacobian(1, 2) * fy) / det_jacobian
 !     delta_eta = (jacobian(1, 1) * fy - jacobian(2, 1) * fx) / det_jacobian
 !
 !     ! Update the guesses
 !     xi = xi - delta_xi
 !     eta = eta - delta_eta
 !
 !     ! Check for convergence
 !     if (abs(delta_xi) < tol .and. abs(delta_eta) < tol) exit
 !   end do
 !
 !   ! If maximum iterations reached without convergence, print a warning
 !   if (iter == max_iter) then
 !     print*, 'Warning: Maximum iterations reached without convergence.'
 !   end if
 !
 ! end subroutine newton_raphson_method
 !!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
 !
 !
 ! 
 !!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
 !  
 !! Subroutine to evaluate the bilinear shape functions and their partial derivatives
 ! subroutine evaluate_functions(x, y, x_vertices, y_vertices, xi, eta, fx, fy, dfdxi_x, dfdxi_y, dfdeta_x, dfdeta_y)
 ! 
 !!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
 ! 
 !   implicit none
 !   real(8), intent(in) :: x, y, xi, eta
 !   real(8), intent(in), dimension(4) :: x_vertices, y_vertices
 !   real(8), intent(out) :: fx, fy, dfdxi_x, dfdxi_y, dfdeta_x, dfdeta_y
 !
 !   real(8) :: N1, N2, N3, N4
 !   real(8) :: dN1dxi, dN2dxi, dN3dxi, dN4dxi
 !   real(8) :: dN1deta, dN2deta, dN3deta, dN4deta
 !
 !   ! Bilinear shape functions
 !   N1 = (1.0d0 - xi) * (1.0d0 - eta)
 !   N2 = xi * (1.0d0 - eta)
 !   N3 = (1.0d0 - xi) * eta
 !   N4 = xi * eta
 !
 !   ! Partial derivatives of the shape functions
 !   dN1dxi = -(1.0d0 - eta)
 !   dN2dxi = (1.0d0 - eta)
 !   dN3dxi = -eta
 !   dN4dxi = eta
 !
 !   dN1deta = -(1.0d0 - xi)
 !   dN2deta = -xi
 !   dN3deta = (1.0d0 - xi)
 !   dN4deta = xi
 !
 !   ! Evaluate the functions
 !   fx = N1 * x_vertices(1) + N2 * x_vertices(2) + N3 * x_vertices(3) + N4 * x_vertices(4) - x
 !   fy = N1 * y_vertices(1) + N2 * y_vertices(2) + N3 * y_vertices(3) + N4 * y_vertices(4) - y
 !
 !   ! Evaluate the derivatives
 !   dfdxi_x = dN1dxi * x_vertices(1) + dN2dxi * x_vertices(2) + dN3dxi * x_vertices(3) + dN4dxi * x_vertices(4)
 !   dfdxi_y = dN1dxi * y_vertices(1) + dN2dxi * y_vertices(2) + dN3dxi * y_vertices(3) + dN4dxi * y_vertices(4)
 !   dfdeta_x = dN1deta * x_vertices(1) + dN2deta * x_vertices(2) + dN3deta * x_vertices(3) + dN4deta * x_vertices(4)
 !   dfdeta_y = dN1deta * y_vertices(1) + dN2deta * y_vertices(2) + dN3deta * y_vertices(3) + dN4deta * y_vertices(4)
 !
 ! end subroutine evaluate_functions  
 !  !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
  
 !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
 
SUBROUTINE NewtonRaphson(x, y, x1, y1, x2, y2, x3, y3, x4, y4, xi, eta, max_iter, tol)
! Newton-Raphson method to find local coordinates (xi, eta) from global coordinates (x, y)
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     

    IMPLICIT NONE

    REAL (8):: x, y, x1, y1, x2, y2, x3, y3, x4, y4
    REAL (8):: xi, eta, dx, dy, detJ, xi_update, eta_update
    REAL (8):: tol
    INTEGER :: iter, max_iter
    REAL (8):: fxi, feta, jac11, jac12, jac21, jac22

    ! Newton-Raphson iteration
    DO iter = 1, max_iter
        ! Compute global coordinates (x_new, y_new) for the current (xi, eta)
        CALL GlobalCoordinates(xi, eta, x1, y1, x2, y2, x3, y3, x4, y4, dx, dy)

        ! Compute residuals (differences between actual (x, y) and computed (dx, dy))
        fxi = x - dx
        feta = y - dy

        ! Check for convergence
        IF (ABS(fxi) < tol .AND. ABS(feta) < tol) THEN
            !PRINT*, "Converged in ", iter, " iterations."
            RETURN
        END IF

        ! Compute the Jacobian matrix for the current (xi, eta)
        CALL ComputeJacobian(xi, eta, x1, y1, x2, y2, x3, y3, x4, y4, jac11, jac12, jac21, jac22)

        ! Compute the determinant of the Jacobian
        detJ = jac11 * jac22 - jac12 * jac21

        ! If determinant is too small, the element may be singular or nearly singular
        IF (ABS(detJ) < tol) THEN
            PRINT*, "Jacobian determinant is too small. Element may be distorted."
            RETURN
        END IF

        ! Newton-Raphson update step for (xi, eta)
        xi_update = (jac22 * fxi - jac12 * feta) / detJ
        eta_update = (-jac21 * fxi + jac11 * feta) / detJ

        ! Update guesses for (xi, eta)
        xi = xi + xi_update
        eta = eta + eta_update
    END DO

    PRINT*, "Warning: Newton-Raphson did not converge within the iteration limit."

END SUBROUTINE NewtonRaphson
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     



!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     

SUBROUTINE GlobalCoordinates(xi, eta, x1, y1, x2, y2, x3, y3, x4, y4, x_new, y_new)
! Compute global coordinates (x_new, y_new) given local coordinates (xi, eta)
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     

    IMPLICIT NONE

    REAL(8) :: xi, eta, x1, y1, x2, y2, x3, y3, x4, y4
    REAL (8):: x_new, y_new
    REAL (8):: N1, N2, N3, N4

    ! Shape functions for bilinear quadrilateral element
    N1 = 0.25 * (1.0 - xi) * (1.0 - eta)
    N2 = 0.25 * (1.0 + xi) * (1.0 - eta)
    N3 = 0.25 * (1.0 + xi) * (1.0 + eta)
    N4 = 0.25 * (1.0 - xi) * (1.0 + eta)

    ! Compute global coordinates using shape functions
    x_new = N1 * x1 + N2 * x2 + N3 * x3 + N4 * x4
    y_new = N1 * y1 + N2 * y2 + N3 * y3 + N4 * y4

END SUBROUTINE GlobalCoordinates
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     


!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
SUBROUTINE ComputeJacobian(xi, eta, x1, y1, x2, y2, x3, y3, x4, y4, jac11, jac12, jac21, jac22)
! Compute the Jacobian matrix for mapping from (xi, eta) to (x, y)
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     

    IMPLICIT NONE

    REAL(8) :: xi, eta, x1, y1, x2, y2, x3, y3, x4, y4
    REAL(8) :: jac11, jac12, jac21, jac22
    REAL(8) :: dN1_dxi, dN2_dxi, dN3_dxi, dN4_dxi
    REAL(8) :: dN1_deta, dN2_deta, dN3_deta, dN4_deta

    ! Derivatives of shape functions with respect to xi and eta
    dN1_dxi  = -0.25 * (1.0 - eta)
    dN2_dxi  =  0.25 * (1.0 - eta)
    dN3_dxi  =  0.25 * (1.0 + eta)
    dN4_dxi  = -0.25 * (1.0 + eta)

    dN1_deta = -0.25 * (1.0 - xi)
    dN2_deta = -0.25 * (1.0 + xi)
    dN3_deta =  0.25 * (1.0 + xi)
    dN4_deta =  0.25 * (1.0 - xi)

    ! Compute Jacobian matrix
    jac11 = dN1_dxi * x1 + dN2_dxi * x2 + dN3_dxi * x3 + dN4_dxi * x4
    jac12 = dN1_deta * x1 + dN2_deta * x2 + dN3_deta * x3 + dN4_deta * x4
    jac21 = dN1_dxi * y1 + dN2_dxi * y2 + dN3_dxi * y3 + dN4_dxi * y4
    jac22 = dN1_deta * y1 + dN2_deta * y2 + dN3_deta * y3 + dN4_deta * y4

END SUBROUTINE ComputeJacobian
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
  
  
    
    
    
    
    
    
    
    


    
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
     SUBROUTINE TrilinearInterpolation3D (ThisElem, xo, yo, zo,ilyr, interpolated_u, interpolated_v, interpolated_w)
    !--Purpose: Use trilinear interpolation to caculate velocity at given particle 
    !--location.
    !--Quadlilateral -Bilinear interpolation ; 
    !--Triangular-Linear Interpolation-using barycentric coordinates.
    !-- Barycentric coordinates are weights that sum up to 1
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
        CLASS(ELEMENT_TYPE),  INTENT(IN)   :: ThisElem
          INTEGER, INTENT (IN)                               ::  ilyr
          INTEGER                                                          ::  nID   
          REAL (KIND=IKIND), INTENT (IN)        ::  xo, yo,zo
          REAL (KIND=IKIND), INTENT (OUT)    ::  interpolated_u, interpolated_v, interpolated_w
          REAL (KIND=IKIND) :: frac_x, frac_y, frac_z, aa, bb, cc, dd
          REAL (KIND=IKIND) :: x(8), y(8), z(8), u(8), v(8), w(8)


          
          REAL (KIND=IKIND) ::  lambda1, lambda2, lambda3, alpha, beta
          REAL (KIND=IKIND) ::  detT

         !--Coordinate positions are set up 
         !- https://en.wikipedia.org/wiki/Trilinear_interpolation
          
        nID=ThisElem%NodeID(1)
            x(1)=Node(nID)%Xcord  
            y(1)=Node(nID)%Ycord
            z(1)=Node(nID)%Elev(ilyr+1)
            u(1)=Node(nID)%Uvel(ilyr+1)
            v(1)=Node(nID)%Vvel(ilyr+1)
            w(1)=Node(nID)%Wvel(ilyr+1)
          
             !-node directly above node 1 (has same x and y coordinates)
            x(5)=Node(nID)%Xcord  
            y(5)=Node(nID)%Ycord             
            z(5)=Node(nID)%Elev(ilyr)
            u(5)=Node(nID)%Uvel(ilyr)
            v(5)=Node(nID)%Vvel(ilyr)
            w(5)=Node(nID)%Wvel(ilyr)
         
         
        nID=ThisElem%NodeID(2)
            x(2)=Node(nID)%Xcord  
            y(2)=Node(nID)%Ycord       
            z(2)=Node(nID)%Elev(ilyr+1)            
            u(2)=Node(nID)%Uvel(ilyr+1)
            v(2)=Node(nID)%Vvel(ilyr+1)
            w(2)=Node(nID)%Wvel(ilyr+1)

              !-node directly above node 1 (has same x and y coordinates)
              x(6)=Node(nID)%Xcord  
              y(6)=Node(nID)%Ycord             
              z(6)=Node(nID)%Elev(ilyr)
              u(6)=Node(nID)%Uvel(ilyr)
              v(6)=Node(nID)%Vvel(ilyr)
              w(6)=Node(nID)%Wvel(ilyr)      
             
       nID=ThisElem%NodeID(3)
           x(3)=Node(nID)%Xcord  
           y(3)=Node(nID)%Ycord       
           z(3)=Node(nID)%Elev(ilyr+1)            
          u(3)=Node(nID)%Uvel(ilyr+1)
          v(3)=Node(nID)%Vvel(ilyr+1)
          w(3)=Node(nID)%Wvel(ilyr+1)

          x(7)=Node(nID)%Xcord  
          y(7)=Node(nID)%Ycord             
          z(7)=Node(nID)%Elev(ilyr)
          u(7)=Node(nID)%Uvel(ilyr)
          v(7)=Node(nID)%Vvel(ilyr)
          w(7)=Node(nID)%Wvel(ilyr)           
          
       nID=ThisElem%NodeID(4)
           x(4)=Node(nID)%Xcord  
           y(4)=Node(nID)%Ycord       
           z(4)=Node(nID)%Elev(ilyr+1)            
          u(4)=Node(nID)%Uvel(ilyr+1)
          v(4)=Node(nID)%Vvel(ilyr+1)
          w(4)=Node(nID)%Wvel(ilyr+1)

          x(8)=Node(nID)%Xcord  
          y(8)=Node(nID)%Ycord             
          z(8)=Node(nID)%Elev(ilyr)
          u(8)=Node(nID)%Uvel(ilyr)
          v(8)=Node(nID)%Vvel(ilyr)
          w(8)=Node(nID)%Wvel(ilyr)            
          

                ! Calculate interpolation factors
                frac_x = (xo - x(1)) / (x(2) - x(1))
                frac_y = (yo - y(1)) / (y(3)- y(1))
                frac_z = (zo - z(1)) / (z(5) - z(1))
                            
                interpolated_u= (1.0 - frac_x) * (1.0 - frac_y)*(1.0 - frac_z) * u(1) +frac_x * (1.0 - frac_y)*(1.0 - frac_z) * u(2)+ &
                                                    (1.0 - frac_x) * frac_y *(1.0 - frac_z)* u(3) + frac_x *frac_y* (1.0 - frac_z) * u(4) +   &
                                                    (1.0 - frac_x) *(1.0 - frac_y) * frac_z * u(5) + frac_x  *(1.0 - frac_y)* frac_z *u(6) + &
                                                    (1.0 - frac_x) * frac_y * frac_z *u(7) + frac_x * frac_y * frac_z *u(8)

                             
                interpolated_v=(1.0 - frac_x) * (1.0 - frac_y)*(1.0 - frac_z) *v(1) +frac_x * (1.0 - frac_y)*(1.0 - frac_z) * v(2)+ &
                                                    (1.0 - frac_x) * frac_y *(1.0 - frac_z)* v(3) + frac_x *frac_y* (1.0 - frac_z) * v(4) +   &
                                                    (1.0 - frac_x) *(1.0 - frac_y) * frac_z * v(5) + frac_x  *(1.0 - frac_y)* frac_z *v(6) + &
                                                    (1.0 - frac_x) * frac_y * frac_z *v(7) + frac_x * frac_y * frac_z *v(8)


                interpolated_w= (1.0 - frac_x) * (1.0 - frac_y)*(1.0 - frac_z) * w(1) +frac_x * (1.0 - frac_y)*(1.0 - frac_z) * w(2)+ &
                                                    (1.0 - frac_x) * frac_y *(1.0 - frac_z)* w(3) + frac_x *frac_y* (1.0 - frac_z) * w(4) +   &
                                                    (1.0 - frac_x) *(1.0 - frac_y) * frac_z * w(5) + frac_x  *(1.0 - frac_y)* frac_z *w(6) + &
                                                    (1.0 - frac_x) * frac_y * frac_z *w(7) + frac_x * frac_y * frac_z *w(8)
                
                
                !--Check if the interpolated velocity is reasonable
                !IF ( abs(interpolated_u) >

              
    
    END SUBROUTINE TrilinearInterpolation3D         
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
   
    
    
    
     
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
     SUBROUTINE TrilinearInterpolation3D_Triang (ThisElem, xo, yo, zo,ilyr, interpolated_u, interpolated_v, interpolated_w)
    !--Purpose: Use trilinear interpolation to caculate velocity at given particle 
    !--location.
    !--Quadlilateral -Bilinear interpolation ; 
    !--Triangular-Linear Interpolation-using barycentric coordinates.
    !-- Barycentric coordinates are weights that sum up to 1
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
        CLASS(ELEMENT_TYPE),  INTENT(IN)   :: ThisElem
          INTEGER, INTENT (IN)                               ::  ilyr
          INTEGER                                                          ::  nID   
          REAL (KIND=IKIND), INTENT (IN)        ::  xo, yo,zo
          REAL (KIND=IKIND), INTENT (OUT)    ::  interpolated_u, interpolated_v, interpolated_w
          REAL (KIND=IKIND) :: frac_x, frac_y, frac_z, aa, bb, cc, dd
          REAL (KIND=IKIND) :: x(8), y(8), z(8), u(8), v(8), w(8)



          
          REAL (KIND=IKIND) ::  lambda1, lambda2, lambda3, alpha, beta
          REAL (KIND=IKIND) ::  detT

         !--Coordinate positions are set up 
         !- https://en.wikipedia.org/wiki/Trilinear_interpolation
          
           nID=ThisElem%NodeID(1)
                x(1)=Node(nID)%Xcord  
                y(1)=Node(nID)%Ycord
                z(1)=Node(nID)%Elev(ilyr+1)
                u(1)=Node(nID)%Uvel(ilyr+1)
                v(1)=Node(nID)%Vvel(ilyr+1)
                w(1)=Node(nID)%Wvel(ilyr+1)
       
                 !-node directly above node 1 (has same x and y coordinates)
                   x(4)=Node(nID)%Xcord  
                   y(4)=Node(nID)%Ycord       
                   z(4)=Node(nID)%Elev(ilyr)            
                  u(4)=Node(nID)%Uvel(ilyr)
                  v(4)=Node(nID)%Vvel(ilyr)
                  w(4)=Node(nID)%Wvel(ilyr)
         
         
        nID=ThisElem%NodeID(2)
                x(2)=Node(nID)%Xcord  
                y(2)=Node(nID)%Ycord       
                z(2)=Node(nID)%Elev(ilyr+1)            
                u(2)=Node(nID)%Uvel(ilyr+1)
                v(2)=Node(nID)%Vvel(ilyr+1)
                w(2)=Node(nID)%Wvel(ilyr+1)


                  !-node directly above node 1 (has same x and y coordinates)
                 x(5)=Node(nID)%Xcord  
                y(5)=Node(nID)%Ycord             
                z(5)=Node(nID)%Elev(ilyr)
                u(5)=Node(nID)%Uvel(ilyr)
                v(5)=Node(nID)%Vvel(ilyr)
                w(5)=Node(nID)%Wvel(ilyr)
             
       nID=ThisElem%NodeID(3)
               x(3)=Node(nID)%Xcord  
               y(3)=Node(nID)%Ycord       
               z(3)=Node(nID)%Elev(ilyr+1)            
              u(3)=Node(nID)%Uvel(ilyr+1)
              v(3)=Node(nID)%Vvel(ilyr+1)
              w(3)=Node(nID)%Wvel(ilyr+1)

              x(6)=Node(nID)%Xcord  
              y(6)=Node(nID)%Ycord             
              z(6)=Node(nID)%Elev(ilyr)
              u(6)=Node(nID)%Uvel(ilyr)
              v(6)=Node(nID)%Vvel(ilyr)
              w(6)=Node(nID)%Wvel(ilyr)      

      
                ! Calculate interpolation factors
                frac_x = (xo - x(1)) / (x(2) - x(1))
                frac_y = (yo - y(1)) / (y(3)- y(1))
                frac_z = (zo - z(1)) / (z(4) - z(1)) 
                            

                 interpolated_u= (1.0 - frac_x) * (1.0 - frac_y) * u(1) +frac_x * (1.0 - frac_y) * u(2) + &
                                                    (1.0 - frac_x) * frac_y * u(3) + frac_x * (1.0 - frac_z) * u(4) +   &
                                                    (1.0 - frac_x) * frac_z * u(5) + frac_x * frac_z * u(6) 
                             
                interpolated_v= (1.0 - frac_x) * (1.0 - frac_y) *v(1) +frac_x * (1.0 - frac_y) * v(2) + &
                                                    (1.0 - frac_x) * frac_y * v(3) + frac_x * (1.0 - frac_z) * v(4) +   &
                                                    (1.0 - frac_x) * frac_z * v(5) + frac_x * frac_z *v(6) 

                interpolated_w= (1.0 - frac_x) * (1.0 - frac_y) * w(1) +frac_x * (1.0 - frac_y) * w(2)+ &
                                                    (1.0 - frac_x) * frac_y * w(3) + frac_x * (1.0 - frac_z) * w(4) +   &
                                                    (1.0 - frac_x) * frac_z * w(5) + frac_x * frac_z *w(6)               
                
    
    END SUBROUTINE TrilinearInterpolation3D_Triang         
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
     
!!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
!
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
    SUBROUTINE Get_Elem_LP (ielem, xcrd, ycrd) 
    !--Purpose: Figure out the element number that the particle mp
    !-- is located. Travese through all neighbouring elements (starting 
    !-- from the previous element. Neighbouring elements should be calculated first
    !--In this function, the winding number algorithm is applied by checking the cross products 
    !--of the edges of the quadrilateral and the vectors from the point to the vertices of the quadrilateral. 
    !--If all cross products have the same sign (positive or negative), the point is inside the quadrilateral. Otherwise, it's outside.
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
    INTEGER, INTENT(INOUT)                :: ielem
    REAL (KIND=IKIND), INTENT(IN)  :: xcrd, ycrd
    TYPE (ELEMENT_TYPE), POINTER :: ThisElem
    
    INTEGER                                                  :: nID, jelem
    REAL (KIND=IKIND)                           :: xcrd1, xcrd2, xcrd3, xcrd4, ycrd1, ycrd2, ycrd3, ycrd4, area_ABC, area_PAB, area_PBC, area_PAC
    REAL (KIND=IKIND)                           :: detT, alpha, beta, gamma, cross_product
    LOGICAL                                                   :: is_inside
    
    !ielem= LP(mp)%Elem
    !xcrd=LP(mp)%xcord
    !ycrd=LP(mp)%ycord
    
    ThisElem=>Element(ielem)
    
      !--First check if the particle is still in the current element
             CALL Get_nodal_coord(ielem, 1, xcrd1, ycrd1)
             CALL Get_nodal_coord(ielem, 2, xcrd2, ycrd2)
             CALL Get_nodal_coord(ielem,3, xcrd3, ycrd3)

              
          IF ( ThisElem%NodeID(4)>0) THEN  
                           CALL Get_nodal_coord(ielem, 4, xcrd4, ycrd4)           
                          is_inside=Is_Inside_Quad(xcrd, ycrd, xcrd1, xcrd2, xcrd3, xcrd4, ycrd1, ycrd2, ycrd3, ycrd4)
                          
                          IF (is_inside ==  .false.) THEN
                              jelem= Get_neighbr_elem (ielem ,xcrd, ycrd)    
                              !--Update the new element of the LP
                              ielem=jelem
                          ENDIF
                          
          ELSE     ! Triangular mesh
                        is_inside= Is_Inside_Triangle ( xcrd, ycrd, xcrd1, xcrd2, xcrd3, ycrd1, ycrd2, ycrd3   )
                        IF (is_inside ==  .false.) THEN
                                jelem= Get_neighbr_elem (ielem ,xcrd, ycrd)    
                                  !--Update the new element of the LP
                                 ielem=jelem     
                        ENDIF  
          END IF
          
    END SUBROUTINE Get_Elem_LP
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
    FUNCTION Get_neighbr_elem (ielem,xcord, ycord) RESULT (jelem)
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
       INTEGER, INTENT(IN)                          :: ielem
       REAL (KIND=IKIND), INTENT(IN)   :: xcord, ycord
       INTEGER                      :: jelem
        TYPE (ELEMENT_TYPE), POINTER :: ThisElem, ThisEl
        INTEGER                                                  :: ielm, nID, kelem
        REAL (KIND=IKIND)                           ::  xcrd1, xcrd2, xcrd3, xcrd4, ycrd1, ycrd2, ycrd3, ycrd4
        REAL (KIND=IKIND)                           :: area_ABC, area_PAB, area_PBC, area_PAC
        REAL (KIND=IKIND)                           :: cross_product
        LOGICAL                                                   :: is_inside
    
            is_inside=.FALSE.
           ThisElem=>Element(ielem)
           I=0
          DO WHILE ( NOT(is_inside)  .AND.  I < ThisElem%NegbrElems )
                 I=I+1
                 kelem= ThisElem%NeghbrElemID(I)
                 ThisEl=>Element(kelem)
                 ielm=ThisEl%ElemID
                 CALL Get_nodal_coord(ielm, 1, xcrd1, ycrd1)
                 CALL Get_nodal_coord(ielm, 2, xcrd2, ycrd2)
                 CALL Get_nodal_coord(ielm,3, xcrd3, ycrd3)
                 
                 IF  (ThisEl%NodeID(4)>0 ) THEN
                            CALL Get_nodal_coord(ielm, 4, xcrd4, ycrd4)        
                            is_inside=Is_Inside_Quad(xcord, ycord, xcrd1, xcrd2, xcrd3, xcrd4, ycrd1, ycrd2, ycrd3, ycrd4)
                 ELSE 
                             is_inside= Is_Inside_Triangle ( xcord, ycord, xcrd1, xcrd2, xcrd3, ycrd1, ycrd2, ycrd3   )
                 ENDIF
          END DO
        IF (is_inside ) THEN
               jelem=kelem
        ELSE IF (I==max_neghbr_elems) THEN
               WRITE (*,*) " LP has moved out of neighboring elements"
               WRITE (*,*) "Program stops"
               PAUSE
        ENDIF 
    END FUNCTION  Get_neighbr_elem
    
    
    FUNCTION GetElemArea (ielem) RESULT (area)
    
    
      INTEGER, INTENT(IN)                          :: ielem

      TYPE (ELEMENT_TYPE), POINTER      :: ThisElem
     INTEGER                                                  :: ielm, nID, kelem
    REAL (KIND=IKIND)  :: x1,x2,x3,x4,y1,y2,y3,y4
        
    ThisElem=>Element(ielem)
    
    x1=ThisElem%xcord(1); x2=ThisElem%xcord(2); x3=ThisElem%xcord(3); x4=ThisElem%xcord(4)
    y1=ThisElem%ycord(1); y2=ThisElem%ycord(2); y3=ThisElem%ycord(3); y4=ThisElem%ycord(4)
    
    IF (ThisElem%NodeID(4)==0)    THEN ! Triangle
        area = 0.5 * ABS( x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2) )
     ELSE 
         ! Use Shoelace formula to calculate the area
         area = 0.5 * ABS( (x1*y2 + x2*y3 + x3*y4 + x4*y1) - (y1*x2 + y2*x3 + y3*x4 + y4*x1) )
     END IF
     
        
        END FUNCTION GetElemArea
    
    SUBROUTINE Get_nodal_coord (jelem, nID, x, y)
           INTEGER, INTENT(IN)                              :: jelem, nID
           REAL (KIND=IKIND), INTENT(OUT)   :: x, y
          TYPE (ELEMENT_TYPE), POINTER        :: ThisElem
          INTEGER                                                         :: inode

          
                  ThisElem=>Element(jelem)
                  inode=ThisElem%NodeID(nID)
                  x=Node(inode)%Xcord
                  y=Node(inode)%Ycord


    END SUBROUTINE Get_nodal_coord
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
    
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
    FUNCTION Is_Inside_Quad ( x, y, x1, x2, x3, x4, y1, y2, y3, y4   )  RESULT(is_inside)
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
    REAL (KIND=IKIND) , INTENT (IN):: x, y, x1, x2, x3, x4, y1, y2, y3, y4
    LOGICAL                                                   :: is_inside
    REAL (KIND=IKIND)                           :: cross_product
    
             !-- Usingwinding number algorithm        
                      cross_product = (x2 - x1) * (y - y1) - (x - x1) * (y2 - y1)
                        IF (cross_product > 0) THEN
                            cross_product = (x3 - x2) * (y - y2) - (x - x2) * (y3 - y2)
                            IF  (cross_product > 0) THEN
                                cross_product = (x4 - x3) * (y - y3) - (x - x3) * (y4 - y3)
                                IF (cross_product > 0) THEN
                                    cross_product = (x1 - x4) * (y - y4) - (x - x4) * (y1 - y4)
                                   IF (cross_product > 0) THEN
                                        is_inside = .true.
                                    ELSE
                                        is_inside = .false.
                                    END IF
                                ELSE
                                    is_inside = .false.
                                 END IF
                            ELSE
                                is_inside = .false.
                             END IF
                          ELSE
                            is_inside = .false.
                          END IF
    
                          !RETURN  is_inside
    
    END FUNCTION Is_Inside_Quad
 !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!   
     FUNCTION Is_Inside_Triangle ( x, y, x1, x2, x3, y1, y2, y3   )  RESULT(is_inside)
            !!-- A(x1,y1), B(x2,y2), C(x3,y3) are the 3 nodes and P is the location of the Lagrangian particle
            !!--If P is inside the triangle Area ABC= Area PAB+Area PBC+ Area PAC

           !--Use barycentric coordinate method
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
    REAL (KIND=IKIND) , INTENT (IN):: x, y, x1, x2, x3, y1, y2, y3
    
    REAL (KIND=IKIND) ::denom, lambda1, lambda2,lambda3, area_ABC, area_PAB, area_PAC
    LOGICAL                                                   :: is_inside
            
                !!--Area of the triangle ABC
                 area_ABC=abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2 
                 
                 area_PAB=abs(x*(y2-y1)+x1*(y-y2)+x2*(y1-y))/2   
                 area_PBC=abs(x*(y3-y2)+x2*(y-y3)+x3*(y2-y))/2    
                 area_PAC=abs(x*(y1-y3)+x3*(y-y1)+x1*(y3-y))/2    
                
                
                IF (ABS(area_ABC-(area_PAB+area_PBC+area_PAC))<0.1  ) THEN
                     is_inside = .true.
                ELSE 
                      is_inside = .false.                  
                END IF 
                
                !
                 ! Calculate the denominator
                    !denom = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
                    !
                    !! Calculate barycentric coordinates
                    !lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denom
                    !lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denom
                    !lambda3 = 1.0 - lambda1 - lambda2
                    !
                    !! Check if point is inside the triangle
                    !IF (lambda1 >= 0.0 .and. lambda2 >= 0.0 .and. lambda3 >= 0.0) THEN
                    !     is_inside = .true.
                    !ELSE
                    !    is_inside = .false. 
                    !ENDIF
                    !
                
                

     END FUNCTION Is_Inside_Triangle
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
     
     
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
    FUNCTION Get_Lyr(ielem, elev)  RESULT (Layer) 
    !--Purpose: Calculate the average elevations of all layers of a given element 
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
    INTEGER, INTENT (IN)                         :: ielem
    REAL(KIND=IKIND), INTENT (IN)   :: elev
    
    INTEGER :: nnds
    !REAL (KIND=IKIND) :: avg_lyr_elev(num_Lyrs+1), max_elev
    
    Elem=>Element(ielem)
    nnds=Elem%nnodes

    Layer=1
    DO WHILE (Elem%AvgLyrTopElv(Layer)>elev) 
        Layer=Layer+1
    ENDDO
    Layer=Layer-1

   if (Layer==0 ) then
       write(*,*) 'Release elevation is outside the range'
       WRITE(fid_dbug, *)  'Release elevation is outside the range'
       PAUSE
    endif
    
    END FUNCTION Get_Lyr
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!     
   
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!      
    SUBROUTINE Check_bounds ( ielem, xcord, ycord, zcord)
    !--Purspose: Check if the LP has moved out of domain and if it's moved out of domain then place it back at the origin 
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!      
    INTEGER, INTENT(IN)    ::ielem
     REAL (KIND=IKIND), INTENT(INOUT)  ::xcord, ycord, zcord
     
    REAL (KIND=IKIND)  :: z_offset,  x_offset, y_offset
    
    
                       
            !  --Check if the parcel has crossed the top of the layer 1 and Reset the z coordinate
             IF (  zcord> Element(ielem)%AvgLyrTopElv(1) ) THEN  
                    z_offset=zcord- Element(ielem)%AvgLyrTopElv(1)
                    zcord= Element(ielem)%AvgLyrTopElv(1)-z_offset
             ENDIF
            !--Check if the parcel has crossed the bottom of the bottom layer and Reset the z coordinate
             IF (zcord< Element(ielem)%AvgLyrTopElv(num_Lyrs) ) THEN  
                    z_offset= Element(ielem)%AvgLyrTopElv(num_Lyrs) -zcord
                    zcord= Element(ielem)%AvgLyrTopElv(num_Lyrs)+z_offset
             ENDIF   
             !--Check if the parcel has crossed the western boundary and Reset X coordiante
             IF (xcord<XLwrLeft  ) THEN
                 x_offset=XLwrLeft-xcord
                 xcord=XLwrLeft+x_offset
             END IF
               !--Check if the parcel has crossed the eastern boundary and Reset X coordiante
             IF (xcord>XUprRight ) THEN                 
                 x_offset=xcord-XUprRight
                 xcord=XUprRight-x_offset
             END IF
            !--Check if the parcel has crossed the southern boundary and Reset Y coordiante
             IF ( ycord<YLwrLeft  ) THEN
                 y_offset=YLwrLeft-ycord
                 ycord=YLwrLeft+y_offset
             END IF
             !--Check if the parcel has crossed the nothern boundary and Reset Y coordiante
              IF ( ycord >YUprRight ) THEN
                 y_offset=ycord-YUprRight 
                ycord=YUprRight-y_offset
              END IF


            !              
            !  !--Check if the parcel has crossed the top of the layer 1 and Reset the z coordinate
            ! IF (  LP(nlp)%zcord> Element(ielem)%AvgLyrTopElv(1) ) THEN  
            !        z_offset=fact_bnd_offset*(Element(ielem)%AvgLyrTopElv(1)-Element(ielem)%AvgLyrTopElv(2))         ! 0.1 is arbiteratry should be tested
            !        LP(nlp)%zcord= Element(ielem)%AvgLyrTopElv(1)-z_offset
            ! ENDIF
            !!--Check if the parcel has crossed the bottom of the bottom layer and Reset the z coordinate
            ! IF (  LP(nlp)%zcord< Element(ielem)%AvgLyrTopElv(num_Lyrs) ) THEN  
            !         z_offset=fact_bnd_offset*(Element(ielem)%AvgLyrTopElv(num_Lyrs-1)-Element(ielem)%AvgLyrTopElv(num_Lyrs))         ! 0.1 is arbiteratry should be tested
            !         LP(nlp)%zcord= Element(ielem)%AvgLyrTopElv(num_Lyrs)+z_offset
            ! ENDIF   
            ! !--Check if the parcel has crossed the western boundary and Reset X coordiante
            ! IF ( LP(nlp)%xcord-XLwrLeft<0.00  ) THEN
            !     x_offset=fact_bnd_offset*dx
            !     LP(nlp)%xcord=XLwrLeft+x_offset
            ! END IF
            !   !--Check if the parcel has crossed the eastern boundary and Reset X coordiante
            ! IF (  (XUprRight-LP(nlp)%xcord)<0.00 ) THEN                 
            !     x_offset=fact_bnd_offset*dx
            !     LP(nlp)%xcord=XUprRight-x_offset
            ! END IF
            !!--Check if the parcel has crossed the southern boundary and Reset Y coordiante
            ! IF ( LP(nlp)%ycord-YLwrLeft <0.00 ) THEN
            !     y_offset=fact_bnd_offset*dy
            !     LP(nlp)%ycord=YLwrLeft+y_offset
            ! END IF
            ! !--Check if the parcel has crossed the nothern boundary and Reset Y coordiante
            !  IF ( YUprRight-LP(nlp)%ycord < 0.00 ) THEN
            !     y_offset=fact_bnd_offset*dy
            !     LP(nlp)%ycord=YUprRight-y_offset
            ! END IF
                       
    END SUBROUTINE Check_bounds          
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!  
    
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
    SUBROUTINE InitCond ( ) 
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
        !--INITIAIL CONDITIONS
	    !Model%Con(1,1)=MassRate		              
	    !Model%ConAnal(1,1)=MassRate
        Model%Npar(:,:,:)=0
         num_par_del_PW1=0
	    !npart=2000
	    !npart_prev=1
   END SUBROUTINE InitCond
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
!    
!       SUBROUTINE Make_List(npar)
!!
!!       create a list for multi-species problem, and the link(0)
!!       correspond to the top of the avail list
!!
!!       beg(i), i=1: sus par
!!                     i=2: sus_lnd par
!!             
!       
!!	  USE AdvDispLagrPar_Mod
!      
!      
!        DO  10 i=0,npar
!                 Link(i)=mod(i+1, npar+1)
!10      CONTINUE 
!        
!       DO 20 i=1,n_partype 
!                beg(i)=-1
!20     CONTINUE
!       
!       RETURN
!    END SUBROUTINE   Make_List 
!    
!    
!  SUBROUTINE List_Insert(ista, new)
!
!       !insert a new particle to the top of ista-species,
!       !the subscript of the new particle is equal to "nnew".
!
!        !USE AdvDispLagrPar_Mod
!        INTEGER avail, loc
!
!        avail=Link(0)      ! the top of avail list
!        loc=beg(ista)
!        new=avail
!        avail=Link(avail)
!       IF (avail.eq.0) THEN
!           WRITE (*,*) 'Particle numbers overpass the maximum when the', ista, 'th kind of particles are inserted.'
!           STOP
!        ENDIF
!        Link(new)=loc
!        beg(ista)=new
!        Link(0)=avail
!        
!        RETURN
!        
!        END SUBROUTINE List_Insert
!!
!    
    

    
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!    
        SUBROUTINE Gauss(  s, am, v)
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!
        INTEGER       :: iseedx, iseedy
        REAL           :: a, s, am, y
        REAL           :: v
!
!       s=standard devation, am=mean
!
        a=0.0
        DO 10  i=1,12
           ! CALL randu(iseedx ,iseedy, y)
            CALL RANDOM_NUMBER(y)
 10       a=a+y
             v=(a-6.0)*s+am
        RETURN
        
        END SUBROUTINE Gauss
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------!        
        
        
!--USGS Solution----------------------------------------------------------------------------------------------------------------
        
        SUBROUTINE CNRML2(QM,POR,DK,T, X, Y, DX, DY, VX,CN,NMX) 
        IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
        COMMON /IOUNIT/ IN, IO 
        COMMON /GLPTS/ WN(256), ZN(256) 
        !THIS ROUTINE CALCULATES SOLUTE CONCENTRATION AT X,Y BASED ON 
        !THE ANALYTIC SOLUTION TO THE TWO-DIMENSIONAL ADVECTIVE160 
        !TECHNIQUESOFWATER-RESOURCESINVESTIGATIONS 

        !DISPERSIVE SOLUTE TRANSPORT EQUATION FOR AN AQUIFER OF 
        !INFINITE AREAL EXTENT WITH A CONTINUOUS POINT SOURCE LOCATED 
        !AT X-XC AND Y-YC. THE INTEGRAL FROM 0 TO T IS EVALUATED 
        !USING A GAUSS-LEGENDRE QUADRATURE INTEGRATION TECHNIQUE. 
                !PI=3.14159265358979D0 
                CN=0.0D0 
                !FOR T-O, ALL CONCENTRATIONS EQUAL 0.0 
                IF(T .LE. 0.0D0) RETURN 
                !START NUMERICAL INTEGRATION LOOP 
                ALPHA=X*X/(4.0D0*DX)+Y*Y/(4.0D0*DY) 
                BETA=VX*VX/(4.0D0*DX)+DK 
                VX2D=VX*X/(2.0D0*DX) 
                SUM=0.00D0 
                DO 20 I=1,NMX 
                !SCALE THE GAUSS-LEGENDRE COEFFICIENTS TO ACCOUNT FOR THE 
                !NON-NORMALIZED LIMITS OF INTEGRATION 
                WI=WN(I) 
                ZI=T*(ZN(I)+1.0D0)/2.0D0 
                !TERM1 
                X1=-ALPHA/ZI-BETA*ZI 
                X1=DEXP(X1)/ZI 
                SUM=SUM+ X1*WI 
        20        CONTINUE 
                SUM=SUM*T/2.0D0 
                CN=QM*SUM*DEXP(VX2D)/(4.0D0*POR*PI*DSQRT(DX*DY)) 
        RETURN 
        END

    

        !******************************************************** 
        !* * 
        !* SUBROUTINE GLQPTS * 
        !* * 
        !* VERSION CURRENT AS OF 10/01/87 * 
        !* * 
        !* * 
        !*****************************************************~~** 
        SUBROUTINE GLQPTS (N) 
        IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
        CHARACTER*1 SKIP 
        COMMON /GLPTS/ WN(256),ZN(256) 
        COMMON /IOUNIT/ IN, IO
        !THIS ROUTINE READS THE NORMALIZED ROOTS ZN(1) AND WEIGHTS WN(1) 
        !OF THE LEGENDRE POLYNOMIALS FROM THE DATA FILE 'GLQ.PTS' 
        !N IS THE NUMBER OF INTEGRATION POINTS AND CAN ONLY HAVE A 
        !VALUE OF EITHER 4,20,60,104,0R 256 
        IN2=77 
        OPEN(IN2,FILE='Input\GLQ.PTS',STATUS='OLD') 
        !SKIP LINES IN FILE UNTIL CORRECT COEFFICIENTS ARE REACHED 
        ISKIP=-1
        !ISKIP=-1 
        IF(N .EQ. 4) ISKIP= 7
        IF(N .EQ. 20) ISKIP=9 
        IF(N .EQ. 60) ISKIP= 15
        IF(N .EQ. 104) ISKIP= 31
        IF(N .EQ. 256) ISKIP=57
        IF (ISKIP .EQ. -1) WRITE(IO ,201) 
        IF (ISKIP .EQ. -1) STOP 
        DO 60 I=1,ISKIP 
        60 READ(IN2,101) SKIP 
   
        !READ IN ZN(1) AND WN(I), FOUR VALUES PER LINE 
        NC=N/8 
        IF (MOD(N,8) .NE. 0) NC=NC+1
        DO 80 I=1,NC 
            K=(I-1)*8-1
            READ(IN2,102) (ZN(K+J*2),J=1,4) 
        80    CONTINUE 
        DO 100 I=1,NC 
            K=(I-1)*8-1 
            READ(IN2,102) (WN(K+J*2),J=1,4) 
        100  CONTINUE 
        !FILL IN THE SYMMETRIC TERMS 
        DO 120 J=2,N,2 
                    J1=J-1 
                    ZN(J)=-ZN(J1) 
        120  WN(J)=WN(J1) 
        CLOSE(IN2) 
        RETURN 
        ! 56 
        !1 FORMAT STATEMENTS 57 
        101 FORMAT(A1) 
        102 FORMAT(4D20.0) 
        201 FORMAT(1H0,20X,'***** ERROR IN ROUTINE GWPTS *****'/  &
         1H ,20X,'NO. OF ROOTS SPECIFIED MUST EQUAL 4,20,60,104 OR 256')  
        END
!    

        
!--Wilson and Miller (1978)         
        !--For GAUSS integration of the analytical solution
        FUNCTION GAUSS_INTGRTN (A, B, FUNCTN)
        
        !NUMERICAL INTEGRATION BY 24 POINT GAUSS-LEGENDRE QUAORATURE
        !ZEROS ANO WEIGHTING FACTORS ARE FROM TABLE 25.4. P916. OF
        !ABRAMOWITZ ANO STENGUN(1966)
        
        INTEGER   J
        REAL(KIND=8)  A, B, C, D, FUNCTN, GAUSS,  SUM, W, Z
        DIMENSION Z(12), W(12)
        
        DATA Z/0.064056892862065, 0.191118867473616,  0.315042679696163, &
                     0.433793507626045, 0.545421471388839,  0.648093651936975,  &
                     0.740124191578554, 0.820001985973902,  0.886415527004401,  &
                     0.938274552002732, 0.974728555971309,  0.995187219997021/
        DATA W/0.127938195346752, 0.125837456346828, 0.121670472927803,  &
                       0.115505668053725, 0.107444270115965, 0.097618652104113,  &
                       0.066190161531953, 0.073346481411080, 0.059296584915436, &
                        0.044277438817419, 0.028531388628933, 0.012341229799987/
        
          !--SET UP INITIAL PARAMETERS
            C=(B-A)/2.0D00 
            D= (B+A)/2.0D00
            
          !--ACCUMULATE THE SUM IN THE 24-POINT FORMULA
            SUM= 0.0
            DO  J=1,12
                IF (Z(J) ==0.00) THEN
                    SUM=SUM+W(J)*FUNCTN(D)
                ELSE
                    SUM=SUM+W(J)*( FUNCTN(Z(J)*C+D)+ FUNCTN(-Z(J)*C+D) )
                ENDIF           
            ENDDO
            
         !--MAKE INTERVAL CORRECTION AND RETURN
            GAUSS_INTGRTN=C*SUM
          
        RETURN
        END FUNCTION GAUSS_INTGRTN
        
        FUNCTION FUNCTN (Z) 
        ! Integrand of Hantush Well function
            REAL (KIND=8)   DB, Z, ARG, FUNCTN
           !COMMON BF
           DB=BB
           ARG=DLOG(Z)+Z+DB*DB/(4.0D00*Z)
           FUNCTN = DEXP(-ARG)
           
           RETURN
        END FUNCTION FUNCTN
        
        FUNCTION BKO ( Z)
        ! EVALUATION OF MODIFIEO BESSEL FUNTION OF SECOND KIND
        ! OF ORDER ZERO
        ! POLYNOMIAL APPROXIMATIONS ARE USED FOR KO(Z)
        ! SEE SECTION 9.8 OF ABRAMOWITZ ANO STEGUN (1966)       
        !        
        
        REAL (KIND=8)  Z, T, T2, T4, T6, T8, T10, T12, SUM
        
        IF(Z <= 0.0 ) GO TO 200
        T=Z/2.0
        T2=T*T
        T4=T2*T2
        T6=T2*T4
        IF (Z > 2.0)  GOTO 100
        T8=T2*T6
        T10=T2*T8
        T12=T2*T10
        BKO=-1.0*LOG(T)*BIO(Z)-0.57721566+0.42278420*T2+ 0.23069756*T4 + 0.03488590*T6  &
                   + 0.00262698*T8 + 0.00010750*T10+ 0.0740E-04*T12
        RETURN
100     CONTINUE
        SUM=(1.25331414 - 0.07832358/T + 0.02189568/T2- 0.01062446/(T*T2) + 0.00587872/T4  &
                    - 0.00251540/(T*T4) + 0.00053208/T6)
        BK0LOG=LOG(SUM)-Z-0.5*LOG(Z)
        BKO=EXP(BK0LOG)
        RETURN
200     CONTINUE
        WRITE(NO,205) Z
  205 FORMAT(6X,'ARGUMENT OF 8ESSEL FUNCTION KO(Z) IS LESS THAN',  &
                  ' OR EQUAL TO ZERO' ,/,6X,'Z = ',E12.6,' -- PROGRAM TERMINATED')       
        END FUNCTION BKO
        

        
       FUNCTION BIO (Z)        
    ! EVALUATION OF MOOIFIEO 8ESSEL FUNCTION OF THE FIRST KINO
    ! OF OROER ZERO
    ! POLYNOMIAL APPROXIMATIONS ARE USEO FOR 10(Z)
    ! SEE SECTION 9.6 OF ABRAMOWITZ ANO STEGUN (1966)
       DIMENSION A(9) 
       DATA A/0.9169365, 4.32105045, 6.09540629, 6.45306739, 4.6926023,  &
                    13.66357842, 3.63608323, 4.10583047, 5.5407023531/
        
        INTEGER I
        REAL (KIND=8)  Z, T, T2, T4, T6, T8, T10, T12, SUM, SIGN, ARG, BIOLOG
        
        IF(Z <= 0.0 ) GO TO 200
        T=Z/3.75
        IF (Z > 3.75) GO TO 100
        T2=T*T
        T4=T2*T2
        T6=T2*T4
        T8=T2*T6
        T10=T2*T8
        T12=T2*T10
        BIO=1.0++ 3.5156229*T2 + 3.0899424*T4 + 1.2067492*T6  &
                 + 0.2659732*T8 + 0.0360768*T10+ 0.0045813*T12
        RETURN
100     CONTINUE 
        SUM=0.00
        DO 150  I=1, 9
            SIGN=(-1.0)**(I+1)
            IF (I==2) SIGN=1.0
            ARG=-1.0*A(I)-0.5*LOG(Z)-FLOAT(I-1)*LOG(T)
            SUM=SUM+SIGN*EXP(ARG)
150     CONTINUE 
        BIOLOG= LOG(SUM)+Z
        BIO=EXP(BIOLOG)
        RETURN
200     CONTINUE
        IF (Z <0) GOTO 300
        BIO =1.0
        RETURN
300     WRITE (*,305) Z
305  FORMAT(6X, 'ARGUMENT OF BESSEL FUNCTION 10(Z) IS NEGATIVE' ,/,    &
          16X, 'Z * . ,E12.6,· -- PROGRAM TERMINATED')
        END FUNCTION BIO
     
     
     END MODULE AdvDispLagrPar_Mod
    