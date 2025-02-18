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
     
    MODULE LinkedList_Mod
    IMPLICIT NONE
            !!--Define a node to reprent a Lagrangian Particle
            TYPE:: LPNodeType
                INTEGER                          :: ID, Igridx, Igridy, Igridz, Elem
                REAL (KIND=8)    :: xcord, ycord, zcord
                LOGICAL                          :: Ingrid
                TYPE(LPNodeType), POINTER  :: Next=>NULL()
            END TYPE LPNodeType

            !--Define link list for LPs
            TYPE  LPLinkeListType
                TYPE (LPNodeType), POINTER :: Head=>NULL()
                TYPE (LPNodeType), POINTER :: Tail=>NULL()        ! Pointer to last node
                INTEGER :: TotNumParcels, NumParDelete
            CONTAINS
                PROCEDURE  :: Initialize =>InitalizeList
                PROCEDURE ::  InsertLP => InsertNode
               PROCEDURE ::  RemoveLP=> RemoveNode
               ! PROCEDURE ::  PrintList=>PrintLPList
            
            END TYPE LPLinkeListType
        
        
        !--Create a list for LPs
           TYPE(LPLinkeListType) :: LPList
   
            CONTAINS
       
        !--Initialize the List
        SUBROUTINE InitalizeList(This)
            CLASS(LPLinkeListType), INTENT(INOUT):: This
            This%Head=>NULL()
            This%Tail=>NULL()
            This%TotNumParcels=0     
        END SUBROUTINE InitalizeList
        
        !--Insert parcels to the list 
        SUBROUTINE InsertNode (This, LPid, xcrd, ycrd, zcrd, ielem, iStat) 
            CLASS(LPLinkeListType), INTENT(INOUT):: This
            !TYPE(RELEASE_TYPE), INTENT(IN)::SRC
            INTEGER, INTENT (IN)            :: LPid, ielem
            INTEGER,INTENT(OUT)          :: iStat
            REAL(KIND=8), INTENT(IN)  :: xcrd, ycrd, zcrd
            
            TYPE (LPNodeType), POINTER ::new_node
            INTEGER::     new_id
            
            iStat=0
            !--Allocate New Node for LP
            ALLOCATE(new_node)
            
            IF(.NOT. ASSOCIATED( This%Tail) ) THEN
                 new_id=1
             ELSE
                 new_id=This%Tail%ID+1
            ENDIF
            
            new_node%next=>NULL ()
            new_node%ID=new_id
            !--Insert new node the list at the end using the tail pointer
            IF (.NOT. ASSOCIATED(This%Head)) THEN
                    This%Head=>new_node
                    This%Tail=>new_node
                    This%Head%Xcord=xcrd
                    This%Head%Ycord=ycrd
                    This%Head%Zcord=zcrd
                    This%Head%Elem=ielem  
                    This%Head%Ingrid=.TRUE.           ! within domain Ingrid=T; Outside the domain Ingrid=F         
            ELSE
                This%Tail%Next=>new_node
                This%Tail=>new_node
                    This%Tail%Xcord=xcrd
                    This%Tail%Ycord=ycrd
                    This%Tail%Zcord=zcrd
                    This%Tail%Elem=ielem 
                    This%Tail%Ingrid=.TRUE.          ! within domain Ingrid=T; Outside the domain Ingrid=F                              
            ENDIF

            This%TotNumParcels=This%TotNumParcels+1
            new_id = This%Tail%ID+1
            
        END SUBROUTINE InsertNode
    
   
     !--Remove a LP from the list
        SUBROUTINE RemoveNode (This, ParcelID, found, iStat)
            CLASS(LPLinkeListType), INTENT(INOUT):: This
            INTEGER, INTENT (IN) :: ParcelID
            INTEGER,INTENT(OUT)          :: iStat
            LOGICAL, INTENT(OUT)  :: found
            
            TYPE (LPNodeType), POINTER ::current, previous
            INTEGER :: icount
            
            iStat=0
            
            IF(.NOT. ASSOCIATED(This%Head) ) RETURN
            
            current=>This%Head
            previous=>NULL()
            
            !--Special case: Delete head node
            IF (ParcelID==1) THEN
                This%Head=>current%next
                IF (.NOT. ASSOCIATED (This%Head)) This%Tail=>NULL() ! List now empty
                DEALLOCATE(current)
                This%TotNumParcels=This%TotNumParcels-1
                This%NumParDelete=This%NumParDelete+1                
                found=.TRUE.
                RETURN
            ENDIF
            
            !--Traverse to the specifed node and remove it
            DO WHILE (ASSOCIATED(current) .AND. current%ID< ParcelID)
                previous=>current
                current=>current%Next
                icount=icount+1
            ENDDO
            
            IF (ASSOCIATED(current) .AND. current%ID==ParcelID) THEN
                previous%next=>current%next
                IF (.NOT.ASSOCIATED(previous%next) )This%Tail=> previous ! update tail if the last node was deleted
                DEALLOCATE(current)
                This%TotNumParcels=This%TotNumParcels-1
                This%NumParDelete=This%NumParDelete+1
                found=.TRUE.
            ENDIF

        END SUBROUTINE RemoveNode
        !
    
    END MODULE LinkedList_Mod

    !---------------------------------------------------------------------------------------------------------------------------!    
!    
!    MODULE LinkedList_Mod
!    IMPLICIT none
!    PRIVATE
!    PUBLIC :: LinkedListType, NodeType, LPList, LPBeginLink!, LinkBeginLP
!    INTEGER:: LPBeginLink
!    
!    TYPE :: NodeType
!        INTEGER:: Data
!        TYPE(NodeType), POINTER ::next => NULL()
!    END  TYPE  NodeType
!    
!    TYPE :: LinkedListType
!        !INTEGER, ALLOCATABLE         :: BeginID(:)
!        TYPE(NodeType), POINTER :: head => NULL()
!    CONTAINS
!        PROCEDURE :: Create=>CreateList
!        PROCEDURE :: Insert => InsertNode
!        PROCEDURE :: Delete => DeleteNode
!        PROCEDURE :: print_list => Print_Nodes
!       ! PROCEDURE :: print_begID=> Print_BegID
!
!    END TYPE LinkedListType
!   
!    !--Create a list for LPs
!  TYPE(LinkedListType) :: LPList
!  
!  !--Array to save begining index
!    !TYPE(NodeType) :: LPBeginLink
!    
!    CONTAINS
!    
!    !!--Creat a list for LPs
!    SUBROUTINE CreateList(This, nLPs)
!        CLASS( LinkedListType), INTENT(INOUT) :: This
!        INTEGER, INTENT (IN) :: nLPs
!        TYPE (NodeType), POINTER :: New_Node, current
!       INTEGER     ::i
!       
!       !-Create the head node
!       ALLOCATE(This%head)
!       
!       current=>This%head
!       
!       !-Create the links for the rest of the nodes
!        DO I =1, nLPs
!            ALLOCATE(New_Node)
!            New_Node%next=>NULL()
!            current%Next=>New_Node
!            current=>New_node    
!        ENDDO
!    
!    END SUBROUTINE CreateList
!    
!    !-Insert a new node at the end of the list
!    SUBROUTINE InsertNode( This, value, iStat)
!     CLASS( LinkedListType), INTENT(INOUT) :: This
!     INTEGER, INTENT (IN)  :: value
!     INTEGER,INTENT(OUT)          :: iStat
!    
!     TYPE (NodeType), POINTER :: New_Node, current
!     
!     !Initialize
!    iStat = 0
!
!
!     ALLOCATE (new_node)
!     New_Node%Data= value
!     New_Node%next=>NULL()
!     
!     IF (.NOT. ASSOCIATED(This%head)) THEN
!         This%head=>New_Node
!      ELSE
!         CURRENT => This%head
!         DO WHILE (ASSOCIATED(Current%Next))
!             current =>Current%next
!         ENDDO
!             LPBeginLink=New_Node%Data
!             Current%next=> New_Node
!      ENDIF
!
!     ! LPList%BeginID(value)=value     ! Need to fix this
!
!    END SUBROUTINE InsertNode
!    
!    !-Delete a node
!    SUBROUTINE DeleteNode ( This, Position, Found)
!    CLASS (LinkedListType), INTENT (INOUT):: This
!    INTEGER, INTENT(IN) :: Position
!    LOGICAL                           :: Found
!    TYPE(NodeType), POINTER :: current, prev
!    INTEGER                          :: icount
!    
!        found=.false.
!    
!        IF (.NOT. ASSOCIATED(This%head)) RETURN
!    
!        !-Special Case-Delete Head node
!        IF (Position==1) THEN
!            current=>This%Head
!            This%Head=>current%Next
!            DEALLOCATE(current)
!            Found=.TRUE.
!            RETURN
!        ENDIF
!   
!        !-Travese to the node at the specified position
!        current=>This%Head
!        prev=>NULL()
!        icount=1
!        DO WHILE (ASSOCIATED(current) .AND. icount <position)
!            prev=>current
!            current=>current%next
!            icount=icount+1
!        ENDDO
!    
!        IF (ASSOCIATED(current) .AND. icount==position) THEN
!            prev%next=>current%next
!            DEALLOCATE(current)
!            Found=.TRUE.
!        ENDIF
!
!
!    END SUBROUTINE DeleteNode
!    
!    
!    !-Print all the nodes in the list
!    SUBROUTINE Print_Nodes( This) 
!        CLASS (LinkedListType), INTENT (INOUT):: This
!        TYPE(NodeType), POINTER :: current
!
!        current=>This%head
!        IF( .NOT. ASSOCIATED(current) ) THEN
!            PRINT *, "List is empty"
!            RETURN
!        ENDIF
!        
!        WRITE (*,*) "Printing Linked List Contents"
!        DO WHILE (ASSOCIATED(current))
!            WRITE(*,*) current%data
!            current=>current%next
!        ENDDO
!    END  SUBROUTINE Print_Nodes
!    
!    
!    !    !-Print all the nodes in the list
!    !SUBROUTINE Print_BegID( This) 
!    !    CLASS (LinkedListType), INTENT (INOUT):: This
!    !    TYPE(NodeType), POINTER :: current
!    !
!    !    current=>This%head
!    !    IF( .NOT. ASSOCIATED(current) ) THEN
!    !        PRINT *, "List is empty"
!    !        RETURN
!    !    ENDIF
!    !    
!    !    WRITE (*,*) "Printing Linked List Contents"
!    !    DO WHILE (ASSOCIATED(current))
!    !        WRITE(*,*) current%data
!    !        current=>current%next
!    !    ENDDO
!    !END  SUBROUTINE Print_Nodes
!
!    END MODULE LinkedList_Mod
!!---------------------------------------------------------------------------------------------------------------------------!    
    
    
 
!----------------------------------------------------------------------------------------------------------------------------!
!--Older version of the linked-list
    
    SUBROUTINE Make_List(npar)
!
!       create a list for multi-species problem, and the link(0)
!       correspond to the top of the avail list
!
!       beg(i), i=1: sus par
!                     i=2: sus_lnd par
!             
       
	  USE AdvDispLagrPar_Mod
      
      
        DO  10 i=0,npar  !+nbuf
                 Link(i)=mod(i+1, npar+1)
10      CONTINUE 
        Link(npar)=npar+1     ! added 1/19/2024
        
        
       DO 20 i=1,n_partype 
                beg(i)=-1
20     CONTINUE
       
       RETURN
    END SUBROUTINE   Make_List 
    
    
  SUBROUTINE List_Insert(ista, new)

       !insert a new particle to the top of ista-species,
       !the subscript of the new particle is equal to "nnew".

        USE AdvDispLagrPar_Mod
        INTEGER avail, loc

        avail=Link(0)      ! the top of avail list
        loc=beg(ista)
        new=avail
        avail=Link(avail)
       IF (avail.eq.0) THEN
           WRITE (*,*) 'Particle numbers overpass the maximum when the', ista, 'th kind of particles are inserted.'
           PAUSE
        ENDIF
        Link(new)=loc
        beg(ista)=new
        Link(0)=avail
        
        RETURN
        
        END SUBROUTINE List_Insert
!!
!!
!!
   
        SUBROUTINE List_Delete(ista,loc)
!
!       delete a particle whose subscript is "loc" from
!       ista-species.
!
        USE AdvDispLagrPar_Mod
        INTEGER avail, loc, locp
!
        avail=Link(0)      ! the top of avail list
        if(loc.eq.beg(ista)) then
           beg(ista)=Link(loc)
        else
           locp=beg(ista)
           do 10, while(Link(locp) .ne. loc)
              locp=Link(locp)
10         continue
           Link(locp)=Link(loc)
        end if
        Link(loc)=avail
        Link(0)=loc
        return
        
        
    END SUBROUTINE List_Delete
    
 !---------------------------------------------------------------------------------------------------------------------------!    

    MODULE GenericLinkedList

  IMPLICIT NONE
  
  
  
! ******************************************************************
! ******************************************************************
! ******************************************************************
! ***
! *** VARIABLE DEFINITIONS
! ***
! ******************************************************************
! ******************************************************************
! ******************************************************************


  ! -------------------------------------------------------------
  ! --- GENERIC NODE TYPE
  ! -------------------------------------------------------------
  TYPE LLNodeType
      CLASS(*),ALLOCATABLE     :: Value
      TYPE(LLNodeType),POINTER :: pNext  => NULL()
  END TYPE LLNodeType
  
  
  ! -------------------------------------------------------------
  ! --- GENERIC LIST TYPE
  ! -------------------------------------------------------------
  TYPE,ABSTRACT :: GenericLinkedListType
       !PRIVATE
      INTEGER                  :: iNNodes  =  0        !Number of nodes in the list
      TYPE(LLNodeType),POINTER :: pHead    => NULL()   !Head of the list
      TYPE(LLNodeType),POINTER :: pTail    => NULL()   !Tail of the list
      TYPE(LLNodeType),POINTER :: pCurrent => NULL()   !Current node in the list
  CONTAINS
      PROCEDURE,NON_OVERRIDABLE,PASS :: AddNode         
      PROCEDURE,NON_OVERRIDABLE,PASS :: GetNNodes       
      PROCEDURE,NON_OVERRIDABLE,PASS :: GetCurrentValue 
      PROCEDURE,NON_OVERRIDABLE,PASS :: Reset           
      PROCEDURE,NON_OVERRIDABLE,PASS :: Next            
      PROCEDURE,NON_OVERRIDABLE,PASS :: Delete          
      PROCEDURE,NON_OVERRIDABLE,PASS :: GetArray        => GenericLinkedList_ConvertToIntegerArray
  END TYPE GenericLinkedListType
  
  
 ! -------------------------------------------------------------
  ! --- LAGRANGIAN PARCEL LINKED LIST TYPE
  ! -------------------------------------------------------------
  TYPE,EXTENDS(GenericLinkedListType) :: LPListType
  END TYPE LPListType


  
  ! -------------------------------------------------------------
  ! --- MISC. ENTITIES
  ! -------------------------------------------------------------
  INTEGER,PARAMETER                   :: ModNameLen = 19
  CHARACTER(LEN=ModNameLen),PARAMETER :: ModName    = 'GenericLinkedList::'
  
    
  
  
CONTAINS
    
    
    
    
! ******************************************************************
! ******************************************************************
! ******************************************************************
! ***
! *** DESTRUCTORS
! ***
! ******************************************************************
! ******************************************************************
! ******************************************************************

  ! -------------------------------------------------------------
  ! --- DELETE LIST 
  ! -------------------------------------------------------------
  SUBROUTINE Delete(List)
    CLASS(GenericLinkedListType),TARGET :: List
    
    !Local variables
    INTEGER                  :: indx
    TYPE(LLNodeType),POINTER :: pNext
    
    !Return if list is empty
    IF (List%iNNodes .EQ. 0) RETURN
    
    !Delete nodes one by one
    DO indx=1,List%iNNodes-1
        pNext => List%pHead%pNext
        DEALLOCATE (List%pHead%Value)
        NULLIFY (List%pHead%pNext)
        DEALLOCATE (List%pHead)
        List%pHead => pNext
    END DO
    DEALLOCATE (List%pHead%Value)
    NULLIFY (List%pHead%pNext)
    DEALLOCATE (List%pHead)
    NULLIFY (List%pTail)
    NULLIFY (List%pCurrent)
    
    !Set the number of data nodes to zero
    List%iNNodes = 0
    
  END SUBROUTINE Delete
  
  
  



! ******************************************************************
! ******************************************************************
! ******************************************************************
! ***
! *** GETTERS
! ***
! ******************************************************************
! ******************************************************************
! ******************************************************************

  ! -------------------------------------------------------------
  ! --- GET NUMBER OF NODES IN LIST
  ! -------------------------------------------------------------
  PURE FUNCTION GetNNodes(List) RESULT(iNNodes)
    CLASS(GenericLinkedListType),INTENT(IN) :: List
    INTEGER                                 :: iNNodes
    
    iNNodes = List%iNNodes
    
  END FUNCTION GetNNodes
  
  
  ! -------------------------------------------------------------
  ! --- GET DATA STORED IN THE NODE POINTED BY pCurrent
  ! -------------------------------------------------------------
  FUNCTION GetCurrentValue(List) RESULT(pValue)
    CLASS(GenericLinkedListType),TARGET,INTENT(IN) :: List
    CLASS(*),POINTER                               :: pValue
    
    IF (ASSOCIATED(List%pCurrent)) THEN
        ALLOCATE (pValue , SOURCE=List%pCurrent%Value)
        pValue => List%pCurrent%Value
    ELSE
        NULLIFY(pValue)
    ENDIF

  END FUNCTION GetCurrentValue 
  
  
  

! ******************************************************************
! ******************************************************************
! ******************************************************************
! ***
! *** MISC. METHODS
! ***
! ******************************************************************
! ******************************************************************
! ******************************************************************

  ! -------------------------------------------------------------
  ! --- CONVERT AN INTEGER LINKED-LIST TO ARRAY
  ! -------------------------------------------------------------
  SUBROUTINE GenericLinkedList_ConvertToIntegerArray(List,iArray,iStat)
    CLASS(GenericLinkedListType),TARGET,INTENT(IN) :: List
    INTEGER,ALLOCATABLE,INTENT(INOUT)              :: iArray(:)
    INTEGER,INTENT(OUT)                            :: iStat
    
    !Local variables
    CHARACTER(LEN=ModNameLen+39) :: ThisProcedure = ModName // 'GenericLinkedList_ConvertToIntegerArray'
    INTEGER                      :: indx,iErrorCode
    CHARACTER                    :: cErrorMsg*200
    
    !Initialize
    iStat = 0
    
    !Deallocate array in case it is already allocated
    DEALLOCATE (iArray ,STAT=iErrorCode)
    
    !Make sure that the list is not empty
    IF (List%iNNodes .EQ. 0) THEN
        ALLOCATE (iArray(0))
        RETURN
    END IF
    
    !Allocate return array
    ALLOCATE (iArray(List%iNNodes) , STAT=iErrorCode , ERRMSG=cErrorMsg)
    IF (iErrorCode .NE. 0) THEN
        !CALL SetLastMessage('Error in allocating memory to convert a linked list to an integer array.'//NEW_LINE('x')//TRIM(cErrorMsg),f_iFatal,ThisProcedure)
        iStat = -1
        RETURN
    END IF
    
    !Store integer linked list in the return array
    CALL List%Reset()
    DO indx=1,List%iNNodes
        SELECT TYPE (pValue => List%pCurrent%Value)
           TYPE IS(INTEGER)
              iArray(indx) = pValue 
        END SELECT
        CALL List%Next()
    END DO
    
  END SUBROUTINE GenericLinkedList_ConvertToIntegerArray
  
  
  ! -------------------------------------------------------------
  ! --- ADD A NODE TO THE LIST
  ! -------------------------------------------------------------
  SUBROUTINE AddNode(List,ValueToStore,iStat)
    CLASS(GenericLinkedListType) :: List
    CLASS(*),INTENT(IN)          :: ValueToStore
    INTEGER,INTENT(OUT)          :: iStat
    
    !Initialize
    iStat = 0
    
    
    IF (List%iNNodes .EQ. 0) THEN
        ALLOCATE (List%pHead)
        ALLOCATE (List%pHead%Value , SOURCE=ValueToStore)
        List%pTail => List%pHead
    ELSE
        ALLOCATE (List%pTail%pNext)
        List%pTail => List%pTail%pNext
        ALLOCATE (List%pTail%Value , SOURCE=ValueToStore)
    END IF
    
    List%pCurrent => List%pTail
    List%iNNodes  =  List%iNNodes + 1
    
    
  END SUBROUTINE AddNode
    

  ! -------------------------------------------------------------
  ! --- RESET THE LIST 
  ! -------------------------------------------------------------
  SUBROUTINE Reset(List)
    CLASS(GenericLinkedListType) :: List
    
    List%pCurrent => List%pHead
    List%pHead => NULL()
    List%pTail => NULL()
    List%iNNodes = 0

    
  END SUBROUTINE Reset
  
  

  
  
  ! -------------------------------------------------------------
  ! --- MOVE CURENT POINTER TO NEXT NODE IN LIST 
  ! --- Pre-condition: Initial pCurrent must not point to pTail
  ! -------------------------------------------------------------
  SUBROUTINE Next(List)
    CLASS(GenericLinkedListType) :: List
    
    List%pCurrent => List%pCurrent%pNext
    
  END SUBROUTINE Next

 
    END MODULE    GenericLinkedList
!---------------------------------------------------------------------------------------------------------------------------!    
    
    

