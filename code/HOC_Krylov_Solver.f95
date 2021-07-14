 !     Last change:  IMA  23 Jul 2010   10:30 am
MODULE PARAM
IMPLICIT NONE
INTEGER,PARAMETER                :: IM=33
INTEGER,PARAMETER                :: JM=33
INTEGER,PARAMETER                :: ANCHOR=IM*((JM-1)/2-1)+2
INTEGER,PARAMETER                :: NX=IM
INTEGER,PARAMETER                :: NY=JM
INTEGER,PARAMETER                :: NODE_NUM=IM*JM
INTEGER,PARAMETER                :: NRD=IM*JM            !RHS dimension for Neumann type boundary
INTEGER,PARAMETER                :: DRD=(IM-2)*(JM-2)    !RHS dimension for Dirichlet type boundary
INTEGER,PARAMETER                :: ITOL=1
INTEGER,PARAMETER                :: ITMAX=1000
REAL(8),PARAMETER                :: PI=3.141592653589793238462643
REAL(8),PARAMETER                :: LX=2.D0*PI
REAL(8),PARAMETER                :: LY=2.D0*PI
REAL(8),PARAMETER                :: DXI=LX/(IM-1)
REAL(8),PARAMETER                :: DX=DXI
REAL(8),PARAMETER                :: DETA=LY/(JM-1)
REAL(8),PARAMETER                :: DXI2=DXI*DXI
REAL(8),PARAMETER                :: DETA2=DETA*DETA
REAL(8),PARAMETER                :: PLUS=DETA2+DXI2
REAL(8),PARAMETER                :: BET=DXI/DETA
REAL(8),PARAMETER                :: EPS=1.D-6
REAL(8),PARAMETER                :: TOL=1.D-11!MAX(EPS*(DX**6),1.D-16)
INTEGER,PARAMETER                :: MAXITS=100000
LOGICAL                          :: CONSTITUTE_AND_SAVE_MATRIX=.FALSE.
LOGICAL                          :: SAVE_MATRIX_PATTERN=.FALSE.
LOGICAL                          :: CONSTITUTE_DENSE_MATRIX_AND_EIGEN_STRUCTURE=.FALSE.
LOGICAL                          :: SAVE_ITERATION_NORM=.TRUE.
END MODULE PARAM
!******************************************************************************
!******************************************************************************
!******************************************************************************
PROGRAM MAIN
USE PARAM
IMPLICIT NONE
INTEGER                          :: INODE,I,J,K,MD,NNZERO
INTEGER,DIMENSION(DRD)           :: RHS_POINTER
REAL(8),DIMENSION(DRD)           :: RHS
REAL(8),DIMENSION(IM,JM)         :: X,Y,XI,ETA,U
REAL(8),DIMENSION(NODE_NUM)      :: XUS,YUS,T,F,G,COEFF_BL,COEFF_B,COEFF_BR,COEFF_L,COEFF_C,COEFF_R,&
                                    COEFF_TL,COEFF_T,COEFF_TR,DXR,DXL,DYU,DYD
LOGICAL,DIMENSION(NODE_NUM)      :: BOUNDARY_TYPE,BL_CORNER,BR_CORNER,TL_CORNER,TR_CORNER,T_BMC,L_BMC,R_BMC,B_BMC
INTEGER,DIMENSION(NODE_NUM)      :: BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,T_POINT,TR_POINT,TT_POINT,&
                                    RR_POINT,BB_POINT,LL_POINT
REAL(8),ALLOCATABLE              :: AA_AUX(:),ACO(:)!,DNS(:,:)
INTEGER,ALLOCATABLE              :: JA_AUX(:),IA_AUX(:),JA_AUX_POINTER(:),ICO(:),JCO(:)
!------------------------------------------------------------------------------
!PRINT*, TOL
CALL STRUCT_GRID          (XI,ETA)
CALL UNSTRUCT_GRID        (XI,ETA,XUS,YUS,BOUNDARY_TYPE,BL_CORNER,BR_CORNER,TL_CORNER,TR_CORNER,T_BMC,L_BMC,R_BMC,B_BMC)
CALL NEIGHBORS            (BOUNDARY_TYPE,TT_POINT,RR_POINT,BB_POINT,LL_POINT,BL_POINT,B_POINT,BR_POINT,L_POINT,&
                           R_POINT,TL_POINT,T_POINT,TR_POINT,BL_CORNER,BR_CORNER,TL_CORNER,TR_CORNER,T_BMC,L_BMC,R_BMC,B_BMC)
CALL DIRICHLET_BC_FORCE_FUNC         (XUS,YUS,BOUNDARY_TYPE,T,F)
CALL COMPACT_FDM_RHS (BOUNDARY_TYPE,TT_POINT,RR_POINT,BB_POINT,LL_POINT,L_POINT,R_POINT,T_POINT,&
                      B_POINT,BL_POINT,BR_POINT,TL_POINT,TR_POINT,F,G)
!------------------------------------------------------------------------------
CALL RHS_CONSTITUTION_DIRICHLET (G,COEFF_BL,COEFF_B,COEFF_BR,COEFF_L,COEFF_C,COEFF_R,&
COEFF_TL,COEFF_T,COEFF_TR,BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,T_POINT,&
TR_POINT,BOUNDARY_TYPE,T,RHS_POINTER,RHS,MD,NNZERO)
!------------------------------------------------------------------------------
ALLOCATE (AA_AUX(NNZERO),JA_AUX_POINTER(NNZERO),JA_AUX(NNZERO),ACO(NNZERO),JCO(NNZERO),ICO(NNZERO),IA_AUX(MD+1))!,DNS(MD,MD)
!------------------------------------------------------------------------------
IF (CONSTITUTE_AND_SAVE_MATRIX) THEN
CALL MATRIX_CONSTITUTION_DIRICHLET (BOUNDARY_TYPE,NNZERO,MD,COEFF_BL,COEFF_B,COEFF_BR,COEFF_L,COEFF_C,COEFF_R,&
COEFF_TL,COEFF_T,COEFF_TR,BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,T_POINT,&
TR_POINT,AA_AUX,JA_AUX_POINTER,IA_AUX,JA_AUX)
CALL SAVE_CSR_MATRIX (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX)
ELSE
CALL READ_CSR_MATRIX (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX)
END IF
!------------------------------------------------------------------------------
!IF (SAVE_MATRIX_PATTERN) THEN
!CALL MATRIX_PATTERN_2 (MD,2,AA_AUX,JA_AUX,IA_AUX,NNZERO,ACO,ICO,JCO)
!END IF
!IF (CONSTITUTE_DENSE_MATRIX_AND_EIGEN_STRUCTURE) THEN
!CALL MATRIX_PATTERN (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,DNS)
!END IF
!==============================================================================
CALL SOLVERS (BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,T_POINT,TR_POINT,&
              MD,NNZERO,G,RHS,RHS_POINTER,AA_AUX,JA_AUX,IA_AUX,BOUNDARY_TYPE,T)
CALL WRITE_RESULT         (XI,ETA,T,F)
!------------------------------------------------------------------------------
END PROGRAM MAIN
!******************************************************************************
!******************************************************************************
!******************************************************************************
!UNSTRUCT_GRID Subroutine
!This subroutine which has been written by Iman Farahbakhsh on Jan 2010 enable
!us to extract the node number for each grid point with i and j indices. The
!type of each node (boundary-node or middle-node) is also determined as a
!logical parameter.
!
!on entry
!===========
!X,Y        The positions which are stored by i and j indices
!
!on exit
!===========
!XUS,YUS    The positions which are stored by node index
!
!BOUNDARY_TYPE  A logical array which shows that a node is boundary or not
!
!Last change: Jan 17 2010
!______________________________________________________________________________
SUBROUTINE UNSTRUCT_GRID (X,Y,XUS,YUS,BOUNDARY_TYPE,BL_CORNER,BR_CORNER,TL_CORNER,TR_CORNER,T_BMC,L_BMC,R_BMC,B_BMC)
USE PARAM, ONLY            : IM,JM,NODE_NUM
IMPLICIT NONE
INTEGER                          :: INODE,I,J
REAL(8),DIMENSION(IM,JM)         :: X,Y
REAL(8),DIMENSION(NODE_NUM)      :: XUS,YUS
LOGICAL,DIMENSION(NODE_NUM)      :: BOUNDARY_TYPE,BL_CORNER,BR_CORNER,TL_CORNER,TR_CORNER,T_BMC,L_BMC,R_BMC,B_BMC
!------------------------------------------------------------------------------
BL_CORNER=.FALSE.
BR_CORNER=.FALSE.
TL_CORNER=.FALSE.
TR_CORNER=.FALSE.
B_BMC=.FALSE.
L_BMC=.FALSE.
R_BMC=.FALSE.
T_BMC=.FALSE.
DO J=1,JM
   DO I=1,IM
      IF (J==1) THEN
         IF (I==1) THEN
            INODE=1
            BOUNDARY_TYPE(INODE)=.TRUE.
            BL_CORNER(INODE)=.TRUE.
            XUS(INODE)=X(I,J)
            YUS(INODE)=Y(I,J)
         ELSE IF (I==IM) THEN
            INODE=IM
            BOUNDARY_TYPE(INODE)=.TRUE.
            BR_CORNER(INODE)=.TRUE.
            XUS(INODE)=X(I,J)
            YUS(INODE)=Y(I,J)
         ELSE
            INODE=I
            BOUNDARY_TYPE(INODE)=.TRUE.
            B_BMC(INODE)=.TRUE.
            XUS(INODE)=X(I,J)
            YUS(INODE)=Y(I,J)
         END IF
      ELSE IF (J==JM) THEN
         IF (I==1) THEN
            INODE=IM*(JM-1)+1
            BOUNDARY_TYPE(INODE)=.TRUE.
            TL_CORNER(INODE)=.TRUE.
            XUS(INODE)=X(I,J)
            YUS(INODE)=Y(I,J)
         ELSE IF (I==IM) THEN
            INODE=IM*JM
            BOUNDARY_TYPE(INODE)=.TRUE.
            TR_CORNER(INODE)=.TRUE.
            XUS(INODE)=X(I,J)
            YUS(INODE)=Y(I,J)
         ELSE
            INODE=IM*(JM-1)+I
            BOUNDARY_TYPE(INODE)=.TRUE.
            T_BMC(INODE)=.TRUE.
            XUS(INODE)=X(I,J)
            YUS(INODE)=Y(I,J)
         END IF
      ELSE
         IF (I==1) THEN
            INODE=(J-1)*IM+1
            BOUNDARY_TYPE(INODE)=.TRUE.
            L_BMC(INODE)=.TRUE.
            XUS(INODE)=X(I,J)
            YUS(INODE)=Y(I,J)
         ELSE IF (I==IM) THEN
            INODE=J*IM
            BOUNDARY_TYPE(INODE)=.TRUE.
            R_BMC(INODE)=.TRUE.
            XUS(INODE)=X(I,J)
            YUS(INODE)=Y(I,J)
         ELSE
            INODE=(J-1)*IM+I
            BOUNDARY_TYPE(INODE)=.FALSE.
            XUS(INODE)=X(I,J)
            YUS(INODE)=Y(I,J)
         END IF
      END IF
   END DO
END DO
!------------------------------------------------------------------------------
!WRITE(*,*) '******* Nodes number were assigned *******'
!------------------------------------------------------------------------------
END SUBROUTINE UNSTRUCT_GRID
!******************************************************************************
!******************************************************************************
!******************************************************************************
!STRUCT_GRID Subroutine
!This subroutine which has been written by Iman Farahbakhsh on Oct 2009 enable
!us to generate a uniform grid on a rectangular domain.
!
!on exit
!===========
!X,Y        The positions on the uniform grid which are stored by indices of i,j
!
!Last change: Jan 17 2010
!______________________________________________________________________________
SUBROUTINE  STRUCT_GRID (XI,ETA)
USE PARAM,ONLY             : IM,JM,DXI,DETA
IMPLICIT NONE
INTEGER                          :: I,J
REAL(8),DIMENSION(IM,JM)         :: XI,ETA
!------------------------------------------------------------------------------
!OPEN (UNIT=1,FILE='STRUCT_GRID.PLT',STATUS='UNKNOWN')
!------------------------------------------------------------------------------
XI(1,:)=0.D0
DO J=1,JM
DO I=1,IM-1
XI(I+1,J)=XI(I,J)+DXI
END DO
END DO
ETA(:,1)=0.D0
DO I=1,IM
DO J=1,JM-1
ETA(I,J+1)=ETA(I,J)+DETA
END DO
END DO
!------------------------------------------------------------------------------
!WRITE (1,*) 'VARIABLES=XI,ETA'
!WRITE (1,*) 'ZONE'
!WRITE (1,*) 'F=POINT'
!WRITE (1,*) 'I=',IM
!WRITE (1,*) 'J=',JM
!DO I=1,IM
!DO J=1,JM
!WRITE (1,*) XI(I,J),ETA(I,J)
!END DO
!END DO
!------------------------------------------------------------------------------
!WRITE(*,*) '******* Structured grid was written *******'
!------------------------------------------------------------------------------
END SUBROUTINE STRUCT_GRID
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                               NEIGHBORS
!______________________________________________________________________________
SUBROUTINE NEIGHBORS (BOUNDARY_TYPE,TT_POINT,RR_POINT,BB_POINT,LL_POINT,BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,&
                      T_POINT,TR_POINT,BL_CORNER,BR_CORNER,TL_CORNER,TR_CORNER,T_BMC,L_BMC,R_BMC,B_BMC)
USE PARAM,ONLY                   :  IM,NODE_NUM
IMPLICIT NONE
INTEGER                          :: INODE
LOGICAL,DIMENSION(NODE_NUM)      :: BOUNDARY_TYPE,BL_CORNER,BR_CORNER,TL_CORNER,TR_CORNER,T_BMC,L_BMC,R_BMC,B_BMC
INTEGER,DIMENSION(NODE_NUM)      :: BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,T_POINT,TR_POINT,TT_POINT,&
                                    RR_POINT,BB_POINT,LL_POINT
!------------------------------------------------------------------------------
BL_POINT=0
B_POINT=0
BR_POINT=0
L_POINT=0
R_POINT=0
TL_POINT=0
T_POINT=0
TR_POINT=0
TT_POINT=0
BB_POINT=0
LL_POINT=0
RR_POINT=0
DO INODE=1,NODE_NUM
    IF      (B_BMC(INODE))            THEN
      L_POINT(INODE)=INODE-1
      R_POINT(INODE)=INODE+1
      TL_POINT(INODE)=INODE+IM-1
      T_POINT(INODE)=INODE+IM
      TR_POINT(INODE)=INODE+IM+1
    ELSE IF (L_BMC(INODE))       THEN
      B_POINT(INODE)=INODE-IM
      BR_POINT(INODE)=INODE-IM+1
      R_POINT(INODE)=INODE+1
      T_POINT(INODE)=INODE+IM
      TR_POINT(INODE)=INODE+IM+1
    ELSE IF (T_BMC(INODE))       THEN
      BL_POINT(INODE)=INODE-IM-1
	  B_POINT(INODE)=INODE-IM
	  BR_POINT(INODE)=INODE-IM+1
	  L_POINT(INODE)=INODE-1
	  R_POINT(INODE)=INODE+1
    ELSE IF (R_BMC(INODE))       THEN
      BL_POINT(INODE)=INODE-IM-1
      B_POINT(INODE)=INODE-IM
      L_POINT(INODE)=INODE-1
      TL_POINT(INODE)=INODE+IM-1
      T_POINT(INODE)=INODE+IM
    ELSE IF (BL_CORNER(INODE))       THEN
	  R_POINT(INODE)=INODE+1
	  T_POINT(INODE)=INODE+IM
	  TR_POINT(INODE)=INODE+IM+1
    ELSE IF (TL_CORNER(INODE))       THEN
      B_POINT(INODE)=INODE-IM
	  BR_POINT(INODE)=INODE-IM+1
	  R_POINT(INODE)=INODE+1
    ELSE IF (TR_CORNER(INODE))       THEN
      BL_POINT(INODE)=INODE-IM-1
	  B_POINT(INODE)=INODE-IM
	  L_POINT(INODE)=INODE-1
    ELSE IF (BR_CORNER(INODE))       THEN
	  L_POINT(INODE)=INODE-1
	  TL_POINT(INODE)=INODE+IM-1
	  T_POINT(INODE)=INODE+IM
    ELSE
          BL_POINT(INODE)=INODE-IM-1
          B_POINT(INODE)=INODE-IM
          BR_POINT(INODE)=INODE-IM+1
          L_POINT(INODE)=INODE-1
          R_POINT(INODE)=INODE+1
          TL_POINT(INODE)=INODE+IM-1
          T_POINT(INODE)=INODE+IM
          TR_POINT(INODE)=INODE+IM+1
          TT_POINT(INODE)=INODE+2*IM
          BB_POINT(INODE)=INODE-2*IM
          LL_POINT(INODE)=INODE-2
          RR_POINT(INODE)=INODE+2
	END IF
END DO
!------------------------------------------------------------------------------
!WRITE(*,*) '******* Neighbors were found *******'
!------------------------------------------------------------------------------
END SUBROUTINE NEIGHBORS
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE DIRICHLET_BC_FORCE_FUNC (X,Y,BOUNDARY_TYPE,T,F)
USE PARAM
IMPLICIT NONE
INTEGER :: INODE
REAL(8),DIMENSION(NODE_NUM)             :: T,F,X,Y
LOGICAL,DIMENSION(NODE_NUM)             :: BOUNDARY_TYPE
!------------------------------------------------------------------------------
DO INODE=1,NODE_NUM
   IF (BOUNDARY_TYPE(INODE)) THEN
       T(INODE)=DCOS(X(INODE))*DSIN(3.D0*PI*Y(INODE)/LY)
       !T(INODE)=DCOS(X(INODE)+Y(INODE))
   END IF
       F(INODE)=-(((9.D0*PI**2.D0)/LY**2.D0)+1.D0)*DCOS(X(INODE))*DSIN(3.D0*PI*Y(INODE)/LY)
       !F(INODE)=-2.D0*DCOS(X(INODE)+Y(INODE))
END DO
!------------------------------------------------------------------------------
!WRITE(*,*) '******* Dirichlet BC and force function are made *******'
!------------------------------------------------------------------------------
END SUBROUTINE DIRICHLET_BC_FORCE_FUNC
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE COMPACT_FDM_RHS (BOUNDARY_TYPE,TT_POINT,RR_POINT,BB_POINT,LL_POINT,&
                            L_POINT,R_POINT,T_POINT,B_POINT,BL_POINT,BR_POINT,&
                            TL_POINT,TR_POINT,F,G)
USE PARAM
IMPLICIT NONE
INTEGER :: INODE
REAL(8),DIMENSION(NODE_NUM)   :: G,F
LOGICAL,DIMENSION(NODE_NUM)   :: BOUNDARY_TYPE
INTEGER,DIMENSION(NODE_NUM)   :: L_POINT,R_POINT,T_POINT,B_POINT,BL_POINT,BR_POINT,&
                                 TL_POINT,TR_POINT,TT_POINT,RR_POINT,BB_POINT,LL_POINT
!------------------------------------------------------------------------------
G=0.D0
DO INODE=1,NODE_NUM
    IF (.NOT.BOUNDARY_TYPE(INODE)) THEN
        IF (BOUNDARY_TYPE(TL_POINT(INODE)).OR.BOUNDARY_TYPE(T_POINT(INODE))&
            .OR.BOUNDARY_TYPE(TR_POINT(INODE)).OR.BOUNDARY_TYPE(BL_POINT(INODE))&
            .OR.BOUNDARY_TYPE(B_POINT(INODE)).OR.BOUNDARY_TYPE(BR_POINT(INODE))&
            .OR.BOUNDARY_TYPE(L_POINT(INODE)).OR.BOUNDARY_TYPE(R_POINT(INODE))) THEN
            G(INODE)=0.5D0*DXI2*(F(B_POINT(INODE))+F(L_POINT(INODE))+8.D0*F(INODE)+F(R_POINT(INODE))+F(T_POINT(INODE)))
        ELSE
            G(INODE)=(1.D0/120.D0)*DXI2*(476.D0*F(INODE)&
                                         +56.D0*(F(B_POINT(INODE))+F(L_POINT(INODE))+F(R_POINT(INODE))+F(T_POINT(INODE)))&
                                         +8.D0*(F(BR_POINT(INODE))+F(BL_POINT(INODE))+F(TL_POINT(INODE))+F(TR_POINT(INODE)))&
                                         -3.D0*(F(TT_POINT(INODE))+F(BB_POINT(INODE))+F(LL_POINT(INODE))+F(RR_POINT(INODE))))
        END IF
    END IF
END DO
!------------------------------------------------------------------------------
!WRITE(*,*) '******* Compact FDM RHS function is made *******'
!------------------------------------------------------------------------------
END SUBROUTINE COMPACT_FDM_RHS
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE WRITE_RESULT (X,Y,T,F)
USE PARAM,ONLY         :  IM,JM,NODE_NUM,PI,LY
IMPLICIT NONE
INTEGER                      :: I,J
REAL(8),DIMENSION(IM,JM)     :: X,Y,T,F
CHARACTER*7                  :: ext
CHARACTER*12                 :: fname
!------------------------------------------------------------------------------
       fname='/results/RESULTS.PLT'
!------------------------------------------------------------------------------
  OPEN(10,FILE=fname,STATUS='UNKNOWN')
  WRITE(10,*)'VARIABLES= "X","Y","T","ERROR","EXACT","FORCE"'
  WRITE(10,*)'ZONE'
  WRITE(10,*)'F=POINT'
  WRITE(10,*)'I=',IM
  WRITE(10,*)'J=',JM
!------------------------------------------------------------------------------
      DO J=1,JM
      DO I=1,IM
         WRITE(10,17) X(I,J),Y(I,J),T(I,J),DABS(T(I,J)-DCOS(X(I,J))*DSIN(3.D0*PI*Y(I,J)/LY)),&
                      DCOS(X(I,J))*DSIN(3.D0*PI*Y(I,J)/LY),F(I,J)
      END DO
      END DO
!------------------------------------------------------------------------------
     CLOSE(10)
17   FORMAT(5E14.6)
!------------------------------------------------------------------------------
!      WRITE(*,*) '======================='
!      WRITE(*,*) 'Printing on ',fname
!      WRITE(*,*) '======================='
!------------------------------------------------------------------------------
END SUBROUTINE WRITE_RESULT
!*************************************************************************************************************
!**************************************************************************************  KRYLOV_SS_SOLVERS  **
!*************************************************************************************************************
 SUBROUTINE SOLVERS (BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,&
                     T_POINT,TR_POINT,MD,NNZERO,G,RHS,RHS_POINTER,AA_AUX,&
                     JA_AUX,IA_AUX,BOUNDARY_TYPE,T)
 USE PARAM
 IMPLICIT NONE
 INTEGER                                 :: NNZERO,I,MD
 REAL(8)                                 :: BEGIN_TIME,END_TIME,FIX
 REAL(8),DIMENSION      (NODE_NUM)       :: T,G,DXR,DXL,DYU,DYD
 LOGICAL,DIMENSION      (NODE_NUM)       :: BOUNDARY_TYPE
 INTEGER,DIMENSION      (NODE_NUM)       :: BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,T_POINT,TR_POINT
 REAL(8),DIMENSION      (DRD)            :: RHS
 INTEGER,DIMENSION      (DRD)            :: RHS_POINTER
!REAL(8),DIMENSION      ((IM-2)*(JM-2),(IM-2)*(JM-2)):: DNS
 REAL(8),DIMENSION      (NODE_NUM)       :: COEFF_C,COEFF_B,COEFF_L,COEFF_R,COEFF_T
 REAL(8),DIMENSION(MD)                   :: X0,R0,SOL
 REAL(8),DIMENSION(NNZERO)               :: AA_AUX,ACO
 INTEGER,DIMENSION(NNZERO)               :: JA_AUX_POINTER,JA_AUX,JCO,ICO
 INTEGER,DIMENSION(MD+1)                 :: IA_AUX
!------------------------------------------------------------------------------
CALL INITIAL_GUESS (MD,X0)
!-------------------------
CALL INITIAL_RESIDUAL (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,RHS,X0,R0)
!==============================================================================
! Krylov subspaces methods
!==============================================================================
!CALL CPU_TIME (BEGIN_TIME)
!===========================
!CALL CG                 (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,X0,R0,SOL)
!CALL ILU0_PCG           (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,X0,R0,SOL)
!CALL CGS                (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,X0,R0,SOL)
!CALL ILU0_PCGS          (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,X0,R0,SOL)
!CALL BCGSTAB            (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,X0,R0,SOL)
CALL ILU0_PBCGSTAB      (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,X0,R0,SOL)
!===========================
!CALL CPU_TIME (END_TIME)
!PRINT*,'Total spent time in Krylov solver is',(END_TIME-BEGIN_TIME)/1.5625D-2,'sec'
!---------------------------
DO I=1,MD
T(RHS_POINTER(I))=SOL(I)
END DO
!==============================================================================
!------------------------------------------------------------------------------
END SUBROUTINE SOLVERS
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                           RHS_CONSTITUTION
!------------------------------------------------------------------------------
SUBROUTINE RHS_CONSTITUTION_DIRICHLET (G,COEFF_BL,COEFF_B,COEFF_BR,COEFF_L,COEFF_C,COEFF_R,&
COEFF_TL,COEFF_T,COEFF_TR,BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,T_POINT,&
TR_POINT,BOUNDARY_TYPE,T,RHS_POINTER,RHS,MD,NNZERO)
USE PARAM
IMPLICIT NONE
INTEGER                              :: INODE,K_RHS,NNZERO,MD
INTEGER,DIMENSION(NODE_NUM)          :: BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,T_POINT,TR_POINT
REAL(8),DIMENSION(NODE_NUM)          :: COEFF_BL,COEFF_B,COEFF_BR,COEFF_L,COEFF_C,COEFF_R,COEFF_TL,COEFF_T,COEFF_TR
LOGICAL,DIMENSION(NODE_NUM)          :: BOUNDARY_TYPE
REAL(8),DIMENSION(NODE_NUM)          :: T,G
REAL(8),DIMENSION(DRD)               :: RHS
INTEGER,DIMENSION(DRD)               :: RHS_POINTER
!------------------------------------------------------------------------------
COEFF_BL=0
COEFF_B=0
COEFF_BR=0
COEFF_L=0
COEFF_C=0
COEFF_R=0
COEFF_TL=0
COEFF_T=0
COEFF_TR=0
K_RHS=0
NNZERO=0
DO INODE=1,NODE_NUM
!------------------------------------------------------------------------------
! Determines whether a node is boundary or not
!------------------------------------------------------------------------------
   IF (.NOT.BOUNDARY_TYPE(INODE)) THEN
!------------------------------------------------------------------------------
	     K_RHS=K_RHS+1
!------------------------------------------------------------------------------
! Determines the RHS before any change
!------------------------------------------------------------------------------
	     RHS(K_RHS)=G(INODE)
!------------------------------------------------------------------------------
! Determines the variable coefficient at the center node
!------------------------------------------------------------------------------
         COEFF_C(INODE)=-20.D0
!------------------------------------------------------------------------------
! Determines the variable coefficient at the bottom left node
!------------------------------------------------------------------------------
         COEFF_BL(BL_POINT(INODE))=1.D0
!------------------------------------------------------------------------------
! Constitutes the new RHS on the basis whether bottom left node is
! boundary or not
!------------------------------------------------------------------------------
		 IF (BOUNDARY_TYPE(BL_POINT(INODE))) THEN
		 RHS(K_RHS)=RHS(K_RHS)-COEFF_BL(BL_POINT(INODE))*T(BL_POINT(INODE))
		 ELSE
		 NNZERO=NNZERO+1
		 END IF
!------------------------------------------------------------------------------
! Determines the variable coefficient at the bottom node
!------------------------------------------------------------------------------
         COEFF_B(B_POINT(INODE))=4.D0
!------------------------------------------------------------------------------
! Constitutes the new RHS on the basis whether bottom node is
! boundary or not
!------------------------------------------------------------------------------
         IF (BOUNDARY_TYPE(B_POINT(INODE))) THEN
   		 RHS(K_RHS)=RHS(K_RHS)-COEFF_B(B_POINT(INODE))*T(B_POINT(INODE))
		 ELSE
		 NNZERO=NNZERO+1
		 END IF
!------------------------------------------------------------------------------
! Determines the variable coefficient at the bottom right node
!------------------------------------------------------------------------------
         COEFF_BR(BR_POINT(INODE))=1.D0
!------------------------------------------------------------------------------
! Constitutes the new RHS on the basis whether bottom right node is
! boundary or not
!------------------------------------------------------------------------------
         IF (BOUNDARY_TYPE(BR_POINT(INODE))) THEN
   		 RHS(K_RHS)=RHS(K_RHS)-COEFF_BR(BR_POINT(INODE))*T(BR_POINT(INODE))
		 ELSE
		 NNZERO=NNZERO+1
		 END IF
!------------------------------------------------------------------------------
! Determines the variable coefficient at the left node
!------------------------------------------------------------------------------
         COEFF_L(L_POINT(INODE))=4.D0
!------------------------------------------------------------------------------
! Constitutes the new RHS on the basis whether left node is
! boundary or not
!------------------------------------------------------------------------------
         IF (BOUNDARY_TYPE(L_POINT(INODE))) THEN
   		 RHS(K_RHS)=RHS(K_RHS)-COEFF_L(L_POINT(INODE))*T(L_POINT(INODE))
		 ELSE
		 NNZERO=NNZERO+1
		 END IF
!------------------------------------------------------------------------------
! Determines the variable coefficient at the right node
!------------------------------------------------------------------------------
         COEFF_R(R_POINT(INODE))=4.D0
!------------------------------------------------------------------------------
! Constitutes the new RHS on the basis whether right node is
! boundary or not
!------------------------------------------------------------------------------
         IF (BOUNDARY_TYPE(R_POINT(INODE))) THEN
   		 RHS(K_RHS)=RHS(K_RHS)-COEFF_R(R_POINT(INODE))*T(R_POINT(INODE))
		 ELSE
		 NNZERO=NNZERO+1
		 END IF
!------------------------------------------------------------------------------
! Determines the variable coefficient at the top left node
!------------------------------------------------------------------------------
         COEFF_TL(TL_POINT(INODE))=1.D0
!------------------------------------------------------------------------------
! Constitutes the new RHS on the basis whether top left node is
! boundary or not
!------------------------------------------------------------------------------
         IF (BOUNDARY_TYPE(TL_POINT(INODE))) THEN
   		 RHS(K_RHS)=RHS(K_RHS)-COEFF_TL(TL_POINT(INODE))*T(TL_POINT(INODE))
		 ELSE
		 NNZERO=NNZERO+1
		 END IF
!------------------------------------------------------------------------------
! Determines the variable coefficient at the top node
!------------------------------------------------------------------------------
         COEFF_T(T_POINT(INODE))=4.D0
!------------------------------------------------------------------------------
! Constitutes the new RHS on the basis whether top node is
! boundary or not
!------------------------------------------------------------------------------
         IF (BOUNDARY_TYPE(T_POINT(INODE))) THEN
   		 RHS(K_RHS)=RHS(K_RHS)-COEFF_T(T_POINT(INODE))*T(T_POINT(INODE))
		 ELSE
		 NNZERO=NNZERO+1
		 END IF
!------------------------------------------------------------------------------
! Determines the variable coefficient at the top right node
!------------------------------------------------------------------------------
         COEFF_TR(TR_POINT(INODE))=1.D0
!------------------------------------------------------------------------------
! Constitutes the new RHS on the basis whether top right node is
! boundary or not
!------------------------------------------------------------------------------
         IF (BOUNDARY_TYPE(TR_POINT(INODE))) THEN
   		 RHS(K_RHS)=RHS(K_RHS)-COEFF_TR(TR_POINT(INODE))*T(TR_POINT(INODE))
		 ELSE
		 NNZERO=NNZERO+1
		 END IF
         RHS_POINTER(K_RHS)=INODE
!------------------------------------------------------------------------------
   END IF
END DO
NNZERO=NNZERO+K_RHS
MD=K_RHS
!------------------------------------------------------------------------------
!WRITE(*,*) '******* RHS was constituted *******'
!------------------------------------------------------------------------------
END SUBROUTINE RHS_CONSTITUTION_DIRICHLET
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                            MATRIX_CONSTITUTION
!______________________________________________________________________________
SUBROUTINE MATRIX_CONSTITUTION_DIRICHLET (BOUNDARY_TYPE,NNZERO,MD,COEFF_BL,COEFF_B,COEFF_BR,COEFF_L,COEFF_C,COEFF_R,&
COEFF_TL,COEFF_T,COEFF_TR,BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,T_POINT,&
TR_POINT,AA_AUX,JA_AUX_POINTER,IA_AUX,JA_AUX)
USE PARAM
IMPLICIT NONE
INTEGER                         :: INODE,CIR,IAA,ROW,MD,NNZERO,BCOUNTER,K,MBC,J,I,M,NC,II
INTEGER,DIMENSION(MD+1)         :: IA_AUX
INTEGER,DIMENSION(NNZERO)       :: JA_AUX_POINTER,JA_AUX
REAL(8),DIMENSION(NNZERO)       :: AA_AUX
LOGICAL,DIMENSION(NODE_NUM)     :: BOUNDARY_TYPE
INTEGER,DIMENSION(NODE_NUM)     :: BL_POINT,B_POINT,BR_POINT,L_POINT,R_POINT,TL_POINT,T_POINT,TR_POINT
REAL(8),DIMENSION(NODE_NUM)     :: COEFF_BL,COEFF_B,COEFF_BR,COEFF_L,COEFF_C,COEFF_R,COEFF_TL,COEFF_T,COEFF_TR
!------------------------------------------------------------------------------
IAA=0
CIR=0
ROW=1
IA_AUX(1)=1
DO INODE=1,NODE_NUM
IF (.NOT.BOUNDARY_TYPE(INODE)) THEN
   ROW=ROW+1
!------------------------------------------------------------------------------
   IF (.NOT.BOUNDARY_TYPE(BL_POINT(INODE))) THEN
      IAA=IAA+1
      AA_AUX(IAA)=COEFF_BL(BL_POINT(INODE))
      JA_AUX_POINTER(IAA)=BL_POINT(INODE)
      CIR=CIR+1
   END IF
!------------------------------------------------------------------------------
   IF (.NOT.BOUNDARY_TYPE(B_POINT(INODE))) THEN
      IAA=IAA+1
      AA_AUX(IAA)=COEFF_B(B_POINT(INODE))
      JA_AUX_POINTER(IAA)=B_POINT(INODE)
      CIR=CIR+1
   END IF
!------------------------------------------------------------------------------
   IF (.NOT.BOUNDARY_TYPE(BR_POINT(INODE))) THEN
      IAA=IAA+1
      AA_AUX(IAA)=COEFF_BR(BR_POINT(INODE))
      JA_AUX_POINTER(IAA)=BR_POINT(INODE)
      CIR=CIR+1
   END IF
!------------------------------------------------------------------------------
   IF (.NOT.BOUNDARY_TYPE(L_POINT(INODE))) THEN
      IAA=IAA+1
      AA_AUX(IAA)=COEFF_L(L_POINT(INODE))
      JA_AUX_POINTER(IAA)=L_POINT(INODE)
      CIR=CIR+1
   END IF
!------------------------------------------------------------------------------
   IAA=IAA+1
   AA_AUX(IAA)=COEFF_C(INODE)
   JA_AUX_POINTER(IAA)=INODE
   CIR=CIR+1
!------------------------------------------------------------------------------
   IF (.NOT.BOUNDARY_TYPE(R_POINT(INODE))) THEN
      IAA=IAA+1
      AA_AUX(IAA)=COEFF_R(R_POINT(INODE))
      JA_AUX_POINTER(IAA)=R_POINT(INODE)
      CIR=CIR+1
   END IF
!------------------------------------------------------------------------------
   IF (.NOT.BOUNDARY_TYPE(TL_POINT(INODE))) THEN
      IAA=IAA+1
      AA_AUX(IAA)=COEFF_TL(TL_POINT(INODE))
      JA_AUX_POINTER(IAA)=TL_POINT(INODE)
      CIR=CIR+1
   END IF
!------------------------------------------------------------------------------
   IF (.NOT.BOUNDARY_TYPE(T_POINT(INODE))) THEN
      IAA=IAA+1
      AA_AUX(IAA)=COEFF_T(T_POINT(INODE))
      JA_AUX_POINTER(IAA)=T_POINT(INODE)
      CIR=CIR+1
   END IF
!------------------------------------------------------------------------------
   IF (.NOT.BOUNDARY_TYPE(TR_POINT(INODE))) THEN
      IAA=IAA+1
      AA_AUX(IAA)=COEFF_TR(TR_POINT(INODE))
      JA_AUX_POINTER(IAA)=TR_POINT(INODE)
      CIR=CIR+1
   END IF
!------------------------------------------------------------------------------
END IF
IA_AUX(ROW)=CIR+1
END DO
IA_AUX(MD+1)=IA_AUX(1)+NNZERO
!------------------------------------------------------------------------------
DO K=1,NNZERO
MBC=0
DO II=1,NODE_NUM
IF (BOUNDARY_TYPE(II).AND.II<JA_AUX_POINTER(K)) THEN
MBC=MBC+1
END IF
JA_AUX(K)=JA_AUX_POINTER(K)-MBC
END DO
END DO
!------------------------------------------------------------------------------
!WRITE(*,*) '******* CSR format of matrix was constituted *******'
!------------------------------------------------------------------------------
END SUBROUTINE MATRIX_CONSTITUTION_DIRICHLET
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE SAVE_CSR_MATRIX (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX)
USE PARAM
IMPLICIT NONE
INTEGER                    :: I,MD,NNZERO
INTEGER,DIMENSION(MD+1)    :: IA_AUX
INTEGER,DIMENSION(NNZERO)  :: JA_AUX
REAL(8),DIMENSION(NNZERO)  :: AA_AUX
CHARACTER*10               :: EXT
CHARACTER*8                :: FN1
CHARACTER*25               :: FNAME
!------------------------------------------------------------------------------
FN1='COMP_MAT'
WRITE(EXT,'(I7)') IM
FNAME=FN1//EXT//'.DAT'
OPEN(IM,FILE=FNAME,STATUS="UNKNOWN")
DO I=1,NNZERO
    IF (I<=MD+1) THEN
    WRITE(IM,*) AA_AUX(I),JA_AUX(I),IA_AUX(I)
    ELSE
    WRITE(IM,*) AA_AUX(I),JA_AUX(I)
    END IF
END DO
!------------------------------------------------------------------------------
!WRITE(*,*) '******* CSR format of matrix was saved *******'
!------------------------------------------------------------------------------
END SUBROUTINE SAVE_CSR_MATRIX
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE READ_CSR_MATRIX (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX)
USE PARAM
IMPLICIT NONE
INTEGER                    :: I,MD,NNZERO
INTEGER,DIMENSION(MD+1)    :: IA_AUX
INTEGER,DIMENSION(NNZERO)  :: JA_AUX
REAL(8),DIMENSION(NNZERO)  :: AA_AUX
CHARACTER*10               :: EXT
CHARACTER*8                :: FN1
CHARACTER*25               :: FNAME
!------------------------------------------------------------------------------
FN1='COMP_MAT'
WRITE(EXT,'(I7)') IM
FNAME=FN1//EXT//'.DAT'
OPEN(IM,FILE=FNAME,STATUS="UNKNOWN")
DO I=1,NNZERO
    IF (I<=MD+1) THEN
    READ(IM,*) AA_AUX(I),JA_AUX(I),IA_AUX(I)
    ELSE
    READ(IM,*) AA_AUX(I),JA_AUX(I)
    END IF
END DO
!------------------------------------------------------------------------------
!WRITE(*,*) '******* CSR format of matrix was read *******'
!------------------------------------------------------------------------------
END SUBROUTINE READ_CSR_MATRIX
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                     INITIAL_GUESS_SUBROUTINE
!______________________________________________________________________________
SUBROUTINE INITIAL_GUESS (MD,X0)
IMPLICIT NONE
INTEGER                       :: MD
REAL(8),DIMENSION(MD)         :: X0
!------------------------------------------------------------------------------
X0=0.D0
!------------------------------------------------------------------------------
END SUBROUTINE INITIAL_GUESS
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                    INITIAL_RESIDUAL_SUBROUTINE
!______________________________________________________________________________
SUBROUTINE INITIAL_RESIDUAL (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,RHS,X0,R0)
IMPLICIT NONE
INTEGER                      :: MD,NNZERO,I,MEDIUM1,MEDIUM2,K
REAL(8),DIMENSION(MD)        :: RHS,X0,R0,Y0
REAL(8),DIMENSION(NNZERO)    :: AA_AUX
INTEGER,DIMENSION(NNZERO)    :: JA_AUX
INTEGER,DIMENSION(MD+1)      :: IA_AUX
!------------------------------------------------------------------------------
CALL CSR_MAT_V_PRODUCT (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,X0,Y0)
!------------------------------------------------------------------------------
R0=RHS-Y0
!------------------------------------------------------------------------------
END SUBROUTINE INITIAL_RESIDUAL
!******************************************************************************
!******************************************************************************
!******************************************************************************
!This subroutine has been developed in FORTRAN 90 via Iman Farahbakhsh
!on Dec 03 2009 with some adaptions from Numerical Recipes in FORTRAN 77
!______________________________________________________________________________
SUBROUTINE WRITE_EXE (NROT,D,V,NP,A)
IMPLICIT NONE
INTEGER                    :: J,L,KK,LL,NP,NROT,K
REAL(8),DIMENSION(NP,NP)   :: A,V
REAL(8),DIMENSION(NP)      :: R
COMPLEX(8),DIMENSION(NP)   :: D
REAL(8)                    :: RATIO
!------------------------------------------------------------------------------
WRITE (*,'(1X,A,I10)') 'Number of JACOBI rotations: ',NROT
WRITE (*,'(/1X,A)') 'Eigenvalues:'
DO J=1,NP
WRITE (*,'(1X,5F12.6)') D(J)
END DO
WRITE(*,'(/1X,A)') 'Eigenvectors:'
DO J=1,NP
WRITE(*,'(1X,T5,A,I3)') 'Number',J
WRITE(*,'(1X,5F12.6)') (V(K,J),K=1,NP)
END DO
!eigenvector test
WRITE(*,'(/1X,A)') 'Eigenvector Test'
DO J=1,NP
DO L=1,NP
R(L)=0.D0
DO K=1,NP
IF (K.GT.L) THEN
KK=L
LL=K
ELSE
KK=K
LL=L
END IF
R(L)=R(L)+A(LL,KK)*V(K,J)
END DO
END DO
WRITE(*,'(/1X,A,I3)') 'Vector Number',J
WRITE(*,'(/1X,T7,A,T18,A,T31,A)')&
'Vector','Mtrx*Vec.','Ratio'
DO  L=1,NP
RATIO=R(L)/V(L,J)
WRITE(*,'(1X,3F12.6)') V(L,J),R(L),RATIO
END DO
END DO
WRITE(*,*) 'press RETURN to END DO...'
READ(*,*)
!------------------------------------------------------------------------------
END SUBROUTINE WRITE_EXE
!******************************************************************************
!******************************************************************************
!******************************************************************************
!This subroutine has been developed in FORTRAN 90 via Iman Farahbakhsh
!on Dec 03 2009 with some adaptions from Numerical Recipes in FORTRAN 77
!______________________________________________________________________________
SUBROUTINE WRITE_FILE (D,NP,A)
IMPLICIT NONE
INTEGER                             :: I,NP
REAL(8),DIMENSION(NP,NP)            :: A
COMPLEX(8),DIMENSION(NP)            :: D
!------------------------------------------------------------------------------
OPEN (UNIT=1,FILE='EIGEN_VALUE_DISTRIBUTION.PLT',STATUS='UNKNOWN')
WRITE(1,*) 'VARIABLES=REAL,IMAGINARY'
WRITE(1,*) 'ZONE'
WRITE(1,*) 'F=POINT'
DO I=1,NP
WRITE(1,*) REAL(D(I)),IMAG(D(I))
END DO
!------------------------------------------------------------------------------
END SUBROUTINE WRITE_FILE
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                      MATRIX_BY_VECTOR PRODUCT
!                       MATRIX IN CSR_FORMAT
!______________________________________________________________________________
SUBROUTINE CSR_MAT_V_PRODUCT (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,VE,YY)
IMPLICIT NONE
INTEGER                        :: I,MEDIUM1,MEDIUM2,NNZERO,MD
REAL(8),DIMENSION(MD)          :: VE,YY
REAL(8),DIMENSION(NNZERO)      :: AA_AUX
INTEGER,DIMENSION(MD+1)        :: IA_AUX
INTEGER,DIMENSION(NNZERO)      :: JA_AUX
!------------------------------------------------------------------------------
YY=0.D0
DO I=1,MD
MEDIUM1=IA_AUX(I)
MEDIUM2=IA_AUX(I+1)-1
YY(I)=DOT_PRODUCT(AA_AUX(MEDIUM1:MEDIUM2),VE(JA_AUX(MEDIUM1:MEDIUM2)))
END DO
!------------------------------------------------------------------------------
END SUBROUTINE CSR_MAT_V_PRODUCT
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                      TRANSPOSE MATRIX_BY_VECTOR PRODUCT
!                           MATRIX IN CSR_FORMAT
!______________________________________________________________________________
SUBROUTINE TRANSPOSE_MATRIX_VECTOR_PRODUCT (MD,NNZERO,AA,JA,IA,X,Y)
IMPLICIT NONE
INTEGER                           :: I,MD,NNZERO,K
REAL(8),DIMENSION(NNZERO)         :: AA
REAL(8),DIMENSION(MD)             :: Y,X
INTEGER,DIMENSION(NNZERO)         :: JA
INTEGER,DIMENSION(MD+1)           :: IA
!------------------------------------------------------------------------------
DO I=1,MD
   Y(I) = 0.D0
END DO
DO I=1,MD
DO K=IA(I),IA(I+1)-1
   Y(JA(K))=Y(JA(K))+X(I)*AA(K)
END DO
END DO
!------------------------------------------------------------------------------
END SUBROUTINE
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                         ASOLVE
!______________________________________________________________________________
SUBROUTINE ASOLVE (TL,MD,SA,B,X,ITRNSP)
IMPLICIT NONE
INTEGER                            :: ITRNSP,IJA,I,TL,MD
REAL(8),DIMENSION(MD)              :: X,B
REAL(8),DIMENSION(TL)              :: SA
!------------------------------------------------------------------------------
DO I=1,MD
   X(I)=B(I)/SA(I)
end DO
!------------------------------------------------------------------------------
END SUBROUTINE ASOLVE
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                         ATIMES
!______________________________________________________________________________
SUBROUTINE ATIMES (TL,MD,SA,IJA,X,R,ITRNSP)
IMPLICIT NONE
INTEGER                            :: ITRNSP,TL,MD
REAL(8),DIMENSION(MD)              :: X,R
REAL(8),DIMENSION(TL)              :: SA
INTEGER,DIMENSION(TL)              :: IJA
!------------------------------------------------------------------------------
IF (ITRNSP.EQ.0) THEN
      CALL SPRSAX (TL,MD,SA,IJA,X,R)
   ELSE
      CALL SPRSTX (TL,MD,SA,IJA,X,R)
END IF
!------------------------------------------------------------------------------
END SUBROUTINE ATIMES
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                      SNRM_FUNCTION
!______________________________________________________________________________
FUNCTION SNRM (MD,SX,ITOL)
IMPLICIT NONE
INTEGER                     :: MD,ITOL,I,ISAMAX
REAL(8),DIMENSION(MD)        :: SX
REAL(8)                     :: SNRM
!------------------------------------------------------------------------------
IF (ITOL.LE.3) THEN
   SNRM=0.D0
   DO I=1,MD
      SNRM=SNRM+SX(I)**2
   END DO
      SNRM=SQRT(SNRM)/MD
      ELSE
      ISAMAX=1
       DO I=1,MD
          IF (ABS(SX(I)).GT.ABS(SX(ISAMAX))) ISAMAX=I
       END DO
          SNRM=ABS(SX(ISAMAX))
       END IF
!------------------------------------------------------------------------------
END FUNCTION SNRM
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                         SPRSAX
!______________________________________________________________________________
SUBROUTINE SPRSAX(TL,MD,SA,IJA,X,B)
IMPLICIT NONE
INTEGER                       :: I,K,TL,MD
REAL(8),DIMENSION(MD)         :: B,X
REAL(8),DIMENSION(TL)         :: SA
INTEGER,DIMENSION(TL)         :: IJA
!------------------------------------------------------------------------------
IF (IJA(1).NE.MD+2) THEN
PRINT*, 'mismatched vector and matrix in sprsax'
ELSE
DO I=1,MD
   B(I)=SA(I)*X(I)
   DO K=IJA(I),IJA(I+1)-1
      B(I)=B(I)+SA(K)*X(IJA(K))
   END DO
END DO
END IF
!------------------------------------------------------------------------------
END SUBROUTINE SPRSAX
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                          SPRSTX
!______________________________________________________________________________
SUBROUTINE SPRSTX (TL,MD,SA,IJA,X,B)
IMPLICIT NONE
INTEGER                     :: I,J,K,TL,MD
REAL(8),DIMENSION(MD)       :: B,X
INTEGER,DIMENSION(TL)       :: IJA
REAL(8),DIMENSION(TL)       :: SA
!------------------------------------------------------------------------------
IF (IJA(1).NE.MD+2) THEN
PRINT*, 'mismatched vector and matrix in sprstx'
ELSE
   DO I=1,MD
      B(I)=SA(I)*X(I)
   END DO
   DO I=1,MD
      DO K=IJA(I),IJA(I+1)-1
         J=IJA(K)
         B(J)=B(J)+SA(K)*X(I)
      END DO
   END DO
END IF
!------------------------------------------------------------------------------
END SUBROUTINE SPRSTX
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                       TOTAL_LENGTH
!______________________________________________________________________________
SUBROUTINE TOTAL_LENGTH (MD,IA,TL)
IMPLICIT NONE
INTEGER                           :: TL,MD
INTEGER,DIMENSION(MD+1)           :: IA
!------------------------------------------------------------------------------
TL=IA(MD+1)
!------------------------------------------------------------------------------
END SUBROUTINE TOTAL_LENGTH
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                        CSRMSR SUBROUTINE
!          Compressed Sparse Row to Modified-Sparse Row
!             Sparse row with separate main diagonal
!------------------------------------------------------------------------------
! converts a general sparse matrix A, JA, IA into
! a compressed matrix using a separated diagonal (referred to as
! the bell-labs format as it is used by bell labs semi conductor
! group.
!
!
! In case AO, JAO, are different from AA, JA.
!------------------------------------------------------------------------------
!
! on entry :
!------------------------------------------------------------------------------
! AA, JA, IA = matrix in csr format. note that the
!	     algorithm is in place: AO, JAO can be the same
!            as AA, JA, in which case it will be overwritten on it
!            upon return.
!
! on return :
!------------------------------------------------------------------------------
!
! AO, JAO  = sparse matrix in modified sparse row storage format:
!	   +  AO(1:N) contains the diagonal of the matrix.
!	   +  AO(N+2:TL) contains the nondiagonal elements of the
!             matrix, stored rowwise.
!	   +  JAO(N+2:TL) : their column indices
!	   +  JAO(1:N+1) contains the pointer array for the
!             nondiagonal elements in AO(N+1:TL) and JAO(N+2:TL).
!         i.e., for I .LE. N+1 JAO(I) points to beginning of row I
!	      in arrays AO, JAO.
!	       here TL = number of nonzero elements+1
! work arrays:
!------------------------------------------------------------------------------
! WK	= real work array of length N
! IWK   = integer work array of length N+1
!
!------------------------------------------------------------------------------
!        Algorithm is in place when
!
!        call csrmsr (TL,NNZERO,AA,JA,IA,AO,JAO)
!        (in which  AO, JAO, are different from AA, JA)
!
!------------------------------------------------------------------------------
! coded by Y. Saad Sep. 1989. Rechecked Feb 27, 1990.
!------------------------------------------------------------------------------
! This code has been rewritten in FORTRAN 90 by I. Farahbakhsh in
! october 2009.
! The arguments are reduced and the array allocation concept is added.
!______________________________________________________________________________
SUBROUTINE CSRMSR  (TL,MD,NNZERO,AA,JA,IA,AO,JAO)
IMPLICIT NONE
INTEGER                               :: I,J,K,NNZERO,II,TL,MD
REAL(8),DIMENSION(NNZERO)             :: AA
REAL(8),DIMENSION(TL)                 :: AO
REAL(8),DIMENSION(MD)                 :: WK,Z,V
INTEGER,DIMENSION(MD+1)               :: IA,IWK
INTEGER,DIMENSION(NNZERO)             :: JA
INTEGER,DIMENSION(TL)                 :: JAO
!------------------------------------------------------------------------------
DO  I=1,MD
WK(I) = 0.D0
IWK(I+1) = IA(I+1)-IA(I)
DO  K=IA(I),IA(I+1)-1
IF (JA(K).EQ.I) THEN
WK(I)= AA(K)
IWK(I+1)=IWK(I+1)-1
END IF
END DO
END DO
!------------------------------------------------------------------------------
! Copy backwards (to avoid collisions)
!------------------------------------------------------------------------------
DO II=MD,1,-1
   DO K=IA(II+1)-1,IA(II),-1
      J=JA(K)
      IF (J.NE.II) THEN
         AO(TL) = AA(K)
         JAO(TL) = J
         TL = TL-1
      END IF
   END DO
END DO
!------------------------------------------------------------------------------
! Compute pointer values and copy wk(N)
!------------------------------------------------------------------------------
JAO(1)=MD+2
      DO  I=1,MD
          AO(I)=WK(I)
          JAO(I+1)=JAO(I)+IWK(I+1)
      END DO
!------------------------------------------------------------------------------
END SUBROUTINE CSRMSR
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                     SUCCESSIVE_SOLUTION SUBROUTINE
!______________________________________________________________________________
!                     SUCCESSIVE_SOLUTION SUBROUTINE
!______________________________________________________________________________
SUBROUTINE SUCCESSIVE_SOL (METHOD,ERR,ITER)
IMPLICIT NONE
INTEGER                             :: ITER
CHARACTER(*)                        :: METHOD
CHARACTER*30                        :: FNAME
REAL(8)                             :: ERR
!------------------------------------------------------------------------------
FNAME="/results/CONV-"//METHOD//".PLT"
!------------------------------------------------------------------------------
OPEN(3,FILE=FNAME)
WRITE(3,*) ITER,ERR
!-----------------------------------------------------------------------------
END SUBROUTINE SUCCESSIVE_SOL
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                         CG_SUBROUTINE
!______________________________________________________________________________
SUBROUTINE CG (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,X_OLD,R_OLD,X)
USE PARAM, ONLY                     : TOL,SAVE_ITERATION_NORM
IMPLICIT NONE
INTEGER                            :: K,ITER,TL,MD,NNZERO,IDIAG,IERR
REAL(8),DIMENSION(MD)              :: X_OLD,X,P_OLD,P,R_OLD,R,AP
INTEGER,DIMENSION(MD+1)            :: IA_AUX,IAR
INTEGER,DIMENSION(2*MD-1)          :: IND
REAL(8)                            :: NORM,S,ALPHA,BETA,M,MM,TIME_BEGIN,TIME_END
REAL(8),DIMENSION(NNZERO)          :: AA_AUX,AAR
INTEGER,DIMENSION(NNZERO)          :: JA_AUX,JAR
REAL(8),ALLOCATABLE                :: SA(:),DIAG(:,:),COEF(:,:)
INTEGER,ALLOCATABLE                :: IJA(:),IOFF(:),JCOEF(:,:)
!------------------------------------------------------------------------------
CALL TOTAL_LENGTH          (MD,IA_AUX,TL)
ALLOCATE                   (SA(TL),IJA(TL))
CALL CSRMSR                (TL,MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,SA,IJA)
!CALL TOTAL_LENGTH          (MD,IA_AUX,TL)
CALL INFDIA                (MD,NNZERO,JA_AUX,IA_AUX,IND,IDIAG)
ALLOCATE                   (COEF(MD,IDIAG),JCOEF(MD,IDIAG))
ALLOCATE                   (IOFF(IDIAG),DIAG(MD,IDIAG))
CALL CSRELL                (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,MD,COEF,JCOEF,MD,IDIAG,IERR)
CALL CSRDIA                (MD,NNZERO,IDIAG,10,AA_AUX,JA_AUX,IA_AUX,MD,DIAG,IOFF,AAR,JAR,IAR,IND)
!------------------------------------------------------------------------------
P_OLD=R_OLD
NORM=1.D0
ITER=0
CALL CPU_TIME (TIME_BEGIN)
DO WHILE (NORM.GT.TOL)
ITER=ITER+1
!PRINT*,ITER
S=0.D0
M=0.D0
MM=0.D0
!------------------------------------------------------------------------------
! Stores the matrix in compressed spars row and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL CSR_MAT_V_PRODUCT     (MD,NNZERO,AA_AUX,JA_AUX,IA_AUX,P_OLD,AP)
!------------------------------------------------------------------------------
! Stores the matrix in modified spars row and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL SPRSAX (TL,MD,SA,IJA,P_OLD,AP)
!------------------------------------------------------------------------------
! Stores the matrix in Ellpack/Itpack format and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL AMUXE (MD,P_OLD,AP,MD,IDIAG,COEF,JCOEF)
!------------------------------------------------------------------------------
! Stores the matrix in Diagonal format and multiplies by
! vector
!------------------------------------------------------------------------------
CALL AMUXD (MD,P_OLD,AP,DIAG,MD,IDIAG,IOFF)
!------------------------------------------------------------------------------
DO K=1,MD
M=M+R_OLD(K)**2
MM=MM+AP(K)*P_OLD(K)
END DO
ALPHA=M/MM
X=X_OLD+ALPHA*P_OLD
R=R_OLD-ALPHA*AP
DO K=1,MD
S=S+R(K)**2
END DO
NORM=SQRT(S)/MD
BETA=S/M
P=R+BETA*P_OLD
P_OLD=P
X_OLD=X
R_OLD=R
!PRINT*,NORM
!CALL CPU_TIME (TIME_END)
!CALL CPU_TIME_WRITE (TIME_BEGIN,TIME_END,ITER)
IF (SAVE_ITERATION_NORM) CALL SUCCESSIVE_SOL ("CG",NORM,ITER)
END DO
CALL CPU_TIME (TIME_END)
!PRINT*, (TIME_END-TIME_BEGIN)/1.5625D-2
!------------------------------------------------------------------------------
END SUBROUTINE CG
!______________________________________________________________________________
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE MATRIX_PATTERN (MD,NN,AA,JA,IA,DNS)
IMPLICIT NONE
INTEGER                       :: K,I,J,IERR,NN,MD
REAL(8),DIMENSION(NN)         :: AA
REAL(8),DIMENSION(MD,MD)      :: DNS
INTEGER,DIMENSION(MD+1)       :: IA
INTEGER,DIMENSION(NN)         :: JA,ROW,COL
!------------------------------------------------------------------------------
! Compressed Sparse Row to Dense
!------------------------------------------------------------------------------
!!This subroutine has been written by Y.Saad on Sep 1989
!and has been changed in FORTRAN 90 via Iman Farahbakhsh on Dec 2009.
!------------------------------------------------------------------------------
! converts a row-stored sparse matrix into a densely stored one
!
! On entry:
!----------
!
! MD     = matrix dimension
! AA,
! JA,
! IA    = input matrix in compressed sparse row format.
!         (AA=value array, JA=column array, IA=pointer array)
! DNS   = array where to store dense matrix
!
! on return:
!-----------
! DNS   = the sparse matrix AA, JA, IA has been stored in DNS(MD,MD)
!
! IERR  = integer error indicator.
!         IERR .eq. 0  means normal return
!         IERR .eq. I  means that the code has stopped when processing
!         row number I, because it found a column number .gt. MD.
!------------------------------------------------------------------------------
OPEN(UNIT=1001,FILE='PATTERN.DAT',STATUS='UNKNOWN')
!------------------------------------------------------------------------------
IERR=0
DO  I=1,MD
    DO  J=1,MD
        DNS(I,J)=0.D0
    END DO
END DO
!------------------------------------------------------------------------------
DO  I=1,MD
    DO  K=IA(I),IA(I+1)-1
        J=JA(K)
        IF (J.GT.MD) THEN
           IERR=I
        END IF
	    DNS(I,J)=AA(K)
    END DO
END DO
K=0
DO I=1,MD
DO J=1,MD
IF (DNS(I,J).NE.0.D0) THEN
K=K+1
ROW(K)=I
COL(K)=J
END IF
END DO
END DO
!------------------------------------------------------------------------------
WRITE (1001,*) 'VARIABLES=I,J'
WRITE (1001,*) 'ZONE'
WRITE (1001,*) 'F=POINT'
DO I=1,NN
WRITE (1001,*) ROW(I),COL(NN+1-I)
END DO
!------------------------------------------------------------------------------
END SUBROUTINE MATRIX_PATTERN
!******************************************************************************
!******************************************************************************
!                     MATRIX_PATTERN_2 SUBROUTINE
!******************************************************************************
!******************************************************************************
!  Compressed Sparse Row      to      Coordinate
!------------------------------------------------------------------------------
! This subroutine was originally written by Y. Saad and was revised by
! Iman Farahbakhsh in Feb 2010.
!------------------------------------------------------------------------------
! converts a matrix that is stored in coordinate format
!  A, IR, JC into a row general sparse AO, JAO, IAO format.
!
! on entry:
!---------
! NROW	= dimension of the matrix.
! JOB   = integer serving as a JOB indicator.
!         if JOB = 1 fill in only the array IR, ignore JC, and AO.
!         if JOB = 2 fill in IR, and JC but not AO
!         if JOB = 3 fill in everything.
!         The reason why these options are provided is that on return
!         AO and JC are the same as A, JA. So when JOB = 3, A and JA are
!         simply copied into AO, JC.  When JOB=2, only JC and IR are
!         returned. With JOB=1 only the array IR is returned. Moreover,
!         the algorithm is in place:
!	      CALL CSRCOO (NROW,1,A,JA,IA,NNZ,A,IA,JA)
!         will write the output matrix in coordinate format on A, JA,IA.
!
! A,
! JA,
! IA    = matrix in compressed sparse row format.
!
! on return:
!-----------
! AO, IR, JC = matrix in coordinate format.
!
! NNZ        = number of nonzero elements in matrix.
!
! NOTES: 1)This routine is PARTIALLY in place: CSRCOO can be called with
!         AO being the same array as as A, and JC the same array as JA.
!         but IR CANNOT be the same as IA.
!         2) note the order in the output arrays,
!------------------------------------------------------------------------------
! Last change: Feb 4 2010 By Iman Farahbakhsh
! This subroutine revised in FORTRAN 90 instead of FORTRAN 77 and enable us to
! write the matrix pattern in Tecplot software.
!______________________________________________________________________________
SUBROUTINE MATRIX_PATTERN_2 (NROW,JOB,A,JA,IA,NNZ,AO,IR,JC)
!------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER                     :: NROW,JOB,NNZ,I,K,K1,K2
REAL(8),DIMENSION(NNZ)      :: A,AO
INTEGER,DIMENSION(NNZ)      :: JA,JC,IR
INTEGER,DIMENSION(NROW+1)   :: IA

!------------------------------------------------------------------------------
OPEN(UNIT=1001,FILE='PATTERN.DAT',STATUS='UNKNOWN')
!------------------------------------------------------------------------------
  GO TO (3,2,1) JOB
1 DO  K=1,NNZ
         AO(K)=A(K)
  END DO
2 DO  K=1,NNZ
         JC(K)=JA(K)
  ENDDO
!
!     copy backward to allow for in-place processing.
!
3 DO I=NROW,1,-1
         K1=IA(I+1)-1
         K2=IA(I)
         DO K=K1,K2,-1
            IR(K)=I
         ENDDO
  ENDDO
!------------------------------------------------------------------------------
WRITE (1001,*) 'VARIABLES=I,J'
WRITE (1001,*) 'ZONE'
WRITE (1001,*) 'F=POINT'
DO I=1,NNZ
WRITE (1001,*) IR(I),JC(NNZ+1-I)
END DO
!------------------------------------------------------------------------------
END SUBROUTINE MATRIX_PATTERN_2
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE ILU0 (N,NNZERO,A,JA,IA,LUVAL,UPTR,IW,ICODE)
IMPLICIT NONE
INTEGER                    :: N,NNZRO,I,J,K,J1,J2,ICODE,JROW,NNZERO,JJ,JW
INTEGER,DIMENSION(NNZERO)  :: JA
INTEGER,DIMENSION(N+1)     :: IA
INTEGER,DIMENSION(N)       :: IW,UPTR
REAL(8),DIMENSION(NNZERO)  :: A
REAL(8),DIMENSION(NNZERO)  :: LUVAL
REAL(8)                    :: T
!------------------------------------------------------------------------------
DO I=1,IA(N+1)-1
   LUVAL(I)=A(I)
END DO
DO I=1,N
   IW(I)=0
END DO
!------------------------------------------------------------------------------
DO K=1,N
   J1=IA(K)
   J2=IA(K+1)-1
   DO J=J1,J2
      IW(JA(J))=J
   END DO
      J=J1
150   JROW=JA(J)
!------------------------------------------------------------------------------
   IF (JROW.GE.K) GO TO 200
!------------------------------------------------------------------------------
   T=LUVAL(J)*LUVAL(UPTR(JROW))
   LUVAL(J)=T
!------------------------------------------------------------------------------
   DO JJ=UPTR(JROW)+1,IA(JROW+1)-1
      JW=IW(JA(JJ))
	  IF (JW.NE.0) LUVAL(JW)=LUVAL(JW)-T*LUVAL(JJ)
   END DO
   J=J+1
   IF (J.LE.J2) GO TO 150
!------------------------------------------------------------------------------
200   UPTR(K)=J
   IF (JROW.NE.K.OR.LUVAL(J).EQ.0.D0) GO TO 600
   LUVAL(J)=1.D0/LUVAL(J)
!------------------------------------------------------------------------------
   DO I=J1,J2
      IW(JA(I))=0
   END DO
END DO
!------------------------------------------------------------------------------
ICODE=0
RETURN
!------------------------------------------------------------------------------
600  ICODE=K
RETURN
!------------------------------------------------------------------------------
END SUBROUTINE ILU0
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE LUSOL (N,NNZERO,RHS,SOL,LUVAL,LUCOL,LUPTR,UPTR)
IMPLICIT NONE
INTEGER                       :: I,K,N,NNZERO
REAL(8),DIMENSION(NNZERO)     :: LUVAL
REAL(8),DIMENSION(N)          :: SOL,RHS
INTEGER,DIMENSION(N)          :: UPTR
INTEGER,DIMENSION(N+1)        :: LUPTR
INTEGER,DIMENSION(NNZERO)     :: LUCOL
!------------------------------------------------------------------------------
DO I=1,N
   SOL(I)=RHS(I)
   DO K=LUPTR(I),UPTR(I)-1
      SOL(I)=SOL(I)-LUVAL(K)*SOL(LUCOL(K))
   END DO
END DO
!------------------------------------------------------------------------------
DO I=N,1,-1
   DO K=UPTR(I)+1,LUPTR(I+1)-1
      SOL(I)=SOL(I)-LUVAL(K)*SOL(LUCOL(K))
   END DO
   SOL(I)=LUVAL(UPTR(I))*SOL(I)
END DO
!------------------------------------------------------------------------------
END SUBROUTINE LUSOL
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE ILU0_PCG (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,X_OLD,R_OLD,X)
USE PARAM, ONLY                     : TOL,SAVE_ITERATION_NORM
IMPLICIT NONE
INTEGER                            :: K,ITER,TL,N,NNZERO,ICODE,IDIAG,IERR
REAL(8),DIMENSION(N)               :: X_OLD,X,P_OLD,P,R_OLD,R,AP,Z_OLD,Z
INTEGER,DIMENSION(N+1)             :: IA_AUX,IAR
INTEGER,DIMENSION(2*N-1)           :: IND
REAL(8)                            :: NORM,S,ALPHA,BETA,M,MM
REAL(8),DIMENSION(NNZERO)          :: AA_AUX,AAR,LUVAL
INTEGER,DIMENSION(NNZERO)          :: JA_AUX,JAR
REAL(8),ALLOCATABLE                :: SA(:),COEF(:,:),DIAG(:,:)
INTEGER,ALLOCATABLE                :: IJA(:),JCOEF(:,:),IOFF(:)
INTEGER,DIMENSION(N)               :: UPTR,IW
REAL(8)                            :: SN,TIME_BEGIN,TIME_END
!------------------------------------------------------------------------------
CALL TOTAL_LENGTH          (N,IA_AUX,TL)
ALLOCATE                   (SA(TL),IJA(TL))
CALL CSRMSR                (TL,N,NNZERO,AA_AUX,JA_AUX,IA_AUX,SA,IJA)
!CALL TOTAL_LENGTH          (N,IA_AUX,TL)
CALL INFDIA                (N,NNZERO,JA_AUX,IA_AUX,IND,IDIAG)
ALLOCATE                   (COEF(N,IDIAG),JCOEF(N,IDIAG))
ALLOCATE                   (IOFF(IDIAG),DIAG(N,IDIAG))
CALL CSRELL                (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,N,COEF,JCOEF,N,IDIAG,IERR)
CALL CSRDIA                (N,NNZERO,IDIAG,10,AA_AUX,JA_AUX,IA_AUX,N,DIAG,IOFF,AAR,JAR,IAR,IND)
!------------------------------------------------------------------------------
CALL ILU0   (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,LUVAL,UPTR,IW,ICODE)
CALL LUSOL  (N,NNZERO,R_OLD,Z_OLD,LUVAL,JA_AUX,IA_AUX,UPTR)
!------------------------------------------------------------------------------
P_OLD=Z_OLD
NORM=1.D0
ITER=0
!CALL CPU_TIME (TIME_BEGIN)
DO WHILE (NORM.GT.TOL)
ITER=ITER+1
!PRINT*,ITER
SN=0.D0
S=0.D0
M=0.D0
MM=0.D0
!------------------------------------------------------------------------------
! Stores the matrix in compressed spars row and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL CSR_MAT_V_PRODUCT     (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,P_OLD,AP)
!------------------------------------------------------------------------------
! Stores the matrix in modified spars row and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL SPRSAX (TL,N,SA,IJA,P_OLD,AP)
!------------------------------------------------------------------------------
! Stores the matrix in Ellpack/Itpack format and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL AMUXE (N,P_OLD,AP,N,IDIAG,COEF,JCOEF)
!------------------------------------------------------------------------------
! Stores the matrix in Diagonal format and multiplies by
! vector
!------------------------------------------------------------------------------
CALL AMUXD (N,P_OLD,AP,DIAG,N,IDIAG,IOFF)
!------------------------------------------------------------------------------
DO K=1,N
M=M+R_OLD(K)*Z_OLD(K)
MM=MM+AP(K)*P_OLD(K)
END DO
ALPHA=M/MM
X=X_OLD+ALPHA*P_OLD
R=R_OLD-ALPHA*AP
!------------------------------------------------------------------------------
CALL LUSOL  (N,NNZERO,R,Z,LUVAL,JA_AUX,IA_AUX,UPTR)
!------------------------------------------------------------------------------
DO K=1,N
S=S+R(K)*Z(K)
END DO
BETA=S/M
P=Z+BETA*P_OLD
DO K=1,N
SN=SN+R(K)**2
END DO
NORM=SQRT(SN)/N
P_OLD=P
X_OLD=X
R_OLD=R
Z_OLD=Z
!PRINT*,NORM
!CALL CPU_TIME (TIME_END)
!CALL CPU_TIME_WRITE (TIME_BEGIN,TIME_END,ITER)
IF (SAVE_ITERATION_NORM) CALL SUCCESSIVE_SOL ("ILU0_PCG",NORM,ITER)
END DO
!------------------------------------------------------------------------------
END SUBROUTINE ILU0_PCG
!______________________________________________________________________________
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                     CPU_TIME_WRITE SUBROUTINE
!______________________________________________________________________________
SUBROUTINE  CPU_TIME_WRITE (TIME_BEGIN,TIME_END,ITER)
IMPLICIT NONE
INTEGER                             :: ITER,COUNTER
CHARACTER*10                        :: EXT
CHARACTER*5                         :: FN1
CHARACTER*25                        :: FNAME
REAL(8)                             :: TIME_BEGIN,TIME_END
!------------------------------------------------------------------------------
FN1='CPU_T'
COUNTER=0
!------------------------------------------------------------------------------
COUNTER=COUNTER+10
WRITE(EXT,'(I7)') ITER
FNAME=FN1//EXT//'.DAT'
!------------------------------------------------------------------------------
OPEN(COUNTER,FILE=FNAME,POSITION='REWIND')
!WRITE(COUNTER,*)'VARIABLES= "ITER","TIME"'
!WRITE(COUNTER,*)'ZONE,F=POINT'
WRITE(COUNTER,*) ITER,TIME_END-TIME_BEGIN
CLOSE(COUNTER)
!------------------------------------------------------------------------------
!WRITE(*,*) '========================'
!WRITE(*,*) 'PRINTING ON ',FNAME
!WRITE(*,*) '========================'
!------------------------------------------------------------------------------
END SUBROUTINE CPU_TIME_WRITE
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE INFDIA (N,NN,JA,IA,IND,IDIAG)
IMPLICIT NONE
INTEGER                          :: N,NN,N2,I,K,J,IDIAG
INTEGER,DIMENSION(N+1)           :: IA
INTEGER,DIMENSION(2*N-1)         :: IND
INTEGER,DIMENSION(NN)            :: JA
!------------------------------------------------------------------------------
!     obtains information on the diagonals of A.
!------------------------------------------------------------------------------
! this subroutine finds the lengths of each of the 2*N-1 diagonals of A
! it also outputs the number of nonzero diagonals found.
!------------------------------------------------------------------------------
! on entry:
!----------
! N	= dimension of the matrix A.
!
! A,    ..... not needed here.
! JA,
! IA    = matrix stored in csr format
!
! on return:
!-----------
!
! IDIAG = integer. number of nonzero diagonals found.
!
! IND   = integer array of length at least 2*N-1. The K-th entry in
!         ind contains the number of nonzero elements in the diagonal
!         number K, the numbering beeing from the lowermost diagonal
!         (bottom-left). In other words ind(K) = length of diagonal
!         whose offset wrt the main diagonal is = - N + K.
!------------------------------------------------------------------------------
!           Y. Saad, Sep. 21 1989
!           with some revisions by I. Farahbakhsh, Dec 2009
!------------------------------------------------------------------------------
      N2=N+N-1
      DO 1 I=1,N2
         IND(I)=0
 1    CONTINUE
      DO 3 I=1,N
         DO 2 K=IA(I),IA(I+1)-1
            J=JA(K)
            IND(N+J-I)=IND(N+J-I)+1
 2       CONTINUE
 3    CONTINUE
!     count the nonzero ones.
      IDIAG=0
      DO 41 K=1,N2
         IF (IND(K).NE.0) IDIAG=IDIAG+1
 41   CONTINUE
      RETURN
! done
!------end-of-INFDIA-----------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE INFDIA
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE CSRELL (NROW,NN,A,JA,IA,MAXCOL,COEF,JCOEF,NCOEF,NDIAG,IERR)
IMPLICIT NONE
INTEGER                         :: NROW,NN,NCOEF,NDIAG,IERR,I,J,K,K1,K2,MAXCOL
INTEGER,DIMENSION(NROW+1)       :: IA
INTEGER,DIMENSION(NN)           :: JA
INTEGER,DIMENSION(NCOEF,NDIAG)  :: JCOEF
REAL(8),DIMENSION(NN)           :: A
REAL(8),DIMENSION(NCOEF,NDIAG)  :: COEF
!------------------------------------------------------------------------------
! Compressed Sparse Row	    to    Ellpack - Itpack format
!------------------------------------------------------------------------------
! this subroutine converts  matrix stored in the general A, JA, IA
! format into the COEF, JCOEF itpack format.
!
!------------------------------------------------------------------------------
! on entry:
!----------
! NROW 	  = row dimension of the matrix A.
!
! A,
! IA,
! JA     = input matrix in compressed sparse row format.
!
! NCOEF  = first dimension of arrays COEF, and JCOEF.
!
! MAXCOL = integer equal to the number of columns available in COEF.
!
! on RETURN:
!----------
! COEF	= real array containing the values of the matrix A in
!         itpack-ellpack format.
! JCOEF = integer array containing the column indices of COEF(i,j)
!         in A.
! NDIAG = number of active 'diagonals' found.
!
! IERR 	= error message. 0 = correct return. If ierr .ne. 0 on
!	  return this means that the number of diagonals found
!         (ndiag) exceeds maxcol.
!------------------------------------------------------------------------------
!           Y. Saad, Sep. 21 1989
!           with some revisions by I. Farahbakhsh, Dec 2009
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! first determine the length of each row of lower-part-of(A)
      IERR=0
      NDIAG=0
      DO 3 I=1,NROW
         K=IA(I+1)-IA(I)
         NDIAG=MAX0(NDIAG,K)
 3    CONTINUE
!----- check whether sufficient columns are available. ------------------------
      IF (NDIAG.GT.MAXCOL) THEN
         IERR=1
         RETURN
      ENDIF
!
! fill COEF with zero elements and jcoef with row numbers.---------------------
!
      DO 4 J=1,NDIAG
         DO 41 I=1,NROW
            COEF(I,J)=0.D0
            JCOEF(I,J)=I
 41      CONTINUE
 4    CONTINUE
!
!------- copy elements row by row.---------------------------------------------
!
      DO 6 I=1,NROW
         K1=IA(I)
         K2=IA(I+1)-1
         DO 5 K=K1,K2
            COEF(I,K-K1+1)=A(K)
            JCOEF(I,K-K1+1)=JA(K)
 5       CONTINUE
 6    CONTINUE
      RETURN
!----end-of-CSRELL-------------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE CSRELL
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE AMUXE (N,X,Y,NA,NCOL,A,JA)
IMPLICIT NONE
INTEGER                              :: N,NA,NCOL,I,J
REAL(8),DIMENSION(N)                 :: X,Y
REAL(8),DIMENSION(NA,NCOL)           :: A
INTEGER,DIMENSION(NA,NCOL)           :: JA
!------------------------------------------------------------------------------
!        A times a vector in Ellpack Itpack format (ELL)
!------------------------------------------------------------------------------
! multiplies a matrix by a vector when the original matrix is stored
! in the ellpack-itpack sparse format.
!------------------------------------------------------------------------------
!
! on entry:
!----------
! N     = row dimension of A
! X     = real array of length equal to the column dimension of
!         the A matrix.
! NA    = integer. The first dimension of arrays A and JA
!         as declared by the calling program.
! NCOL  = integer. The number of active columns in array A.
!         (i.e., the number of generalized diagonals in matrix.)
! A,JA  = the real and integer arrays of the itpack format
!         (A(I,K),K=1,NCOL contains the elements of row I in matrix
!          JA(I,K),K=1,NCOL contains their column numbers)
!
! on return:
!-----------
! Y     = real array of length N, containing the product Y=A*X
!------------------------------------------------------------------------------
!           Y. Saad, Sep. 21 1989
!           with some revision by I. Farahbakhsh, Dec 2009
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      DO 1 I=1,N
         Y(I)=0.D0
 1    CONTINUE
      DO 10 J=1,NCOL
         DO 25 I=1,N
            Y(I)=Y(I)+A(I,J)*X(JA(I,J))
 25      CONTINUE
 10   CONTINUE
!
      RETURN
!--------end-of-AMUXE----------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE AMUXE
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE CSRDIA (N,NN,IDIAG,JOB,A,JA,IA,NDIAG,DIAG,IOFF,AO,JAO,IAO,IND)
IMPLICIT NONE
INTEGER                             :: N,NN,IDIAG,NDIAG,JOB1,JOB2,JOB,N2,IDUM,II,JMAX,K,I,J,L,KO
REAL(8),DIMENSION(NDIAG,IDIAG)      :: DIAG
REAL(8),DIMENSION(NN)               :: A,AO
INTEGER,DIMENSION(NN)               :: JA,JAO
INTEGER,DIMENSION(N+1)              :: IA,IAO
INTEGER,DIMENSION(2*N-1)            :: IND
INTEGER,DIMENSION(IDIAG)            :: IOFF
!------------------------------------------------------------------------------
! Compressed sparse row     to    diagonal format
!------------------------------------------------------------------------------
! this subroutine extracts  IDIAG diagonals  from the  input matrix A,
! A, IA, and puts the rest of  the matrix  in the  output matrix AO,
! JAO, IAO.  The diagonals to be extracted depend  on the  value of JOB
! (see below for details.)  In  the first  case, the  diagonals to be
! extracted are simply identified by  their offsets  provided in IOFF
! by the caller.  In the second case, the  code internally determines
! the IDIAG most significant diagonals, i.e., those  diagonals of the
! matrix which  have  the  largest  number  of  nonzero elements, and
! extracts them.
!------------------------------------------------------------------------------
! on entry:
!----------
! N	= dimension of the matrix A.
! IDIAG = integer equal to the number of diagonals to be extracted.
!         Note: on return IDIAG may be modified.
! A, JA,
!    IA = matrix stored in A, JA, IA, format
! JOB	= integer. serves as A JOB indicator.  JOB is better thought
!         of as A two-digit number JOB=xy. If the first (x) digit
!         is one on entry then the diagonals to be extracted are
!         internally determined. In this case csrdia exctracts the
!         IDIAG most important diagonals, i.e. those having the largest
!         number on nonzero elements. If the first digit is zero
!         then csrdia assumes that IOFF(*) contains the offsets
!         of the diagonals to be extracted. there is no verification
!         that IOFF(*) contains valid entries.
!         The second (y) digit of JOB determines whether or not
!         the remainder of the matrix is to be written on AO,JAO,IAO.
!         If it is zero  then AO, JAO, IAO is not filled, i.e.,
!         the diagonals are found  and put in array DIAG and the rest is
!         is discarded. if it is one, AO, JAO, IAO contains matrix
!         of the remaining elements.
!         Thus:
!         JOB= 0 means do not select diagonals internally (pick those
!                defined by IOFF) and do not fill AO,JAO,IAO
!         JOB= 1 means do not select diagonals internally
!                      and fill AO,JAO,IAO
!         JOB=10 means  select diagonals internally
!                      and do not fill AO,JAO,IAO
!         JOB=11 means select diagonals internally
!                      and fill AO,JAO,IAO
!
! NDIAG = integer equal to the first dimension of array DIAG.
!
! on return:
!-----------
!
! IDIAG = number of diagonals found. This may be smaller than its value
!         on entry.
! DIAG  = real array of size (NDIAG x IDIAG) containing the diagonals
!         of A on return
!
! IOFF  = integer array of length IDIAG, containing the offsets of the
!   	  diagonals to be extracted.
! AO, JAO
!  IAO  = remainder of the matrix in A, JA, IA format.
! work arrays:
!------------
! IND   = integer array of length 2*N-1 used as integer work space.
!         needed only when JOB.ge.10 i.e., in case the diagonals are to
!         be selected internally.
!
! Notes:
!-------
!    1) The algorithm is in place: AO, JAO, IAO can be overwritten on
!       A, JA, IA if desired
!    2) When the code is required to select the diagonals (JOB .ge. 10)
!       the selection of the diagonals is done from left to right
!       as A result if several diagonals have the same weight (number
!       of nonzero elemnts) the leftmost one is selected first.
!------------------------------------------------------------------------------
      JOB1 = JOB/10
      JOB2 = JOB-JOB1*10
      IF (JOB1 .EQ. 0) GOTO 50
      N2 = N+N-1
      CALL INFDIA(N,NN,JA,IA,IND,IDUM)
!----------- determine diagonals to  accept.-----------------------------------
!------------------------------------------------------------------------------
      II=0
 4    II=II+1
      JMAX=0
      DO 41 K=1,N2
         J=IND(K)
         IF (J.LE.JMAX) GOTO 41
         I=K
         JMAX=J
 41   CONTINUE
      IF (JMAX.LE.0) THEN
         II=II-1
         GOTO 42
      ENDIF
      IOFF(II)=I-N
      IND(I)=-JMAX
      IF (II.LT.IDIAG) GOTO 4
 42   IDIAG=II
!---------------- initialize diago to zero ------------------------------------
 50   CONTINUE
      DO 55 J=1,IDIAG
         DO 54 I=1,N
            DIAG(I,J)=0.D0
 54      CONTINUE
 55   CONTINUE
!------------------------------------------------------------------------------
      KO=1
!------------------------------------------------------------------------------
! extract diagonals and accumulate remaining matrix.
!------------------------------------------------------------------------------
      DO 6 I=1,N
         DO 51 K=IA(I),IA(I+1)-1
            J=JA(K)
            DO 52 L=1,IDIAG
               IF (J-I.NE.IOFF(L)) GOTO 52
               DIAG(I,L)=A(K)
               GOTO 51
 52         CONTINUE
!--------------- append element not in any diagonal to AO,JAO,IAO -------------
            IF (JOB2.EQ.0) GOTO 51
            AO(KO)=A(K)
            JAO(KO)=J
            KO=KO+1
 51      CONTINUE
         IF (JOB2.NE.0) IND(I+1)=KO
 6    CONTINUE
      IF (JOB2.EQ.0) RETURN
!     finish with IAO
      IAO(1)=1
      DO 7 I=2,N+1
         IAO(I)=IND(I)
 7    CONTINUE
      RETURN
!----------- end of CSRDIA ----------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE CSRDIA
!******************************************************************************
!******************************************************************************
!******************************************************************************
SUBROUTINE AMUXD (N,X,Y,DIAG,NDIAG,IDIAG,IOFF)
IMPLICIT NONE
INTEGER                                    :: N,J,K,IO,I1,I2,NDIAG,IDIAG
INTEGER,DIMENSION(IDIAG)                   :: IOFF
REAL(8),DIMENSION(N)                       :: X,Y
REAL(8),DIMENSION(NDIAG,IDIAG)             :: DIAG
!------------------------------------------------------------------------------
!        A times a vector in Diagonal storage format (DIA)
!------------------------------------------------------------------------------
! multiplies a matrix by a vector when the original matrix is stored
! in the diagonal storage format.
!------------------------------------------------------------------------------
!
! on entry:
!----------
! N     = row dimension of A
! X     = real array of length equal to the column dimension of
!         the A matrix.
! NDIAG  = integer. The first dimension of array adiag as declared in
!         the calling program.
! IDIAG  = integer. The number of diagonals in the matrix.
! DIAG   = real array of size (NDIAG x IDIAG) containing the diagonals
!
! IOFF   = integer array of length IDIAG, containing the offsets of the
!   	   diagonals of the matrix:
!          DIAG(I,K) contains the element A(I,I+IOFF(K)) of the matrix.
!
! on return:
!-----------
! Y     = real array of length N, containing the product Y=A*X
!------------------------------------------------------------------------------
      DO 1 J=1,N
         Y(J)=0.D0
 1    CONTINUE
      DO 10 J=1,IDIAG
         IO=IOFF(J)
         I1=MAX0(1,1-IO)
         I2=MIN0(N,N-IO)
         DO 9 K=I1,I2
            Y(K)=Y(K)+DIAG(K,J)*X(K+IO)
 9       CONTINUE
 10   CONTINUE

      RETURN
!----------end-of-AMUXD--------------------------------------------------------
!------------------------------------------------------------------------------
END SUBROUTINE AMUXD
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                         CGS_SUBROUTINE
!______________________________________________________________________________
SUBROUTINE CGS (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,X_OLD,R_OLD,X)
USE PARAM, ONLY                    : TOL,SAVE_ITERATION_NORM
IMPLICIT NONE
INTEGER                            :: I,K,ITER,N,NNZERO,TL,ICODE,IDIAG,IERR
REAL(8),DIMENSION(NNZERO)          :: AA_AUX,AAR,LUVAL
INTEGER,DIMENSION(NNZERO)          :: JA_AUX,JAR
INTEGER,DIMENSION(2*N-1)           :: IND
INTEGER,DIMENSION(N+1)             :: IA_AUX,IAR
REAL(8),ALLOCATABLE                :: SA(:),COEF(:,:),DIAG(:,:)
INTEGER,ALLOCATABLE                :: IJA(:),JCOEF(:,:),IOFF(:)
REAL(8),DIMENSION(N)               :: X_OLD,X,P_OLD,P,R_OLD,U,U_OLD,R,AP,AUQ,RS,Q
REAL(8)                            :: NORM,S,ALPHA,BETA,M,MM,MN,HARVEST,TIME_BEGIN,TIME_END
!------------------------------------------------------------------------------
CALL TOTAL_LENGTH          (N,IA_AUX,TL)
ALLOCATE                   (SA(TL),IJA(TL))
CALL CSRMSR                (TL,N,NNZERO,AA_AUX,JA_AUX,IA_AUX,SA,IJA)
CALL TOTAL_LENGTH          (N,IA_AUX,TL)
CALL INFDIA                (N,NNZERO,JA_AUX,IA_AUX,IND,IDIAG)
ALLOCATE                   (COEF(N,IDIAG),JCOEF(N,IDIAG))
ALLOCATE                   (IOFF(IDIAG),DIAG(N,IDIAG))
CALL CSRELL                (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,N,COEF,JCOEF,N,IDIAG,IERR)
CALL CSRDIA                (N,NNZERO,IDIAG,10,AA_AUX,JA_AUX,IA_AUX,N,DIAG,IOFF,AAR,JAR,IAR,IND)
!------------------------------------------------------------------------------
!DO I=1,N
!CALL RANDOM_NUMBER (HARVEST)
!RS(I)=HARVEST
!END DO
RS=R_OLD
P_OLD=R_OLD
U_OLD=R_OLD
NORM=1.D0
ITER=0
!CALL CPU_TIME (TIME_BEGIN)
DO WHILE (NORM.GT.TOL)
ITER=ITER+1
!PRINT*,ITER
S=0.D0
M=0.D0
MM=0.D0
MN=0.D0
!------------------------------------------------------------------------------
! Stores the matrix in compressed spars row and multiplies by
! vector
!------------------------------------------------------------------------------
CALL CSR_MAT_V_PRODUCT     (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,P_OLD,AP)
!------------------------------------------------------------------------------
! Stores the matrix in modified spars row and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL SPRSAX (TL,N,SA,IJA,P_OLD,AP)
!------------------------------------------------------------------------------
! Stores the matrix in Ellpack/Itpack format and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL AMUXE (N,P_OLD,AP,N,IDIAG,COEF,JCOEF)
!------------------------------------------------------------------------------
! Stores the matrix in Diagonal format and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL AMUXD (N,P_OLD,AP,DIAG,N,IDIAG,IOFF)
!------------------------------------------------------------------------------
DO K=1,N
M=M+R_OLD(K)*RS(K)
MM=MM+AP(K)*RS(K)
END DO
ALPHA=M/MM
Q=U_OLD-ALPHA*AP
X=X_OLD+ALPHA*(U_OLD+Q)
!------------------------------------------------------------------------------
! Stores the matrix in compressed spars row and multiplies by
! vector
!------------------------------------------------------------------------------
CALL CSR_MAT_V_PRODUCT     (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,U_OLD+Q,AUQ)
!------------------------------------------------------------------------------
! Stores the matrix in modified spars row and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL SPRSAX (TL,N,SA,IJA,U_OLD+Q,AUQ)
!------------------------------------------------------------------------------
! Stores the matrix in Ellpack/Itpack format and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL AMUXE (N,U_OLD+Q,AUQ,N,IDIAG,COEF,JCOEF)
!------------------------------------------------------------------------------
! Stores the matrix in Diagonal format and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL AMUXD (N,U_OLD+Q,AUQ,DIAG,N,IDIAG,IOFF)
!------------------------------------------------------------------------------
R=R_OLD-ALPHA*AUQ
DO K=1,N
S=S+R(K)**2
END DO
NORM=SQRT(S)/N
DO K=1,N
MN=MN+R(K)*RS(K)
END DO
BETA=MN/M
U=R+BETA*Q
P=U+BETA*(Q+BETA*P_OLD)
P_OLD=P
X_OLD=X
R_OLD=R
U_OLD=U
!PRINT*,NORM
!CALL CPU_TIME (TIME_END)
!CALL CPU_TIME_WRITE (TIME_BEGIN,TIME_END,ITER)
IF (SAVE_ITERATION_NORM) CALL SUCCESSIVE_SOL ("CGS",NORM,ITER)
END DO
!------------------------------------------------------------------------------
END SUBROUTINE CGS
!______________________________________________________________________________
!******************************************************************************
!******************************************************************************
!******************************************************************************
!                         BCGSTAB_SUBROUTINE
!______________________________________________________________________________
SUBROUTINE BCGSTAB (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,X_OLD,R_OLD,X)
USE PARAM, ONLY                    : TOL,SAVE_ITERATION_NORM
IMPLICIT NONE
INTEGER                            :: I,K,ITER,N,NNZERO,TL,ICODE,IDIAG,IERR
REAL(8),DIMENSION(NNZERO)          :: AA_AUX,AAR,LUVAL
INTEGER,DIMENSION(NNZERO)          :: JA_AUX,JAR
INTEGER,DIMENSION(2*N-1)           :: IND
INTEGER,DIMENSION(N+1)             :: IA_AUX,IAR
REAL(8),ALLOCATABLE                :: SA(:),COEF(:,:),DIAG(:,:)
INTEGER,ALLOCATABLE                :: IJA(:),JCOEF(:,:),IOFF(:)
REAL(8),DIMENSION(N)               :: X_OLD,X,P_OLD,P,R_OLD,R,AP,AQQ,RS,QQ
REAL(8)                            :: NORM,S,ALPHA,BETA,M,MS,MMS,MM,MN,HARVEST,TIME_BEGIN,TIME_END,OM
!------------------------------------------------------------------------------
CALL TOTAL_LENGTH          (N,IA_AUX,TL)
ALLOCATE                   (SA(TL),IJA(TL))
CALL CSRMSR                (TL,N,NNZERO,AA_AUX,JA_AUX,IA_AUX,SA,IJA)
CALL TOTAL_LENGTH          (N,IA_AUX,TL)
CALL INFDIA                (N,NNZERO,JA_AUX,IA_AUX,IND,IDIAG)
ALLOCATE                   (COEF(N,IDIAG),JCOEF(N,IDIAG))
ALLOCATE                   (IOFF(IDIAG),DIAG(N,IDIAG))
CALL CSRELL                (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,N,COEF,JCOEF,N,IDIAG,IERR)
CALL CSRDIA                (N,NNZERO,IDIAG,10,AA_AUX,JA_AUX,IA_AUX,N,DIAG,IOFF,AAR,JAR,IAR,IND)
!------------------------------------------------------------------------------
!DO I=1,N
!CALL RANDOM_NUMBER (HARVEST)
!RS(I)=HARVEST
!END DO
RS=R_OLD
P_OLD=R_OLD
NORM=1.D0
ITER=0
!CALL CPU_TIME (TIME_BEGIN)
DO WHILE (NORM.GT.TOL)
ITER=ITER+1
!PRINT*,ITER
S=0.D0
M=0.D0
MS=0.D0
MMS=0.D0
MM=0.D0
MN=0.D0
!------------------------------------------------------------------------------
! Stores the matrix in compressed spars row and multiplies by
! vector
!------------------------------------------------------------------------------
CALL CSR_MAT_V_PRODUCT     (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,P_OLD,AP)
!------------------------------------------------------------------------------
! Stores the matrix in modified spars row and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL SPRSAX (TL,N,SA,IJA,P_OLD,AP)
!------------------------------------------------------------------------------
! Stores the matrix in Ellpack/Itpack format and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL AMUXE (N,P_OLD,AP,N,IDIAG,COEF,JCOEF)
!------------------------------------------------------------------------------
! Stores the matrix in Diagonal format and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL AMUXD (N,P_OLD,AP,DIAG,N,IDIAG,IOFF)
!------------------------------------------------------------------------------
DO K=1,N
M=M+R_OLD(K)*RS(K)
MM=MM+AP(K)*RS(K)
END DO
ALPHA=M/MM
QQ=R_OLD-ALPHA*AP
!------------------------------------------------------------------------------
! Stores the matrix in compressed spars row and multiplies by
! vector
!------------------------------------------------------------------------------
CALL CSR_MAT_V_PRODUCT     (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,QQ,AQQ)
!------------------------------------------------------------------------------
! Stores the matrix in modified spars row and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL SPRSAX (TL,N,SA,IJA,QQ,AQQ)
!------------------------------------------------------------------------------
! Stores the matrix in Ellpack/Itpack format and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL AMUXE (N,QQ,AQQ,N,IDIAG,COEF,JCOEF)
!------------------------------------------------------------------------------
! Stores the matrix in Diagonal format and multiplies by
! vector
!------------------------------------------------------------------------------
!CALL AMUXD (N,QQ,AQQ,DIAG,N,IDIAG,IOFF)
!------------------------------------------------------------------------------
DO K=1,N
MS=MS+AQQ(K)*QQ(K)
MMS=MMS+AQQ(K)*AQQ(K)
END DO
OM=MS/MMS
X=X_OLD+ALPHA*P_OLD+OM*QQ
R=QQ-OM*AQQ
DO K=1,N
S=S+R(K)**2
END DO
NORM=SQRT(S)/N
DO K=1,N
MN=MN+R(K)*RS(K)
END DO
BETA=(MN/M)*(ALPHA/OM)
P=R+BETA*(P_OLD-OM*AP)
P_OLD=P
X_OLD=X
R_OLD=R
!PRINT*,NORM
!CALL CPU_TIME (TIME_END)
!CALL CPU_TIME_WRITE (TIME_BEGIN,TIME_END,ITER)
IF (SAVE_ITERATION_NORM) CALL SUCCESSIVE_SOL ("BCGSTAB",NORM,ITER)
END DO
!------------------------------------------------------------------------------
END SUBROUTINE BCGSTAB
!******************************************************************************
!******************************************************************************
!                         ILU0_PBCGSTAB_SUBROUTINE
!******************************************************************************
!******************************************************************************
! This subroutine originally has been written by
! Iman Farahbakhsh on Jan 19 2010.
!==============================================================================
! This subroutine solve the linear system by a preconditioned
! Bi-Conjugate Gradient Stabilized and first time its algorithm developed by
! Henk Van der Vorst in 1992. The preconditioner is ILU(0) which stores in CSR
! format. The matrix-by-vector multiplication has been performed in 4 formats
! (CSR, MSR, Ellpack-Itpack and Diagonal) which can be determined by user.
!==============================================================================
! On Entry,
!============
!
! N          Dimension of matrix
! NNZERO     Number of nonzero elements of matrix
! AA_AUX     Real array of size NNZERO containing the nonzero elements of
!            matrix which are stored row by row
! JA_AUX     Integer array of size NNZERO containing the corresponding column
!            indices of nonzero elements
! IA_AUX     Integer array of size N+1 containing the pointers to the
!            beginning of each row
! X_OLD      Initial guess
! R_OLD      Initial residual
!==============================================================================
! On Exit
!============
! X          Solution
!==============================================================================
! Changes History,
!=================
!
!______________________________________________________________________________
 SUBROUTINE ILU0_PBCGSTAB (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,X_OLD,R_OLD,X)
 USE PARAM, ONLY                    : TOL,SAVE_ITERATION_NORM
 IMPLICIT NONE
 INTEGER                            :: I,K,ITER,N,NNZERO,TL,ICODE,IDIAG,IERR
 REAL(8),DIMENSION(NNZERO)          :: AA_AUX,AAR,LUVAL
 INTEGER,DIMENSION(NNZERO)          :: JA_AUX,JAR
 INTEGER,DIMENSION(2*N-1)           :: IND
 INTEGER,DIMENSION(N+1)             :: IA_AUX,IAR
 REAL(8),ALLOCATABLE                :: SA(:),COEF(:,:),DIAG(:,:)
 INTEGER,ALLOCATABLE                :: IJA(:),JCOEF(:,:),IOFF(:)
 REAL(8),DIMENSION(N)               :: X_OLD,X,P_OLD,P,R_OLD,R,RS,SR,SRHAT,T,V,V_OLD,PHAT,SHAT
 REAL(8)                            :: NORM,TIME_BEGIN,TIME_END,OM,OM_OLD,ALPHA,ALPHA_OLD,BETA,RHO,RHO_OLD,HARVEST
 INTEGER,DIMENSION(N)               :: UPTR,IW
!------------------------------------------------------------------------------
 CALL TOTAL_LENGTH          (N,IA_AUX,TL)
 ALLOCATE                   (SA(TL),IJA(TL))
 CALL CSRMSR                (TL,N,NNZERO,AA_AUX,JA_AUX,IA_AUX,SA,IJA)
 CALL TOTAL_LENGTH          (N,IA_AUX,TL)
 CALL INFDIA                (N,NNZERO,JA_AUX,IA_AUX,IND,IDIAG)
 ALLOCATE                   (COEF(N,IDIAG),JCOEF(N,IDIAG))
 ALLOCATE                   (IOFF(IDIAG),DIAG(N,IDIAG))
 CALL CSRELL                (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,N,COEF,JCOEF,N,IDIAG,IERR)
 CALL CSRDIA                (N,NNZERO,IDIAG,10,AA_AUX,JA_AUX,IA_AUX,N,DIAG,IOFF,AAR,JAR,IAR,IND)
!------------------------------------------------------------------------------
 CALL ILU0   (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,LUVAL,UPTR,IW,ICODE)
!------------------------------------------------------------------------------
 !DO I=1,N
 !CALL RANDOM_NUMBER (HARVEST)
 !RS(I)=HARVEST
 !END DO
 RS=R_OLD
 NORM=1.D0
 ITER=0
 !CALL CPU_TIME (TIME_BEGIN)
 DO WHILE (NORM.GT.TOL)
 ITER=ITER+1
! PRINT*,ITER
 RHO=DOT_PRODUCT(RS,R_OLD)
!------------------------------------------------------------------------------
 IF (RHO.EQ.0.D0) PRINT*, 'ILU0_PBiconjugate gradient stabilized method fails'
!------------------------------------------------------------------------------
 IF (ITER.EQ.1) THEN
 P=R_OLD
 ELSE
 BETA=(RHO/RHO_OLD)*(ALPHA_OLD/OM_OLD)
 P=R_OLD+BETA*(P_OLD-OM_OLD*V_OLD)
 END IF
!------------------------------------------------------------------------------
 CALL LUSOL  (N,NNZERO,P,PHAT,LUVAL,JA_AUX,IA_AUX,UPTR)
!------------------------------------------------------------------------------
! 1- Matrix by vector multiplication when matrix in CSR format
!------------------------------------------------------------------------------
CALL CSR_MAT_V_PRODUCT     (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,PHAT,V)
!------------------------------------------------------------------------------
! 2- Matrix by vector multiplication when matrix in MSR format
!------------------------------------------------------------------------------
!CALL SPRSAX (TL,N,SA,IJA,PHAT,V)
!------------------------------------------------------------------------------
! 3- Matrix by vector multiplication when matrix in Ellpack-Itpack format
!------------------------------------------------------------------------------
!CALL AMUXE (N,PHAT,V,N,IDIAG,COEF,JCOEF)
!------------------------------------------------------------------------------
! 4- Matrix by vector multiplication when matrix in Diagonal format
!------------------------------------------------------------------------------
!CALL AMUXD (N,PHAT,V,DIAG,N,IDIAG,IOFF)
!------------------------------------------------------------------------------
 ALPHA=RHO/DOT_PRODUCT(RS,V)
 SR=R_OLD-ALPHA*V
!------------------------------------------------------------------------------
 CALL LUSOL  (N,NNZERO,SR,SRHAT,LUVAL,JA_AUX,IA_AUX,UPTR)
!------------------------------------------------------------------------------
! 1- Matrix by vector multiplication when matrix in CSR format
!------------------------------------------------------------------------------
CALL CSR_MAT_V_PRODUCT     (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,SRHAT,T)
!------------------------------------------------------------------------------
! 2- Matrix by vector multiplication when matrix in MSR format
!------------------------------------------------------------------------------
!CALL SPRSAX (TL,N,SA,IJA,SRHAT,T)
!------------------------------------------------------------------------------
! 3- Matrix by vector multiplication when matrix in Ellpack-Itpack format
!------------------------------------------------------------------------------
!CALL AMUXE (N,SRHAT,T,N,IDIAG,COEF,JCOEF)
!------------------------------------------------------------------------------
! 4- Matrix by vector multiplication when matrix in Diagonal format
!------------------------------------------------------------------------------
!CALL AMUXD (N,SRHAT,T,DIAG,N,IDIAG,IOFF)
!------------------------------------------------------------------------------
 OM=DOT_PRODUCT(T,SR)/DOT_PRODUCT(T,T)
 X= X_OLD+ALPHA*PHAT+OM*SRHAT
 R=SR-OM*T
 NORM=DSQRT(DOT_PRODUCT(R,R))/N
 !PRINT*,ITER,NORM
 X_OLD=X
 R_OLD=R
 RHO_OLD=RHO
 P_OLD=P
 OM_OLD=OM
 V_OLD=V
 ALPHA_OLD=ALPHA
 !CALL CPU_TIME (TIME_END)
 !CALL CPU_TIME_WRITE (TIME_BEGIN,TIME_END,ITER)
 IF (SAVE_ITERATION_NORM) CALL SUCCESSIVE_SOL ("ILU0_PBCGSTAB",NORM,ITER)
 END DO
!------------------------------------------------------------------------------
 END SUBROUTINE ILU0_PBCGSTAB
!******************************************************************************
!******************************************************************************
!                         ILU0_PCGS_SUBROUTINE
!******************************************************************************
!******************************************************************************
! This subroutine originally has been written by
! Iman Farahbakhsh on Jan 19 2010.
!==============================================================================
! This subroutine solve the linear system by a preconditioned
! Conjugate Gradient Squared and first time its algorithm developed by
! Sonneveld in 1980. The preconditioner is ILU(0) which stores in CSR
! format. The matrix-by-vector multiplication has been performed in 4 formats
! (CSR, MSR, Ellpack-Itpack and Diagonal) which can be determined by user.
!==============================================================================
! On Entry,
!============
!
! N          Dimension of matrix
! NNZERO     Number of nonzero elements of matrix
! AA_AUX     Real array of size NNZERO containing the nonzero elements of
!            matrix which are stored row by row
! JA_AUX     Integer array of size NNZERO containing the corresponding column
!            indices of nonzero elements
! IA_AUX     Integer array of size N+1 containing the pointers to the
!            beginning of each row
! X_OLD      Initial guess
! R_OLD      Initial residual
!==============================================================================
! On Exit
!============
! X          Solution
!==============================================================================
! Changes History,
!=================
!
!______________________________________________________________________________
 SUBROUTINE ILU0_PCGS (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,X_OLD,R_OLD,X)
 USE PARAM, ONLY                    : TOL,SAVE_ITERATION_NORM
 IMPLICIT NONE
 INTEGER                            :: I,K,ITER,N,NNZERO,TL,ICODE,IDIAG,IERR
 REAL(8),DIMENSION(NNZERO)          :: AA_AUX,AAR,LUVAL
 INTEGER,DIMENSION(NNZERO)          :: JA_AUX,JAR
 INTEGER,DIMENSION(2*N-1)           :: IND
 INTEGER,DIMENSION(N+1)             :: IA_AUX,IAR
 REAL(8),ALLOCATABLE                :: SA(:),COEF(:,:),DIAG(:,:)
 INTEGER,ALLOCATABLE                :: IJA(:),JCOEF(:,:),IOFF(:)
 REAL(8),DIMENSION(N)               :: X_OLD,X,P_OLD,P,R_OLD,R,RS,Q,Q_OLD,UHAT,V,PHAT,U,QHAT
 REAL(8)                            :: NORM,ALPHA,BETA,HARVEST,TIME_BEGIN,TIME_END,RHO,RHO_OLD
 INTEGER,DIMENSION(N)               :: UPTR,IW
!------------------------------------------------------------------------------
 CALL TOTAL_LENGTH          (N,IA_AUX,TL)
 ALLOCATE                   (SA(TL),IJA(TL))
 CALL CSRMSR                (TL,N,NNZERO,AA_AUX,JA_AUX,IA_AUX,SA,IJA)
 CALL TOTAL_LENGTH          (N,IA_AUX,TL)
 CALL INFDIA                (N,NNZERO,JA_AUX,IA_AUX,IND,IDIAG)
 ALLOCATE                   (COEF(N,IDIAG),JCOEF(N,IDIAG))
 ALLOCATE                   (IOFF(IDIAG),DIAG(N,IDIAG))
 CALL CSRELL                (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,N,COEF,JCOEF,N,IDIAG,IERR)
 CALL CSRDIA                (N,NNZERO,IDIAG,10,AA_AUX,JA_AUX,IA_AUX,N,DIAG,IOFF,AAR,JAR,IAR,IND)
!------------------------------------------------------------------------------
 CALL ILU0   (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,LUVAL,UPTR,IW,ICODE)
 !-----------------------------------------------------------------------------
! DO I=1,N
! CALL RANDOM_NUMBER (HARVEST)
! RS(I)=HARVEST
! END DO
 RS=R_OLD
 NORM=1.D0
 ITER=0
 !CALL CPU_TIME (TIME_BEGIN)
DO WHILE (NORM.GT.TOL)
 ITER=ITER+1
 !PRINT*,ITER
 RHO=DOT_PRODUCT(RS,R_OLD)
!------------------------------------------------------------------------------
 IF (RHO.EQ.0.D0) PRINT*, 'conjugate gradient squared method fails'
!------------------------------------------------------------------------------
 IF (ITER.EQ.1) THEN
    U=R_OLD
    P=U
 ELSE
    BETA=RHO/RHO_OLD
    U=R_OLD+BETA*Q_OLD
	P=U+BETA*(Q_OLD+BETA*P_OLD)
 ENDIF
!------------------------------------------------------------------------------
 CALL LUSOL  (N,NNZERO,P,PHAT,LUVAL,JA_AUX,IA_AUX,UPTR)
!------------------------------------------------------------------------------
! 1- Matrix by vector multiplication when matrix in CSR format
!------------------------------------------------------------------------------
CALL CSR_MAT_V_PRODUCT     (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,PHAT,V)
!------------------------------------------------------------------------------
! 2- Matrix by vector multiplication when matrix in MSR format
!------------------------------------------------------------------------------
!CALL SPRSAX (TL,N,SA,IJA,PHAT,V)
!------------------------------------------------------------------------------
! 3- Matrix by vector multiplication when matrix in Ellpack-Itpack format
!------------------------------------------------------------------------------
!CALL AMUXE (N,PHAT,V,N,IDIAG,COEF,JCOEF)
!------------------------------------------------------------------------------
! 4- Matrix by vector multiplication when matrix in Diagonal format
!------------------------------------------------------------------------------
!CALL AMUXD (N,PHAT,V,DIAG,N,IDIAG,IOFF)
!------------------------------------------------------------------------------
 ALPHA=RHO/(DOT_PRODUCT(RS,V))
 Q=U-ALPHA*V
 CALL LUSOL  (N,NNZERO,U+Q,UHAT,LUVAL,JA_AUX,IA_AUX,UPTR)
 X=X_OLD+ALPHA*UHAT
!------------------------------------------------------------------------------
! 1- Matrix by vector multiplication when matrix in CSR format
!------------------------------------------------------------------------------
CALL CSR_MAT_V_PRODUCT     (N,NNZERO,AA_AUX,JA_AUX,IA_AUX,UHAT,QHAT)
!------------------------------------------------------------------------------
! 2- Matrix by vector multiplication when matrix in MSR format
!------------------------------------------------------------------------------
!CALL SPRSAX (TL,N,SA,IJA,UHAT,QHAT)
!------------------------------------------------------------------------------
! 3- Matrix by vector multiplication when matrix in Ellpack-Itpack format
!------------------------------------------------------------------------------
!CALL AMUXE (N,UHAT,QHAT,N,IDIAG,COEF,JCOEF)
!------------------------------------------------------------------------------
! 4- Matrix by vector multiplication when matrix in Diagonal format
!------------------------------------------------------------------------------
!CALL AMUXD (N,UHAT,QHAT,DIAG,N,IDIAG,IOFF)
!------------------------------------------------------------------------------
R=R_OLD-ALPHA*QHAT
NORM=DSQRT(DOT_PRODUCT(R,R))/N
P_OLD=P
X_OLD=X
R_OLD=R
Q_OLD=Q
RHO_OLD=RHO
!PRINT*,NORM
!CALL CPU_TIME (TIME_END)
!CALL CPU_TIME_WRITE (TIME_BEGIN,TIME_END,ITER)
IF (SAVE_ITERATION_NORM) CALL SUCCESSIVE_SOL ("ILU0_PCGS",NORM,ITER)
END DO
!------------------------------------------------------------------------------
END SUBROUTINE ILU0_PCGS
!******************************************************************************
!******************************************************************************




