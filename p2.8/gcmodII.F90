!-----------------------------------------------------------------------
!                                                                       
!   COMPUTER            - HP9000/DOUBLE                                 
!                                                                       
!   LATEST REVISION     - NOVEMBER 1, 1979                              
!                                                                       
!   PURPOSE             - A CONJUGATE GRADIENT ALGORITHM FOR FINDING    
!                           THE MINIMUM OF A FUNCTION OF N VARIABLES    
!                                                                       
!   USAGE               - CALL ZXCGR (FUNCT,N,ACC,MAXFN,DFPRED,X,G,F,W, 
!                           IER)                                        
!                                                                       
!   ARGUMENTS    FUNCT  - A USER SUPPLIED SUBROUTINE WHICH CALCULATES   
!                           THE OBJECTIVE FUNCTION AND ITS GRADIENT     
!                           FOR GIVEN PARAMETER VALUES                  
!                           X(1),X(2),...,X(N).                         
!                           THE CALLING SEQUENCE HAS THE FOLLOWING FORM 
!                           CALL FUNCT (N,X,F,G)                        
!                           WHERE X AND G ARE VECTORS OF LENGTH N.      
!                           THE SCALAR F IS FOR THE OBJECTIVE FUNCTION. 
!                           G(1), G(2), ..., G(N) ARE FOR THE COMPONENTS
!                           OF THE GRADIENT OF F.                       
!                           FUNCT MUST APPEAR IN AN EXTERNAL STATEMENT  
!                           IN THE CALLING PROGRAM. FUNCT MUST NOT      
!                           ALTER THE VALUES OF X(I),I=1,...,N OR N.    
!                N      - THE NUMBER OF PARAMETERS OF THE OBJECTIVE     
!                           FUNCTION. (INPUT) (I.E.,THE LENGTH OF X)    
!                ACC    - CONVERGENCE CRITERION. (INPUT)                
!                           THE CALCULATION ENDS WHEN THE SUM OF SQUARES
!                           OF THE COMPONENTS OF G IS LESS THAN ACC.    
!                MAXFN  - MAXIMUM NUMBER OF FUNCTION EVALUATIONS (I.E., 
!                           CALLS TO SUBROUTINE FUNCT) ALLOWED. (INPUT) 
!                           IF MAXFN IS SET TO ZERO, THEN THERE IS      
!                           NO RESTRICTION ON THE NUMBER OF FUNCTION    
!                           EVALUATIONS.                                
!                DFPRED - A ROUGH ESTIMATE OF THE EXPECTED REDUCTION    
!                           IN F, WHICH IS USED TO DETERMINE THE SIZE   
!                           OF THE INITIAL CHANGE TO X. (INPUT)         
!                           NOTE THAT DFPRED IS THE EXPECTED REDUCTION  
!                           ITSELF, AND DOES NOT DEPEND ON ANY RATIOS.  
!                           A BAD VALUE OF DFPRED CAUSES AN ERROR       
!                           MESSAGE, WITH IER=129, AND A RETURN ON THE  
!                           FIRST ITERATION. (SEE THE DESCRIPTION OF    
!                           IER BELOW)                                  
!                X      - VECTOR OF LENGTH N CONTAINING PARAMETER       
!                           VALUES.                                     
!                         ON INPUT, X MUST CONTAIN THE INITIAL          
!                           PARAMETER ESTIMATES.                        
!                         ON OUTPUT, X CONTAINS THE FINAL PARAMETER     
!                           ESTIMATES AS DETERMINED BY ZXCGR.           
!                G      - A VECTOR OF LENGTH N CONTAINING THE           
!                           COMPONENTS OF THE GRADIENT OF F AT THE      
!                           FINAL PARAMETER ESTIMATES. (OUTPUT)         
!                F      - A SCALAR CONTAINING THE VALUE OF THE FUNCTION 
!                           AT THE FINAL PARAMETER ESTIMATES. (OUTPUT)  
!                W      - WORK VECTOR OF LENGTH 6*N.                    
!                IER    - ERROR PARAMETER. (OUTPUT)                     
!                         IER = 0 IMPLIES THAT CONVERGENCE WAS          
!                           ACHIEVED AND NO ERRORS OCCURRED.            
!                         TERMINAL ERROR                                
!                           IER = 129 IMPLIES THAT THE LINE SEARCH OF   
!                             AN INTEGRATION WAS ABANDONED. THIS        
!                             ERROR MAY BE CAUSED BY AN ERROR IN THE    
!                             GRADIENT.                                 
!                           IER = 130 IMPLIES THAT THE CALCULATION      
!                             CANNOT CONTINUE BECAUSE THE SEARCH        
!                             DIRECTION IS UPHILL.                      
!                           IER = 131 IMPLIES THAT THE ITERATION WAS    
!                             TERMINATED BECAUSE MAXFN WAS EXCEEDED.    
!                           IER = 132 IMPLIES THAT THE CALCULATION      
!                             WAS TERMINATED BECAUSE TWO CONSECUTIVE    
!                             ITERATIONS FAILED TO REDUCE F.            
!		criterion	0 for sum < acc, 1 for max(g(I)**2) < acc
! 		a1, a2, a3	lattice vectors passed to FUNCT
!		N_A		number of atoms passed to FUNCT
!		M_A		masses of atoms passed to FUNCT
!		N_CONSTRAINTS	number of constraints passed to FUNCT
!		constraints	constraint array passed to FUNCT
!                                                                       
!   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
!                       - SINGLE/H36,H48,H60                            
!                                                                       
!   REQD. IMSL ROUTINES - UERTST,UGETIO                                 
!   (the call is commented out)
!                                                                       
!   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
!                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
!                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
!                                                                       
!   REMARKS  1.  THE ROUTINE INCLUDES NO THOROUGH CHECKS ON THE PART    
!                OF THE USER PROGRAM THAT CALCULATES THE DERIVATIVES    
!                OF THE OBJECTIVE FUNCTION. THEREFORE, BECAUSE          
!                DERIVATIVE CALCULATION IS A FREQUENT SOURCE OF         
!                ERROR, THE USER SHOULD VERIFY INDEPENDENTLY THE        
!                CORRECTNESS OF THE DERIVATIVES THAT ARE GIVEN TO       
!                THE ROUTINE.                                           
!            2.  BECAUSE OF THE CLOSE RELATION BETWEEN THE CONJUGATE    
!                GRADIENT METHOD AND THE METHOD OF STEEPEST DESCENTS,   
!                IT IS VERY HELPFUL TO CHOOSE THE SCALE OF THE          
!                VARIABLES IN A WAY THAT BALANCES THE MAGNITUDES OF     
!                THE COMPONENTS OF A TYPICAL DERIVATE VECTOR. IT        
!                CAN BE PARTICULARLY INEFFICIENT IF A FEW COMPONENTS    
!                OF THE GRADIENT ARE MUCH LARGER THAN THE REST.         
!            3.  IF THE VALUE OF THE PARAMETER ACC IN THE ARGUMENT      
!                LIST OF THE ROUTINE IS SET TO ZERO, THEN THE           
!                SUBROUTINE WILL CONTINUE ITS CALCULATION UNTIL IT      
!                STOPS REDUCING THE OBJECTIVE FUNCTION. IN THIS CASE    
!                THE USUAL BEHAVIOUR IS THAT CHANGES IN THE             
!                OBJECTIVE FUNCTION BECOME DOMINATED BY COMPUTER        
!                ROUNDING ERRORS BEFORE PRECISION IS LOST IN THE        
!                GRADIENT VECTOR. THEREFORE, BECAUSE THE POINT OF       
!                VIEW HAS BEEN TAKEN THAT THE USER REQUIRES THE         
!                LEAST POSSIBLE VALUE OF THE FUNCTION, A VALUE OF       
!                THE OBJECTIVE FUNCTION THAT IS SMALL DUE TO            
!                COMPUTER ROUNDING ERRORS CAN PREVENT FURTHER           
!                PROGRESS. HENCE THE PRECISION IN THE FINAL VALUES      
!                OF THE VARIABLES MAY BE ONLY ABOUT HALF THE            
!                NUMBER OF SIGNIFICANT DIGITS IN THE COMPUTER           
!                ARITHMETIC, BUT THE LEAST VALUE OF F IS USUALLY        
!                FOUND TO QUITE HIGH ACCURACY.                          
!                                                                       
!   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
!                                                                       
!   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
!                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
!                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
!                                                                       
!-----------------------------------------------------------------------
!                   G                                                    

module gcmodII_mod


#ifdef PARA
!  use mpi
  USE mod_para,only:MPI_COMM_space,ierr,nprocspace,status

#endif  
  implicit none
#ifdef PARA
  include 'mpif.h'
#endif  
contains

  SUBROUTINE ZXCGRII(FUNCT,N,ACC,MAXFN,X,G,F,W,IER,criterion,NCALLS)
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:dfpred,rang
    !  USE gen_com_m, ONLY:
    !                                  SPECIFICATIONS FOR ARGUMENTS         

    INTEGER            N,MAXFN,IER,iopt                              
    DOUBLE PRECISION   ACC,X(N),G(N),F,W(6*N)                    
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    !-----------------------------------------------
    integer criterion
    integer NCALLS
    logical do_print, do_print_verbose
    !                                  SPECIFICATIONS FOR LOCAL VARIABLES   
    INTEGER            MAXLIN,MXFCON,I,IGINIT,IGOPT,IRETRY,IRSDG,         &
         IRSDX,ITERC,ITERFM,ITERRS,IXOPT,NFBEG,NFOPT                                          
    DOUBLE PRECISION   BETA,DDSPLN,DFPR,FCH,FINIT,FMIN,GAMDEN,GAMA,       &  
         GINIT,GMIN,GNEW,GSPLN,GSQRD,SBOUND,STEP,STEPCH,    &  
         STMIN,SUM,WORK                                 
    double precision g_elem_max
    integer ::ip
    external funct
    DATA               MAXLIN/5/,MXFCON/2/                            
    !                                  FIRST EXECUTABLE STATEMENT           
    IER = 0                                                           
    !                                  THE WORKING SPACE ARRAY IS SPLIT     
    !                                    INTO SIX VECTORS OF LENGTH N. THE  
    !                                    FIRST PART IS USED FOR THE SEARCH  
    !                                    DIRECTION OF AN ITERATION. THE     
    !                                    SECOND AND THIRD PARTS CONTAIN THE 
    !                                    INFORMATION THAT IS REQUIRED BY    
    !                                    THE CONJUGACY CONDITIONS OF THE    
    !                                    RESTART PROCEDURE. THE FOURTH PART 
    !                                    CONTAINS THE GRADIENT AT THE START 
    !                                    OF AN ITERATION. THE FIFTH PART    
    !                                    CONTAINS THE PARAMETERS THAT GIVE  
    !                                    THE LEAST CALCULATED VALUE OF F.   
    !                                    THE SIXTH PART CONTAINS THE        
    !                                    GRADIENT VECTOR WHERE F IS LEAST.  
    IRSDX = N                                                         
    IRSDG = IRSDX+N                                                   
    IGINIT = IRSDG+N                                                  
    IXOPT = IGINIT+N                                                  
    IGOPT = IXOPT+N                                                   
    !                                  SET SOME PARAMETERS TO BEGIN THE     
    !                                    CALCULATION. ITERC AND             
    !                                    NCALLS COUNT THE NUMBER OF         
    !                                    ITERATIONS AND CALLS OF FUNCT.     
    !                                    ITERFM IS THE NUMBER OF THE MOST   
    !                                    RECENT ITERATION THAT DECREASES F. 
    ITERC = 0                                                         
    NCALLS = 0                                                        
    ITERFM = ITERC                                                    
    !                                  CALL SUBROUTINE FUNCT. LET THE       
    !                                    INITIAL SEARCH DIRECTION BE MINUS  
    !                                    THE GRADIENT VECTOR. USUALLY THE   
    !                                    PARAMETER ITERRS GIVES THE         
    !                                    ITERATION NUMBER OF THE MOST       
    !                                    RECENT RESTART, BUT IT IS SET TO   
    !                                    ZERO WHEN THE STEEPEST DESCENT     
    !                                    DIRECTION IS USED.   

    do_print=.true.      



5   NCALLS = NCALLS+1 
!    if (rang==0) then
!       write(6,*)'X1',x
!       write(6,*)'G1',G
!    end if
    CALL FUNCT (N,X,F,G,NCALLS)                                



    if (rang==0) then
!       write(6,*)'X',x
!       write(6,*)'G',G
       
       if (NCALLS .eq. 1) then
          !cos	if (do_print) print 3000, F
          !cos 3000    FORMAT ("ZXCGR Starting energy = ", F20.10)
       end if
       IF (NCALLS.GE.2) GO TO 20                                         
10     DO I=1,N                                                       
          W(I) = -G(I)
       end DO
       ITERRS = 0                                                        
       IF (ITERC.GT.0) GO TO 80                                          
       !                                  SET SUM TO G SQUARED. GMIN AND GNEW  
       !                                    ARE THE OLD AND THE NEW            
       !                                    DIRECTIONAL DERIVATIVES ALONG THE  
       !                                    CURRENT SEARCH DIRECTION. LET FCH  
       !                                    BE THE DIFFERENCE BETWEEN F AND    
       !                                    THE PREVIOUS BEST VALUE OF THE     
       !                                    OBJECTIVE FUNCTION.                
20     GNEW = 0.0D0                                                      
       SUM = 0.0D0                                                       
       g_elem_max = 0.0
       DO I=1,N                                                       
          GNEW = GNEW+W(I)*G(I)                                          
          if (G(I)**2 .gt. g_elem_max) g_elem_max = G(I)**2
25        SUM = SUM+G(I)**2
       end DO
       IF (NCALLS.EQ.1) GO TO 35                                         
       FCH = F-FMIN                                                      
       !                                  STORE THE VALUES OF X, F AND G, IF   
       !                                    THEY ARE THE BEST THAT HAVE BEEN   
       !                                    CALCULATED SO FAR, AND NOTE G      
       !                                    SQUARED AND THE VALUE OF NCALLS.   
       !                                    TEST FOR CONVERGENCE.              
       IF (FCH) 35,30,50                                                 
30     IF (GNEW/GMIN.LT.-1.0D0) GO TO 45                                 
35     FMIN = F                                                          
       GSQRD = SUM                                                       
       NFOPT = NCALLS                                                    
       DO I=1,N                                                       
          W(IXOPT+I) = X(I)                                              
40        W(IGOPT+I) = G(I)
       end DO
45     continue !cos    if (do_print) print 3010, ITERC,NCALLS, sqrt(SUM), F
       !cos 3010 FORMAT ("ZXCGR:",2I5,2F20.10)
       IF ((SUM.LE.ACC .and. criterion.eq.0) .or.                 &
            (g_elem_max.LE.ACC .and. criterion.eq.1)) GO TO 9005      
       !                                  TEST IF THE VALUE OF MAXFN ALLOWS    
       !                                    ANOTHER CALL OF FUNCT.             
50     IF (NCALLS.NE.MAXFN) GO TO 55                                     
       IER = 131                                                         
       GO TO 9000                                                        
55     IF (NCALLS.GT.1) GO TO 100                                        
       !                                  SET DFPr TO THE ESTIMATE OF THE      
       !                                    REDUCTION IN F GIVEN IN THE        
       !                                    ARGUMENT LIST, IN ORDER THAT THE   
       !                                    INITIAL CHANGE TO THE PARAMETERS   
       !                                    IS OF A SUITABLE SIZE. THE VALUE   
       !                                    OF STMIN IS USUALLY THE            
       !                                    STEP-LENGTH OF THE MOST RECENT     
       !                                    LINE SEARCH THAT GIVES THE LEAST   
       !                                    CALCULATED VALUE OF F.             
       DFPR = dfpred                                                    
       STMIN = dfpred/GSQRD                                              
       !      if (do_print) print *, "DFPR, STMIN = ", DFPR, STMIN
       !                                  BEGIN THE ITERATION                  
80     ITERC = ITERC+1                                                   
       !                                  STORE THE INITIAL FUNCTION VALUE AND 
       !                                    GRADIENT, CALCULATE THE INITIAL    
       !                                    DIRECTIONAL DERIVATIVE, AND BRANCH 
       !                                    IF ITS VALUE IS NOT NEGATIVE. SET  
       !                                    SBOUND TO MINUS ONE TO INDICATE    
       !                                    THAT A BOUND ON THE STEP IS NOT    
       !                                    KNOWN YET, AND SET NFBEG TO THE    
       !                                    CURRENT VALUE OF NCALLS. THE       
       !                                    PARAMETER IRETRY SHOWS THE NUMBER  
       !                                    OF ATTEMPTS AT SATISFYING THE BETA 
       !                                    CONDITION.                         
       FINIT = F                                                         
       GINIT = 0.0D0                                                     
       DO  I=1,N                                                       
          W(IGINIT+I) = G(I)                                             
85        GINIT = GINIT+W(I)*G(I)
       end DO
       !      if (do_print) print *, "GINIT = ", GINIT
       IF (GINIT.GE.0.0D0) GO TO 165                                     
       GMIN = GINIT                                                      
       SBOUND = -1.0D0                                                   
       NFBEG = NCALLS                                                    
       IRETRY = -1                                                       
       !                                  SET STEPCH SO THAT THE INITIAL       
       !                                    STEP-LENGTH IS CONSISTENT WITH THE 
       !                                    PREDICTED REDUCTION IN F, SUBJECT  
       !                                    TO THE CONDITION THAT IT DOES NOT  
       !                                    EXCEED THE STEP-LENGTH OF THE      
       !                                    PREVIOUS ITERATION. LET STMIN BE   
       !                                    THE STEP TO THE LEAST CALCULATED   
       !                                    VALUE OF F.                        
       !      if (do_print) print *, "setting STEPCH ", STMIN, DFPR/GINIT
       STEPCH = DMIN1(STMIN,DABS(DFPR/GINIT))                            
       STMIN = 0.0D0                                                     
       !                                  CALL SUBROUTINE FUNCT AT THE VALUE   
       !                                    OF X THAT IS DEFINED BY THE NEW    
       !                                    CHANGE TO THE STEP-LENGTH, AND LET 
       !                                    THE NEW STEP-LENGTH BE STEP. THE   
       !                                    VARIABLE WORK IS USED AS WORK      
       !                                    SPACE.                             
90     STEP = STMIN+STEPCH                                               
       WORK = 0.0D0                                                      
       DO  I=1,N                                                       
          X(I) = W(IXOPT+I)+STEPCH*W(I)                                  
95        WORK = DMAX1(WORK,DABS(X(I)-W(IXOPT+I)))
       end DO
       IF (WORK.GT.0.0D0) then
#ifdef PARA
          if (nprocspace.gt.1) then
             iopt=1
             do ip=1,nprocspace-1
                call MPI_SEND(iopt,  1, MPI_INTEGER, ip, 10001, MPI_COMM_space, status, ierr)
             end do
          end if
#endif                         

          GO TO 5     !!!!!!!!!!!      !!!!!!!!!!!!                                      
       end IF
       !                                  TERMINATE THE LINE SEARCH IF STEPCH  
       !                                    IS EFFECTIVELY ZERO.               
       IF (NCALLS.GT.NFBEG+1) GO TO 115                                  
       IF (DABS(GMIN/GINIT)-0.2D0) 170,170,115                           
       !                                  LET SPLN BE THE QUADRATIC SPLINE     
       !                                    THAT INTERPOLATES THE CALCULATED   
       !                                    FUNCTION VALUES AND DIRECTIONAL    
       !                                    DERIVATIVES AT THE POINTS STMIN    
       !                                    AND STEP OF THE LINE SEARCH, WHERE 
       !                                    THE KNOT OF THE SPLINE IS AT       
       !                                    0.5*(STMIN+STEP). REVISE STMIN,    
       !                                    GMIN AND SBOUND, AND SET DDSPLN TO 
       !                                    THE SECOND DERIVATIVE OF SPLN AT   
       !                                    THE NEW STMIN. HOWEVER, IF FCH IS  
       !                                    ZERO, IT IS ASSUMED THAT THE       
       !                                    MAXIMUM ACCURACY IS ALMOST         
       !                                    ACHIEVED, SO DDSPLN IS CALCULATED  
       !                                    USING ONLY THE CHANGE IN THE       
       !                                    GRADIENT.                          
100    WORK = (FCH+FCH)/STEPCH-GNEW-GMIN                                 
       DDSPLN = (GNEW-GMIN)/STEPCH                                       
       IF (NCALLS.GT.NFOPT) SBOUND = STEP                                
       IF (NCALLS.GT.NFOPT) GO TO 105                                    
       IF (GMIN*GNEW.LE.0.0D0) SBOUND = STMIN                            
       STMIN = STEP                                                      
       GMIN = GNEW                                                       
       STEPCH = -STEPCH                                                  
105    IF (FCH.NE.0.0D0) DDSPLN = DDSPLN+(WORK+WORK)/STEPCH              
       !                                                                       
       !                                  TEST FOR CONVERGENCE OF THE LINE     
       !                                    SEARCH, BUT FORCE AT LEAST TWO     
       !                                    STEPS TO BE TAKEN IN ORDER NOT TO  
       !                                    LOSE QUADRATIC TERMINATION.        
       IF (GMIN.EQ.0.0D0) GO TO 170                                      
       IF (NCALLS.LE.NFBEG+1) GO TO 120                                  
       IF (DABS(GMIN/GINIT).LE.0.2D0) GO TO 170                          
       !                                  APPLY THE TEST THAT DEPENDS ON THE   
       !                                    PARAMETER MAXLIN.                  
110    IF (NCALLS.LT.NFOPT+MAXLIN) GO TO 120                             
115    IER = 129                                                         
       GO TO 170                                                         
       !                                  SET STEPCH TO THE GREATEST CHANGE TO 
       !                                    THE CURRENT VALUE OF STMIN THAT IS 
       !                                    ALLOWED BY THE BOUND ON THE LINE   
       !                                    SEARCH. SET GSPLN TO THE GRADIENT  
       !                                    OF THE QUADRATIC SPLINE AT         
       !                                    (STMIN+STEPCH). HENCE CALCULATE    
       !                                    THE VALUE OF STEPCH THAT MINIMIZES 
       !                                    THE SPLINE FUNCTION, AND THEN      
       !                                    OBTAIN THE NEW FUNCTION AND        
       !                                    GRADIENT VECTOR, FOR HIS VALUE OF 
       !                                    THE CHANGE TO THE STEP-LENGTH.     
120    STEPCH = 0.5D0*(SBOUND-STMIN)                                     
       IF (SBOUND.LT.-0.5D0) STEPCH = 9.0D0*STMIN                        
       GSPLN = GMIN+STEPCH*DDSPLN                                        
       IF (GMIN*GSPLN.LT.0.0D0) STEPCH = STEPCH*GMIN/(GMIN-GSPLN)        
       GO TO 90                                                          
       !                                  CALCULATE THE VALUE OF BETA THAT     
       !                                    OCCURS IN THE NEW SEARCH           
       !                                    DIRECTION.                         
125    SUM = 0.0D0                                                       
       g_elem_max = 0.0D0
       DO  I=1,N                                                      
          if (G(I)*W(IGINIT+I) .gt. g_elem_max) g_elem_max = G(I)*W(IGINIT+I)
130       SUM = SUM+G(I)*W(IGINIT+I)
       end DO
       BETA = (GSQRD-SUM)/(GMIN-GINIT)                                   
       !                                  TEST THAT THE NEW SEARCH DIRECTION   
       !                                    CAN BE MADE DOWNHILL. IF IT        
       !                                    CANNOT, THEN MAKE ONE ATTEMPT TO   
       !                                    IMPROVE THE ACCURACY OF THE LINE   
       !                                    SEARCH.                            
       IF (DABS(BETA*GMIN).LE.0.2D0*GSQRD) GO TO 135                     
       IRETRY = IRETRY+1                                                 
       IF (IRETRY.LE.0) GO TO 110                                        
       !                                  APPLY THE TEST THAT DEPENDS ON THE   
       !                                    PARAMETER MXFCON.                  
       !                                    SET DFPR TO THE PREDICTED          
       !                                    REDUCTION IN F ON THE NEXT         
       !                                    ITERATION.                         
135    IF (F.LT.FINIT) ITERFM = ITERC                                    
       IF (ITERC.LT.ITERFM+MXFCON) GO TO 140                             
       IER = 132                                                         
       GO TO 9000                                                        
140    DFPR = STMIN*GINIT                                                
       !                                  BRANCH IF A RESTART PROCEDURE IS     
       !                                    REQUIRED DUE TO THE ITERATION      
       !                                    NUMBER OR DUE TO THE SCALAR        
       !                                    PRODUCT OF CONSECUTIVE GRADIENTS.  
       IF (IRETRY.GT.0) GO TO 10                                         
       IF (ITERRS.EQ.0) GO TO 155                                        
       IF (ITERC-ITERRS.GE.N) GO TO 155                                  
       IF (DABS(SUM).GE.0.2D0*GSQRD) GO TO 155                           
       !                                  CALCULATE THE VALUE OF GAMA THAT     
       !                                    OCCURS IN THE NEW SEARCH           
       !                                    DIRECTION, AND SET SUM TO A SCALAR 
       !                                    PRODUCT FOR THE TEST BELOW. THE    
       !                                    VALUE OF GAMDEN IS SET BY THE      
       !                                    RESTART PROCEDURE.                 
       GAMA = 0.0D0                                                      
       SUM = 0.0D0                                                       
       g_elem_max = 0.0D0
       DO  I=1,N                                                      
          if (G(I)*W(IRSDX+I) .gt. g_elem_max) g_elem_max = G(I)*W(IRSDX+I)
          GAMA = GAMA+G(I)*W(IRSDG+I)                                    
145       SUM = SUM+G(I)*W(IRSDX+I)
       end DO
       GAMA = GAMA/GAMDEN                                                
       !                                  RESTART IF THE NEW SEARCH DIRECTION  
       !                                    IS NOT SUFFICIENTLY DOWNHILL.      
       !                                                                       
       IF (DABS(BETA*GMIN+GAMA*SUM).GE.0.2D0*GSQRD) GO TO 155            
       !                                                                       
       !                                  CALCULATE THE NEW SEARCH DIRECTION.  
       DO I=1,N                                                      
150       W(I) = -G(I)+BETA*W(I)+GAMA*W(IRSDX+I)
       end DO
       GO TO 80                                                          
       !                                  APPLY THE RESTART PROCEDURE.         
155    GAMDEN = GMIN-GINIT                                               
       DO  I=1,N                                                      
          W(IRSDX+I) = W(I)                                              
          W(IRSDG+I) = G(I)-W(IGINIT+I)                                  
160       W(I) = -G(I)+BETA*W(I)
       end DO
       ITERRS = ITERC                                                    
       GO TO 80                                                          
       !                                  SET IER TO INDICATE THAT THE SEARCH  
       !                                    DIRECTION IS UPHILL.               
165    IER = 130                                                         
       !                                  ENSURE THAT F, X AND G ARE OPTIMAL.  
170    IF (NCALLS.EQ.NFOPT) GO TO 180                                    
       F = FMIN                                                          
       DO  I=1,N                                                      
          X(I) = W(IXOPT+I)                                              
175       G(I) = W(IGOPT+I)
       end DO
180    IF (IER.EQ.0) GO TO 125                                           
9000   CONTINUE                                                          
       !      CALL UERTST (IER,'ZXCGR ')                                        
       if (do_print) print 3020
3020   FORMAT ("came from 9000")
9005   if (do_print) print 3030,NCALLS
3030   FORMAT ("NCALLS",I5)
#ifdef PARA
       if (nprocspace.gt.1) then

          iopt=0
          do ip=1,nprocspace-1
             call MPI_SEND(iopt,  1, MPI_INTEGER, ip, 10001, MPI_COMM_space, status, ierr)
          end do
       end if
#endif                                     
       RETURN                                                            
       !END DO
    else
#ifdef PARA
       if (nprocspace.gt.1) then

          call MPI_RECV(iopt,  1, MPI_INTEGER, 0, 10001, MPI_COMM_space, status, ierr)
          select case (iopt)
          case (0)
             return
          case(1)
             goto 5
          end select
       end if
#endif
    endif
  end subroutine ZXCGRII
   end module
