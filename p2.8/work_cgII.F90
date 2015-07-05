module work_cgII

  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  implicit none


contains

  subroutine FUNCT(N,X,F,G,NCALLS,                      &
       xp, xpp, vp, ax, fp,  ielat, iwmax, ityp)
    double precision X(N),G(N),F,forctot,formax  
    integer i,N,NCALLS
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: ielat(imm)
    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    real(double)  :: fp(3,imm)
    !-----------------------------------------------
    INTEGER :: IGC

    

          IF (3*imm.NE.N) THEN
                  WRITE(0,'(a,i0)') "3*imm = ", 3*imm
                  WRITE(0,'(a,i0)') "N     = ", N
                  STOP "< work_cgII >"
          END IF
          IF (dmtype.EQ.30) THEN ! Variables = reduced coordinates
                  do i=1,im
                     xp(1:3,i) = MatMul( at, X(3*i-2:3*i) )
                  end do
          ELSE ! Variables = cartesian coordinates (in A)
                  do i=1,im
                     xp(1:3,i)=X(3*i-2:3*i)*inv_angst
                  end do
          END IF

    !back to internal units and JP world.......................................

    it=NCALLS-1

    if (lperiod)          call period 

    call controle 
    call calfo
    call analyse 
   
    
    if (it.ne.0) then
    if (rang==0) then
!           write(6,*)'work_cg_II analyse -> sauvegarde',it
       if (itesauv.GT.0) then
          if (mod(it,itesauv)==0) call sauvegarde 
       endif

!            write(6,*)'work_cg_II analyse -> sauveposition',it
       if (itesauvposition.GT.0) then
          if (mod(it,itesauvposition)==0) call sauveposition ( it)
       endif
       if (itesauvforce.GT.0) then
          if (mod(it,itesauvforce)==0) call sauveforce ( it)
       endif
!            write(6,*)'work_cg_II sauvposition -> control',it
    endif                                   ! fin rang=0
    end if
    !go to into eV, ang and GC world............................................      

    
    F=potist*erg2eV
            IF (dmtype.EQ.30) THEN ! Variables = reduced coordinates
                    do i=1,im
                       G(3*i-2:3*i)=-MatMul(fp(:,i), at)*erg2eV
                    end do
            ELSE ! Variables = cartesian coordinates (in A)
                    do i=1,im
                       G(3*i-2:3*i)=-fp(1:3,i)*erg2eV/angst
                    end do
            END IF
            G(3*im+1:N)=0.d0
    

    return
  end subroutine FUNCT



end module work_cgII
