module work_cgII

  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  implicit none


contains

  subroutine FUNCT(N,X,F,G,NCALLS,                      &
       xp, xpp, vp, ax, fp,  ielat, iwmax, ityp)
    double precision X(N),G(N),F,forctot,formax  
    integer i,N,ncalls
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
    REAL(double) :: inv_angst
    INTEGER :: igc


    inv_angst = 1.d0/angst

    IF (lFrozen) THEN   ! Some atoms are fixed
          iGC=0
          do i=1 ,im
             IF (.Not.Free(i)) Cycle
             iGC=iGC+1
             xp(1:3,i) = X(3*iGC-2:3*iGC)*inv_angst
          end do
          IF (3*iGC.NE.N) THEN
                  WRITE(0,'(a,i0)') "3*iGC = ", 3*iGC
                  WRITE(0,'(a,i0)') "N     = ", N
                  STOP "< work_cgII >"
          END IF
    ELSE                ! All atoms are free to relax
          IF (3*imm.NE.N) THEN
                  WRITE(0,'(a,i0)') "3*imm = ", 3*imm
                  WRITE(0,'(a,i0)') "N     = ", N
                  STOP "< work_cgII >"
          END IF
            do i=1,N/3
               xp(1:3,i)=X(3*i-2:3*i)*inv_angst
            end do
    END IF

!    write(*,*) im,N/3,imm,N
!    stop
    !back to internal units and JP world.......................................

    it=NCALLS-1
    if (lperiod)          call period 
   
    call controle 
    
    call calfo
    call analyse 
    
    if (rang==0) then
       !     write(6,*)'analyse -> sauvegarde'
       if (itesauv/=0) then
          if (mod(it,itesauv)==0) call sauvegarde 
       endif

       !     write(6,*)'analyse -> sauveposition'
       if (itesauvposition/=0) then
          if (mod(it,itesauvposition)==0) call sauveposition ( it)
       endif
       !     write(6,*)'sauvposition -> control'
    endif                                   ! fin rang=0

    !go to into eV, ang and GC world............................................      

    F=potist*erg2eV
    do i=1,N/3
      !some contraintes
       !if ((i.eq.1).or.(i.eq.117)) then
       !G(3*i-2)= 0
       !G(3*i-1)= 0
       !G(3*i)  = 0
       !else 
       G(3*i-2)=-fp(1,i)*erg2eV/angst
       G(3*i-1)=-fp(2,i)*erg2eV/angst     
       G(3*i)  =-fp(3,i)*erg2eV/angst
       !end if       
    end do
    

    return
  end subroutine FUNCT



end module work_cgII
