module work_cgII

  USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im, imm,at, inv_angst, lperiod, rang, &
                       it, itesauv, itesauvposition, itesauvforce, &
                       inv_angst, erg2ev, angst, &
                       dmtype, potist
  use controle_mod
  use calfo_mod
  use analyse_mod
  use sauvegarde_mod
  use sauveposition_mod
  use sauveforce_mod
  implicit none


contains

  subroutine FUNCT(N,X,F,G,NCALLS,                      &
       xp_local, xpp, vp, ax, fp_local,  ielat, iwmax, ityp)
    use tab_imm_m, only : xp, fp
    double precision X(N),G(N),F
    integer i,N,NCALLS
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    integer  :: ielat(imm)
    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp_local(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    real(double)  :: fp_local(3,imm)
    !-----------------------------------------------

    

          IF (3*imm.NE.N) THEN
                  WRITE(0,'(a,i0)') "3*imm = ", 3*imm
                  WRITE(0,'(a,i0)') "N     = ", N
                  STOP "< work_cgII >"
          END IF
          IF (dmtype.EQ.30) THEN ! Variables = reduced coordinates
                  do i=1,im
                     xp_local(1:3,i) = MatMul( at, X(3*i-2:3*i) )
                  end do
          ELSE ! Variables = cartesian coordinates (in A)
                  do i=1,im
                     xp_local(1:3,i)=X(3*i-2:3*i)*inv_angst
                  end do
          END IF
          xp(:,:) =  xp_local(:,:)
    !back to internal units and JP world.......................................

    it=NCALLS-1

    if (lperiod)          call period 

    call controle 
    !write(*,*) 'inside FUNCT debug_incg1', xp_local(1,1), fp_local(1,1)
    call calfo
    fp_local(:,:) =  fp(:,:)
    !write(*,*) 'inside FUNCT debug_incg2', xp_local(1,1), fp_local(1,1)
    call analyse 
    !write(*,*) 'inside FUNCT debug_incg3', xp_local(1,1), fp_local(1,1)
   
    
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
                       G(3*i-2:3*i)=-MatMul(fp_local(:,i), at)*erg2eV
                    end do
            ELSE ! Variables = cartesian coordinates (in A)
                    do i=1,im
                       G(3*i-2:3*i)=-fp_local(1:3,i)*erg2eV/angst
                    end do
            END IF
            G(3*im+1:N)=0.d0
    

    return
  end subroutine FUNCT



end module work_cgII
