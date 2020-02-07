module dyn_vverlet_mod
        use calfo_mod
        use calfoberend_mod 
        use caltabt_mod
#ifdef PARA
        use layer_mod 
#endif
        implicit none
        contains
! *************************************************************
subroutine dyn_vverlet
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
use tab_imm_m
use jqmod
use suivinonpbc
use elec_cell,only: dynelec
#ifdef PARA
  use mod_para
#endif

  use elec_cell, only:TTlangevin
  use Parrinello_Rahman

  implicit none
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, ia,ic,il
  real(double), dimension(ntyp) :: aux
  real(double), save :: tmoyinst, imesureT
  !-----------------------------------------------
  !real(double), external :: tempinst
  real(double)::eatommoy
#ifdef PARA
  real(double)::jq_tot(3)
#endif

  !      write (*,*) 'sub dynvverlet'

  timel = timel+tstep
  aux(:ntyp) = tstep/cm(:ntyp)/2.d0


! step 1 First half-step velocities update, v(t) -> v(t+dt/2)



  if (dmtype==2) then
     do i = 1, im
        do ic = 1, 3
           if (vp(ic,i)*fp(ic,i)<0) then
              vp(ic,i)=0.
           end if
        end do
     end do
 end if

 if (lLangevin) then
    il=2*(ilangevin-1)+1
    call dynlangevin(xp,vp,fp,ityp,il,Gl)
 elseif (l2T) then
    il=2*(ilangevin-1)+1
    call TTlangevin(xp,vp,fp,ityp,il,Gl)
 else
   DO i=1, imd
       vp(1:3,i) = vp(1:3,i) + aux(iTyp(i))*fp(1:3,i)
    END DO
 end if


!step 2  Coordinate update, x(t)-> x(t+dt)

  DO i=1, imd
     xpp(1:3,i)=xp(1:3,i)
     xp(1:3,i) = xp(1:3,i) + tstep*vp(1:3,i)
  END DO


  !conditions periodiques
  if (lsuivinonpbc) then
     DO i=1,imd
       tmpsuivi(1:3,i)=tmpsuivi(1:3,i)+tstep*vp(1:3,i)
     END DO 
  end if
  if (lperiod)  call period 

  ! repartition des atomes dans la nouvelle boite


  ! cell dispatching

  if (.not.lpr) then
     if (itab/=0) then
        if (mod(it,itab)==0) then
           call caltabt
        endif
     endif
  end if


#ifdef PARA
  temps_debpara=MPI_Wtime()
	   ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
	   call maj_atomes_frt_ftm
  temps_para=temps_para+MPI_Wtime()-temps_debpara
  if (ltranche) call layer
#endif


  if (l2T) then
     call dynelec
  end if

  ! Force calculation

  call calfo   ! F(t+dt)

  if (lnemd) then
     eatommoy=0.
     do i=1,imd
        eatommoy=eatommoy+eatom(i)/float(imd)
     end do
     do i=1,imd
        fp(1,i)=fp(1,i)+(eatom(i)-eatommoy)*Fnemd
     end do
  end if



 
  ! Second half-step velocities update, v(t+1/2dt) -> v(t+dt)
  if (llangevin.eqv..true.) then
    il=2*(ilangevin-1)+2
     call dynlangevin(xp,vp,fp,ityp,il,Gl)
  elseif (l2T) then
     il=2*(ilangevin-1)+2
     call TTlangevin(xp,vp,fp,ityp,il,Gl)
  else
     DO i=1, imd
        vp(1:3,i) = vp(1:3,i) + aux(iTyp(i))*fp(1:3,i)
     END DO
  end if

  if (associated(eatom))  eatom(1:im)=eatom(1:im)+0.5*cm(ityp(1:im))*(vp(1,1:im)**2+vp(2,1:im)**2+vp(3,1:im)**2)

  if (lcalcjq) then
     jqp=jq ; jqk=0.0 !; expvect(:)=0.0
     do i=1,imd
!        eatom(i)=eatom(i)+0.5*cm(ityp(i))*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
        expvect(:)=expvect(:)+eatom(i)*xpp(:,i)
        do ic=1,3
           jqk(ic)=jqk(ic)+eatom(i)*vp(ic,i)
           !               jq(ic)=jq(ic)+eatom(i)*vp(ic,i)
        end do
     end do
     jq=jqp+jqk

#ifdef PARA
     call MPI_ALLREDUCE(jq,jq_tot,3,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     jq=jq_tot
#endif



     if(rang==0)write(65,'(I8,3D16.8)')it-1,jq(1),jq(2),jq(3)
     !         write(66,'(I8,3D15.6)')it-1,expvect(1),expvect(2),expvect(3)
  end if



  return
end subroutine dyn_vverlet
end module
