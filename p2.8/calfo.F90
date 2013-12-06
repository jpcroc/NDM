! ************************************************
!           Sous-programme calfo
!routine d'appel des routines de forces
! ************************************************

subroutine calfo

  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use contrainte
  use tab_imm_m
  use jqmod
#if(PARA)
  use mod_mpi
#endif

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real(double), dimension(3) :: fptot
  integer :: i,ilocal,ipot
  real(double)::vn,ic,v1,f1
  integer::nv1

#if(PARA)
  real(double), dimension(3,3) :: sig_tot,sigkine_tot
  real(double),dimension (3):: fptot_tot
  real(double)::elosselectot,elosselec1tot
#endif
  !  if (rang==0) write(6,*) 'PARA-T entree calfo'
  sig(:,:)=0.d0 ; if (ltpcel.EQV..true.) sigc=0
  potist=0.
  potis1=0. ; potis2=0.; potis3=0.; potis0=0. ; potcp=0.; potisP=0.
  potisTersoff=0.
  potisrep=0.; potisglue=0.; potiseam=0.
  if(lsigat) sigat(:,:,:)=0. ; 


  if (lsigtyp) then
     sigtyp=0. ; sigtyptyp=0.
#if(PARA)
     sigtyp_loc=0.;     sigtyptyp_loc=0.
#endif 
  end if



  fp(:,:) = zero
  jq=0.0 
  if (associated(eatom)) eatom(:)=0


  do ipot=0,npotmax
     if (lpotentiel(ipot).EQV..true.) then
        ipotentiel=ipot
        if(ipotentiel.lt.10) then
           select case (ipotentiel)
           case(0,1,3,4,5,6,7)
              if (ltabvois) then
                 call calfo2ctabvois (xp, vp,  fp, iwmax, ityp)
              else
                 call calfo2ccel 

              endif


              if (iewald.ge.1) call calfoew

              ! Potentiel total
              potisP = potis0+potis1+potis2+potis3
              potist=potist+potisP

           case(2)
              ! !!! le cas parallele n'est pas pris en compte !!!

              if (.not.parallele) call calfow(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)

           case default
           end select

           ! !!! le cas parallele n'est pas pris en compte !!!
           if (.not.parallele.and.l3c) call calfo3c (xp,  vp,  fp, ielat, iwmax, ityp)

           ! !!! le cas parallele n'est pas pris en compte !!!
           !potentiels EAM
        else
           select case (ipotentiel)
           case(12)
              ! !!! le cas parallele n'est pas pris en compte !!!
              if (.not.parallele) call calfojuli(xp,  vp, fp, ielat, iwmax, ityp)
           case(13,14,15)
              if (ltabvois) then
                 ! !!! le cas parallele n'est pas pris en compte !!!
                 if (.not.parallele) call force_tersoff (xp,  vp, fp, iwmax, ityp)
              else
                 call force_tersoff_cel
              endif
              potist=potist+potisTersoff
           case (10,11)
              if (ltabvois) then
                 ! !!! le cas parallele n'est pas pris en compte !!!
                 if (.not.parallele) then
                    IF(ldecal_bc==.FALSE.) THEN
                       call calfoeamtabvois(xp,  vp,  fp, ielat, iwmax, ityp)
                    ELSE IF (ldecal_bc==.TRUE.) THEN !*!
                       call calfo_decalage(xp,  vp,  fp, ielat, iwmax, ityp)
                    END IF
		 end if
              else
                 call calfoeamcel
              endif
              potist=potist+potiseam

           end select
        end if
     end if
  end do

  ! !!! le cas parallele n'est pas pris en compte !!!

  if (lTberendsen) call calfoberend(xp,vp,fp,ityp)


  ! calcul de sigtot
  if (dmtype.ne.4.) then
     sigkine=0.
     do ilocal = 1, im
        sigkine(1:3,1) = sigkine(1:3,1) + &
             cm(ityp(ilocal))*vp(1:3,ilocal)*vp(1,ilocal)
        sigkine(1:3,2) = sigkine(1:3,2) + &
             cm(ityp(ilocal))*vp(1:3,ilocal)*vp(2,ilocal)
        sigkine(1:3,3) = sigkine(1:3,3) + &
             cm(ityp(ilocal))*vp(1:3,ilocal)*vp(3,ilocal)
     end do
     sigkine(1:3,1:3) = sigkine(1:3,1:3)/volu




#if(PARA)

     !  call MPI_ALLREDUCE(sig,sig_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     !  sig=sig_tot
     call MPI_ALLREDUCE(sigkine,sigkine_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     sigkine=sigkine_tot

#endif

     sigtot = sigkine+sig
  end if

  !      if (ldislo) call forcedislo(fp)


  if (.not.parallele) then
     if(lcontr) call contr (xp,vp,fp,ityp)
     if(ltranche) fp(:,imd+1:im)=0.0
  end if

  if (lFrozen.EQV..true.) then
     WHERE (Frozen(:,1:im)) fp(:,1:im)=0.d0
  endif


  if (ldesinteg) then
     fptot=0
     do i=1,im
        fptot(:)=fptot+fp(:,i)
     enddo
#if(PARA)
     call MPI_ALLREDUCE(fptot,fptot_tot,3,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     fptot=fptot_tot
#endif
     fptot=fptot/im_glob
     do i=1,im
        fp(:,i)=fp(:,i)-fptot(:)
     enddo
  endif

  !If (ibound==1 .OR. ibound==2 .OR. ibound==3) call surf_calc	!*!
  If (ibound == 1)                             call strain_bc	!*!
  If (ibound == 2 .OR. ibound == 3)            call stress_bc	!*!

  !stop

  if(ibrake.gt.0) then
     do i=1,im

        vn= vp(1,i)**2+vp(2,i)**2+vp(3,i)**2
        if (vn.ne.0) then
           vn=sqrt(vn)
           v1=elstopforce(ityp(i),1,1)
           nv1=1+INT(vn/v1)
           if (nv1.gt.ngrdel) then
              write(6,*)'elstop velocity > 49, rebuild elstop.in'
              stop
           end if
           f1=elstopforce(ityp(i),2,nv1)-(elstopforce(ityp(i),2,nv1)-elstopforce(ityp(i),2,nv1-1))*(nv1-vn/v1)
           !                 write (6,*)'felstop',f1,vn, vn/v1
           if (f1.le.0) then
              write(6,*)'f1<0 ?', f1
              stop
           end if
           do ic=1,3
              fp(ic,i)=fp(ic,i)-vp(ic,i)*f1/vn
              Elosselec=Elosselec+(vp(ic,i)*f1/vn)*(xp(ic,i)-xpp(ic,i))*erg2ev
              if (i==iko)then 
!                 write (6,*)'felstop',f1,vn,vp(ic,i)*f1/vn,fp(ic,i)                        
!                write(6,*)'elfp',fp(ic,i)
                 Elosselec1=Elosselec1+(vp(ic,i)*f1/vn)*(xp(ic,i)-xpp(ic,i))*erg2ev
              end if
           end do
        end if
     end do
#if(PARA)
     call MPI_ALLREDUCE(elosselectot,elosselec,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     elosselec=elosselectot
     call MPI_ALLREDUCE(elosselec1tot,elosselec1,1,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
     elosselec1=elosselec1tot
#endif



  end if











  return
end subroutine calfo


