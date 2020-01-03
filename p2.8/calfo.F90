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
  use eloss, only : calceloss,ibrake !, tcelec,ecelec,ibrake,elstopforce,elosselectot,elosselectot1,elosselec1,ngrdel,elosselec
  use elec_cell, only :i2t
#if(PARA)
  use mod_mpi
#endif

#ifdef LAMMPS_VERSION

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
  integer :: i,ilocal,ipot,ic
  !  real(double)::vn,v1,f1,ekin
  !  integer::nv1,koo
  logical:: test_sigma

#if(PARA)
  real(double), dimension(3,3) :: sig_tot,sigkine_tot
  real(double),dimension (3):: fptot_tot



#endif
  !   if (rang==0) write(6,*) 'ldemintab',ldemitab

  interface 

     subroutine calcforce_lammps2
     end subroutine calcforce_lammps2

  end interface

  
  potist=0.
  potis1=0. ; potis2=0.; potis3=0.; potis0=0. ; potcp=0.; potisP=0.
  potisTersoff=0.; potiszbl=0
  potisrep=0.; potisglue=0.; potiseam=0.

  test_sigma=(mod(it,itesigma)==0)
  if (test_sigma) then
     sig(:,:)=0.d0 ; if (ltpcel.EQV..true.) sigc=0
     if(lSigat) sigat(:,:,:)=0. ; 
     if (lsigtyp) then
        sigtyp=0. ; sigtyptyp=0.
#if(PARA)
        sigtyp_loc=0.;     sigtyptyp_loc=0.
#endif 
     end if
  end if


  fp(:,:) = zero
  jq=0.0 
  if (associated(eatom)) eatom(:)=0

#ifdef LAMMPS_VERSION
  if ((ipotentiel==-10).or.(ipotentiel==-11)) then
     do i=1,im
        posa(i)=xp(1,i)/A2cm
        posa(im+i)=xp(2,i)/A2cm
        posa(2*im+i)=xp(3,i)/A2cm
     end do


     boxl(1)=at(1,1)/A2cm
     boxl(2)=at(2,2)/A2cm
     boxl(3)=at(3,3)/A2cm

     call calcforce_lammps2

  else
#endif  



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
                 if (ltabvois) then 
                    call calfojuli(xp,  vp, fp, ielat, iwmax, ityp)
                 else
                    call calfojulicel
                 end if
              case(13,14,15)
                 if (ltabvois) then
                    ! !!! le cas parallele n'est pas pris en compte !!!
                    if (.not.parallele) call force_tersoff (xp,  vp, fp, iwmax, ityp)
                 else
                    call force_tersoff_cel
                 endif
                 potist=potist+potisTersoff+potiszbl
              case (10,11)
                 if (ltabvois) then
                    ! !!! le cas parallele n'est pas pris en compte !!!
                    if (.not.parallele) then
                       IF(ldecal_bc.EQV..FALSE.) THEN
                          !write(*,*) 'NDM eam calfo1', xp(1,1)
                          call calfoeamtabvois(xp,  vp,  fp, ielat, iwmax, ityp)
                          !write(*,*) 'NDM eam calfo2', fp(1,1), maxval(fp)
                       ELSE IF (ldecal_bc.EQV..TRUE.) THEN !*!
                          call calfo_decalage(xp,  vp,  fp, ielat, iwmax, ityp)
                       END IF
                    end if
                 else
                    call calfoeamcel
                 endif
                 potist=potist+potiseam
#if(ML)
              case (20)
                 !write(*,*) 'NDM ml calfo1', xp(1,1)
                 call md_calfo_ml
                 !write(*,*) 'NDM ml calfo2', xp(1,1), fp(1,1)
                 !stop 'ndm'
#endif          
              end select
           end if
        end if
     end do

#ifdef LAMMPS_VERSION
  endif
#endif  

  ! !!! le cas parallele n'est pas pris en compte !!!
  if (lTberendsen) call calfoberend(xp,vp,fp,ityp)


  ! calcul de sigtot
  if (dmtype.ne.4.) then
     if (test_sigma) then                   
        sigkine=0.
        do ilocal = 1, im
           sigkine(1:3,1) = sigkine(1:3,1) + &
                cm(ityp(ilocal))*vp(1:3,ilocal)*vp(1,ilocal)
           sigkine(1:3,2) = sigkine(1:3,2) + &
                cm(ityp(ilocal))*vp(1:3,ilocal)*vp(2,ilocal)
           sigkine(1:3,3) = sigkine(1:3,3) + &
                cm(ityp(ilocal))*vp(1:3,ilocal)*vp(3,ilocal)
           if (lsigat) then 
              sigat(1:3,1,ilocal) = sigat(1:3,1,ilocal) +  cm(ityp(ilocal))*vp(1:3,ilocal)*vp(1,ilocal)
              sigat(1:3,2,ilocal) = sigat(1:3,2,ilocal) +  cm(ityp(ilocal))*vp(1:3,ilocal)*vp(2,ilocal)
              sigat(1:3,3,ilocal) = sigat(1:3,3,ilocal) +  cm(ityp(ilocal))*vp(1:3,ilocal)*vp(3,ilocal)
           end if
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
  if (l2t)then
     if (i2t==1)  call calceloss
  else
     if(ibrake.gt.0) call calceloss
  end if




  !do i=1,im
  !   write(96,'(2I3,6G15.7)')i,ityp(i),xp(1,i),xp(2,i),xp(3,i),fp(1,i),fp(2,i),fp(3,i)
  !end do
  !stop
  !write(*,*) 'NDM calfo end debug', xp(1,1), fp(1,1) 

  return
end subroutine calfo


