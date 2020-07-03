module calfo_mod
#ifdef ML
  USE calfo_ml_mod, ONLY : md_calfo_ml
#endif 
  USE calfoew_mod,only:calfoew
  USE calfoberend_mod,only:calfoberend
  USE calfo2ctabvois_mod,only:calfo2ctabvois
  USE calfo2ccel_mod,only:calfo2ccel
  USE calfo3c_mod,only:calfo3c
  USE calfow_mod,only:calfow
  USE calfo_decalage_mod,only:calfo_decalage
  USE calfoeamtabvois_mod,only:calfoeamtabvois
  USE calfoeamcel_mod,only:calfoeamcel
  USE calfojuli_mod,only:calfojuli
  USE calfojulicel_mod,only:calfojulicel
  USE force_tersoff_cel_mod,only:force_tersoff_cel
  use var_pot, only: iewald,l3c,npotmax,potiseam,lpotentiel,cm,ipotentiel,potisglue,potisrep,potiseam

  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:ibound,im_glob,lcontr,ldecal_bc,ldesinteg,lsigtyp,&
       &ltberendsen,ltranche,parallele,potis0,potis2,potisp,sigkine,sigtot,sigtyp,sigtyptyp,l2t,&
       &sigat,lsigat,eatom,volu,zero,dmtype,it,itesigma,ltpcel,potistersoff,sigc,potiszbl,&
       &potiszbl,potcp,potis1,potis3,sigkine,sigtot,sigtyp,sigtyptyp,zero,sigtyp_loc,sigtyptyp_loc

  USE contrainte,only:initcontr,contr
  USE jqmod,only:jq
  USE eloss, ONLY : calceloss,ibrake !, tcelec,ecelec,ibrake,elstopforce,elosselectot,elosselectot1,elosselec1,ngrdel,elosselec
  USE elec_cell, ONLY :i2t
  USE strain_bc_mod,only:strain_bc
  USE stress_bc_mod,only:stress_bc
  USE force_tersoff_mod,only:force_tersoff
!  USE tab_imm_m,only::
!  USE tab_imm_m
  USE atomconfig,only : atom_config_d,atom_config_d,atom_config_e
  USE calfocommon ! stocke des variables LOCALES sig et potist

  USE cellconfig, only : cell_config
#ifdef PARA
  USE mod_para
#endif
#ifdef LAMMPS_VERSION
  use lammps_util_mod,only: read_lammps,calcforce_lammps2
  use vars_lammps
#endif
  implicit none
contains
  ! ************************************************
  !           Sous-programme calfo
  !routine d'appel des routines de forces
  ! ************************************************

  subroutine calfo (sigcf,potistcf,atcf,celcf)
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    type(atom_config_d),intent(inout)::atcf
    type(cell_config),intent(in)::celcf
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------

    integer::im
!    integer, intent(in) ::imm
    real(double),intent(out)::potistcf,sigcf(3,3)
    real(double),dimension(:,:),allocatable:: xp,vp,fp,xpp
!    integer,dimension(:),allocatable::ityp,ielat
!    logical::ltabvois
!    integer,allocatable ::iwmax(:)
!    integer,allocatable:: indi(:)
    
    real(double), dimension(3) :: fptot
    integer :: i,ilocal,ipot,ic
    !  real(double)::vn,v1,f1,ekin
    !  integer::nv1,koo
    logical:: test_sigma
 
#ifdef PARA
    real(double), dimension(3,3) :: sig_tot,sigkine_tot
    real(double),dimension (3):: fptot_tot



#endif
    !   if (rang==0) write(6,*) 'ldemintab',ldemitab

!    write(6,*)'dml potist1 ',potist,atcf%potist
!    call config2ndm(atcf,im,imm,potist,sig,xp,atcf%fp,vp,xpp,ityp,ielat,ltabvois,iwmax,indi)
!    write(6,*)'dml potist 2',potist,atcf%potist
    !    write(6,*)'apres call'
    if(celcf%icaltabt.ne.atcf%icaltabt) then
       write (6,*)'incoherence dans icaltabt'
       stop
    end if
    
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
#ifdef PARA
          sigtyp_loc=0.;     sigtyptyp_loc=0.
#endif 
       end if
    end if


    atcf%fp(:,:) = zero
    jq=0.0 
    if (allocated(eatom)) eatom(:)=0

#ifdef LAMMPS_VERSION
   if ((ipotentiel==-10).or.(ipotentiel==-11)) then
!       do i=1,atcf%im
!          posa(i)=atcf%xp(1,i)/A2cm
!          posa(atcf%im+i)=atcf%xp(2,i)/A2cm
!          posa(2*atcf%im+i)=atcf%xp(3,i)/A2cm
!       end do


       !     boxl(1)=at(1,1)/A2cm
       !     boxl(2)=at(2,2)/A2cm
       !     boxl(3)=at(3,3)/A2cm

       call calcforce_lammps2(atcf%im,atcf%xp,atcf%ityp,atcf%fp,potist)

    else
#endif  


       do ipot=0,npotmax
          if (lpotentiel(ipot).EQV..true.) then
             ipotentiel=ipot
             if(ipotentiel.lt.10) then
                select case (ipotentiel)
                case(0,1,3,4,5,6,7)
                   if (atcf%ltabvois) then
                      call calfo2ctabvois (atcf%im,atcf%xp, atcf%vp,  atcf%fp, atcf%iwmax, atcf%ityp,atcf%indi)
                   else
                      call calfo2ccel (atcf%im,atcf%xp, atcf%vp,  atcf%fp, atcf%ityp,atcf%ielat,atcf%num_at_glob,&
                           &celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,celcf%deltadist)

                   endif


                   if (iewald.ge.1) call calfoew(atcf%im,atcf%xp,atcf%fp,atcf%ityp,celcf%noxyz)

                   ! Potentiel total
                   potisP = potis0+potis1+potis2+potis3
                   potist=potist+potisP

                case(2)
                   ! !!! le cas parallele n'est pas pris en compte !!!

                   if (.not.parallele) call calfow(atcf%im,atcf%xp,  atcf%vp, atcf%fp, atcf%ielat, atcf%ityp,&
                        &celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,celcf%deltadist)

                case default
                end select

                ! !!! le cas parallele n'est pas pris en compte !!!
                if (.not.parallele.and.l3c) call calfo3c (atcf%im,atcf%xp,  atcf%vp,  atcf%fp, atcf%ielat, atcf%ityp,&
                        &celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,celcf%deltadist)

                ! !!! le cas parallele n'est pas pris en compte !!!
                !potentiels EAM
             else
                select case (ipotentiel)
                case(12)
                   ! !!! le cas parallele n'est pas pris en compte !!!
                   if (atcf%ltabvois) then 
                      call calfojuli(atcf%im,atcf%xp,  atcf%vp, atcf%fp, atcf%iwmax, atcf%ityp,atcf%indi)
                   else
                      call calfojulicel(atcf%im,atcf%xp,  atcf%vp, atcf%fp, atcf%ielat, atcf%ityp,&
                      &celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,celcf%deltadist)

                   end if
                case(13,14,15)
                   if (atcf%ltabvois) then
                      ! !!! le cas parallele n'est pas pris en compte !!!
                      if (.not.parallele) call force_tersoff (atcf%im,atcf%xp,  atcf%vp, atcf%fp, atcf%iwmax, atcf%ityp,atcf%indi)
                   else
                      call force_tersoff_cel(atcf%im,atcf%xp,  atcf%vp, atcf%fp, atcf%ielat, atcf%ityp,&
                      &celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,celcf%deltadist)
                   endif
                   potist=potist+potisTersoff+potiszbl
                case (10,11)
                   if (atcf%ltabvois) then
                      ! !!! le cas parallele n'est pas pris en compte !!!
                      if (.not.parallele) then
                         IF(ldecal_bc.EQV..FALSE.) THEN
                            !write(*,*) 'NDM eam calfo1', xp(1,1)
                            call calfoeamtabvois(atcf%im,atcf%xp,  atcf%vp,  atcf%fp, atcf%iwmax, atcf%ityp)
                            !write(*,*) 'NDM eam calfo2', fp(1,1), maxval(fp)
                         ELSE IF (ldecal_bc.EQV..TRUE.) THEN !*!
                            call calfo_decalage(atcf%im,atcf%xp,  atcf%vp,  atcf%fp, atcf%iwmax, atcf%ityp,atcf%indi)
                         END IF
                      end if
                   else
                      call calfoeamcel(atcf%im,atcf%xp,  atcf%vp,  atcf%fp, atcf%ielat, atcf%ityp,atcf%num_at_glob,&
                           &celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,celcf%deltadist,&
                           &celcf%nox,celcf%noy,celcf%noz)
     
                   endif
                   potist=potist+potiseam
#ifdef ML
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

    if (lTberendsen) call calfoberend(atcf%im,atcf%xp,atcf%vp,atcf%fp,atcf%ityp)


    ! calcul de sigtot
    if (dmtype.ne.4.) then
       if (test_sigma) then                   
          sigkine=0.
          do ilocal = 1, atcf%im
             sigkine(1:3,1) = sigkine(1:3,1) + &
                  cm(atcf%ityp(ilocal))*atcf%vp(1:3,ilocal)*atcf%vp(1,ilocal)
             sigkine(1:3,2) = sigkine(1:3,2) + &
                  cm(atcf%ityp(ilocal))*atcf%vp(1:3,ilocal)*atcf%vp(2,ilocal)
             sigkine(1:3,3) = sigkine(1:3,3) + &
                  cm(atcf%ityp(ilocal))*atcf%vp(1:3,ilocal)*atcf%vp(3,ilocal)
             if (lsigat) then 
                sigat(1:3,1,ilocal) = sigat(1:3,1,ilocal) +  cm(atcf%ityp(ilocal))*atcf%vp(1:3,ilocal)*atcf%vp(1,ilocal)
                sigat(1:3,2,ilocal) = sigat(1:3,2,ilocal) +  cm(atcf%ityp(ilocal))*atcf%vp(1:3,ilocal)*atcf%vp(2,ilocal)
                sigat(1:3,3,ilocal) = sigat(1:3,3,ilocal) +  cm(atcf%ityp(ilocal))*atcf%vp(1:3,ilocal)*atcf%vp(3,ilocal)
             end if
          end do
          sigkine(1:3,1:3) = sigkine(1:3,1:3)/volu




#ifdef PARA

          !  call MPI_ALLREDUCE(sig,sig_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
          !  sig=sig_tot
          call MPI_ALLREDUCE(sigkine,sigkine_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
          sigkine=sigkine_tot

#endif

          sigtot = sigkine+sig
       end if
    end if

    !      if (ldislo) call forcedislo(atcf%fp)


    if (.not.parallele) then
       if(lcontr) call contr (atcf%xp,atcf%vp,atcf%fp,atcf%ityp)
!       if(ltranche) atcf%fp(:,imd+1:im)=0.0
    end if

!    if (lFrozen.EQV..true.) then
!       WHERE (Frozen(:,1:im)) atcf%fp(:,1:im)=0.d0
!    endif


    if (ldesinteg) then
       fptot=0
       do i=1,atcf%im
          fptot(:)=fptot+atcf%fp(:,i)
       enddo
#ifdef PARA
       call MPI_ALLREDUCE(fptot,fptot_tot,3,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
       fptot=fptot_tot
#endif
       fptot=fptot/im_glob
       do i=1,atcf%im
          atcf%fp(:,i)=atcf%fp(:,i)-fptot(:)
       enddo
    endif

    If (ibound==1 .OR. ibound==2 .OR. ibound==3) then
       write(6,*) 'strain stress BC doit être traité en dehors de calfo'
       stop
!    If (ibound == 1)                             call strain_bc	!*!
!    If (ibound == 2 .OR. ibound == 3)            call stress_bc	!*!
    end If
    !stop
    if ((l2t).or.(ibrake.gt.0)) then
       write(6,*) 'calceloss  doit être traité en dehors de calfo'
       stop
    end if
!    if (l2t)then
!       if (i2t==1)  call calceloss
!    else
!       if(ibrake.gt.0) call calceloss
!    end if
    sigcf=sig;potistcf=potist
    return
  end subroutine calfo

end module calfo_mod
