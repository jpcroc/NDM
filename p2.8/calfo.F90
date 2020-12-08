module calfo_mod
#ifdef ML
  USE calfo_ml_mod, ONLY : md_calfo_ml
#endif 
  USE calfoew_mod,only:calfoew
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
  USE gen_com_m, ONLY:ldecal_bc,parallele,potis0,potis2,potisp&
       &,potistersoff,potiszbl,potcp,potis1,potis3,zero,rang

  USE contrainte,only:initcontr,contr
!  USE jqmod,only:jq
!  USE strain_bc_mod,only:strain_bc
!  USE stress_bc_mod,only:stress_bc
  USE force_tersoff_mod,only:force_tersoff
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE calfocommon ! stocke des variables LOCALES sig et potist eat sigat etc.
  USE cellconfig, only : cell_config
  use boxconfig,only: box_config,ndm2boxconfig,boxconfig2ndm
!#ifdef PARA
!  use mpi
!  USE mod_para,only:MPI_COMM_space,nprocspace,NDM_MPI_real_double

!#else
!    USE mod_para,only:nprocspace
!#endif
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

  subroutine calfo (sigcf,potistcf,atcf,celcf,boxcf,t_sigma)
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    class(atom_config),intent(inout),target::atcf
    type(cell_config),intent(in),target::celcf
    type(box_config),intent(in)::boxcf
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------

!    integer::im
!    integer, intent(in) ::imm
    real(double),intent(out)::potistcf,sigcf(3,3)
    real(double),dimension(:,:),allocatable:: xp,fp,xpp
!    integer,dimension(:),allocatable::ityp,ielat
!    logical::ltabvois
!    integer,allocatable ::iwmax(:)
!    integer,allocatable:: indi(:)
    
    integer :: i,ilocal,ipot,ic
    !  real(double)::vn,v1,f1,ekin
    !  integer::nv1,koo
    logical,optional, intent(in)  ::t_sigma
    ltpcel=.false.
    test_sigma=.false.
    if (present(t_sigma))test_sigma=t_sigma
    if((test_sigma).and.(celcf%ltpcel))then
       ltpcel=.true.
       sigc=>celcf%sigc
    end if
    if(celcf%icaltabt.ne.atcf%icaltabt) then
       write (6,*)'incoherence dans icaltabt'
       stop
    end if
    
    potist=0.
    potis1=0. ; potis2=0.; potis3=0.; potis0=0. ; potcp=0.; potisP=0.
    potisTersoff=0.; potiszbl=0
    potisrep=0.; potisglue=0.; potiseam=0.


    lprteat=.false.
    lsigat=.false.
    if (test_sigma) then
       sig(:,:)=0.d0 ; if (ltpcel.EQV..true.) sigc=0
    end if
    select type(atcf)
    type is (atom_config_e)
       if (atcf%lsigat)then
          lsigat=.true.
          sigat=> atcf%sigat
          sigat=0
       end if
       if (atcf%lprteat) then
          lprteat=.true.
          eat=>atcf%eat(:)
          eat=0
       end if
       
    end select
    atcf%fp(:,:) = zero


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

       call calcforce_lammps2(atcf%im,atcf%imm,atcf%xp,atcf%ityp,atcf%fp,potist)

    else
#endif  


       do ipot=0,npotmax
          if (lpotentiel(ipot).EQV..true.) then
             ipotentiel=ipot
             if(ipotentiel.lt.10) then
                select case (ipotentiel)
                case(0,1,3,4,5,6,7)
                   if (atcf%ltabvois) then
                      call calfo2ctabvois (atcf%im,atcf%imm,atcf%xp,   atcf%fp, atcf%iwmax, atcf%ityp,atcf%indi,&
                           &boxcf%at,boxcf%bg,boxcf%volu)
                   else
                      call calfo2ccel (atcf%im,atcf%imm,atcf%xp,  atcf%fp, atcf%ityp,atcf%ielat,atcf%num_at_glob,&
                           &celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,celcf%deltadist,&
                           &boxcf%at,boxcf%bg,boxcf%volu)

                   endif


                   if (iewald.ge.1) call calfoew(atcf%im,atcf%imm,atcf%xp,atcf%fp,atcf%ityp,celcf%noxyz,boxcf%at,&
                        &boxcf%bg,boxcf%volu)

                   ! Potentiel total
                   potisP = potis0+potis1+potis2+potis3
                   potist=potist+potisP

                case(2)
                   ! !!! le cas parallele n'est pas pris en compte !!!

                   if (.not.parallele) call calfow(atcf%im,atcf%imm,atcf%xp,  atcf%fp, atcf%ielat, atcf%ityp,&
                        &celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,celcf%deltadist,&
                           &boxcf%at,boxcf%bg,boxcf%volu)

                case default
                end select

                ! !!! le cas parallele n'est pas pris en compte !!!
                if (.not.parallele.and.l3c) call calfo3c (atcf%im,atcf%imm,atcf%xp,   atcf%fp, atcf%ielat, atcf%ityp,&
                        &celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,celcf%deltadist,&
                           &boxcf%at,boxcf%bg,boxcf%volu,celcf%sigc)

                ! !!! le cas parallele n'est pas pris en compte !!!
                !potentiels EAM
             else
                select case (ipotentiel)
                case(12)
                   ! !!! le cas parallele n'est pas pris en compte !!!
                   if (atcf%ltabvois) then 
                      call calfojuli(atcf%im,atcf%imm,atcf%xp,  atcf%fp, atcf%iwmax, atcf%ityp,atcf%indi,&
                           &boxcf%at,boxcf%bg,boxcf%volu)
                   else
                      call calfojulicel(atcf%im,atcf%imm,atcf%xp,  atcf%fp, atcf%ielat, atcf%ityp,&
                      &celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,celcf%deltadist,&
                           &boxcf%at,boxcf%bg,boxcf%volu)

                   end if
                case(13,14,15)
                   if (atcf%ltabvois) then
                      ! !!! le cas parallele n'est pas pris en compte !!!
                      if (.not.parallele) call force_tersoff (atcf%im,atcf%imm,atcf%xp,   atcf%fp, atcf%iwmax, &
                           &atcf%ityp,atcf%indi,boxcf%at,boxcf%bg,boxcf%volu,boxcf%zl)

                   else
                      call force_tersoff_cel(atcf%im,atcf%imm,atcf%xp,   atcf%fp, &
                           &atcf%ielat, atcf%ityp,celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,&
                           &celcf%deltadist,boxcf%at,boxcf%bg,boxcf%volu)

                   endif
                   potist=potist+potisTersoff+potiszbl
                case (10,11)
                   if (atcf%ltabvois) then
                      ! !!! le cas parallele n'est pas pris en compte !!!
                      if (.not.parallele) then
                         IF(ldecal_bc.EQV..FALSE.) THEN
                            call calfoeamtabvois(atcf%im,atcf%imm,atcf%xp,    atcf%fp, atcf%iwmax, atcf%ityp,atcf%indi,&
                                 &boxcf%at,boxcf%bg,boxcf%volu)
                         ELSE IF (ldecal_bc.EQV..TRUE.) THEN !*!
                            call calfo_decalage(atcf%im,atcf%imm,atcf%xp,   atcf%fp, atcf%iwmax, atcf%ityp,atcf%indi,&
                                 &boxcf%at,boxcf%bg,boxcf%volu)
                         END IF
                      end if
                   else
                      call calfoeamcel(atcf%im,atcf%imm,atcf%xp,   atcf%fp, atcf%ielat, atcf%ityp,atcf%num_at_glob,&
                           &celcf%noxyz,celcf%natperc,celcf%atincel,celcf%nato,celcf%ncel,celcf%deltadist,&
                           &celcf%nox,celcf%noy,celcf%noz,boxcf%at,boxcf%bg,boxcf%volu)
     
                   endif
                   potist=potist+potiseam
#ifdef ML
                case (20)
                   call md_calfo_ml
#endif          
                end select
             end if
          end if
       end do
#ifdef LAMMPS_VERSION
    endif
#endif  
    ! !!! le cas parallele n'est pas pris en compte !!!






    ! A MODULARISER
!    if (.not.parallele) then
!       if(lcontr) call contr (atcf%xp,atcf%vp,atcf%fp,atcf%ityp)
!!       if(ltranche) atcf%fp(:,imd+1:im)=0.0
!    end if

!    if (lFrozen.EQV..true.) then
!       WHERE (Frozen(:,1:im)) atcf%fp(:,1:im)=0.d0
!    endif

    sigcf=sig;potistcf=potist
    nullify(eat);nullify(sigat)
    return
  end subroutine calfo

end module calfo_mod
