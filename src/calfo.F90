module calfo_mod
#ifdef ML
!  USE NDM_ML, ONLY : calfo_ml
   use mld_interface_mod, only: mld_calfo
   !use var_pot !, only: ipotentiel, rue_pot,typ_pot_pair,npair
#endif 
   USE arret_ndm_mod,only:arret_ndm
  USE calfoew_mod,only:calfoew,calfozz
  USE calfo2ctabvois_mod,only:calfo2ctabvois
  USE calfo2ccel_mod,only:calfo2ccel
  USE calfo3c_mod,only:calfo3c
  USE calfow_mod,only:calfow
  USE calfoeamtabvois_mod,only:calfoeamtabvois
  USE calfoeamcel_mod,only:calfoeamcel
  USE calfojuli_mod,only:calfojuli
  USE calfojulicel_mod,only:calfojulicel
  USE force_tersoff_cel_mod,only:force_tersoff_cel
  use var_pot, only: iewald,l3c,npotmax,potiseam,lpotentiel,cm,ipotentiel,potisglue,potisrep,potiseam,zz,potis1, npotentiel

  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:potis2,potisp,erg2ev&
       &,potistersoff,potiszbl,potcp,potis3,zero,rang,lperiod

  USE force_tersoff_mod,only:force_tersoff
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE calfocommon ! stocke des variables LOCALES sig et potist eat sigat etc.
  USE cellconfig, only : cell_config
  use boxconfig,only: box_config
#ifdef LAMMPS_VERSION
  use lammps_util_mod,only: calcforce_lammps2 
  use vars_lammps
#endif
  use Tpara,only:para_space_config,nprocspace
  implicit none
contains
  ! ************************************************
  !           Sous-programme calfo
  !routine d'appel des routines de forces
  ! ************************************************

  subroutine calfo (sigcf,potistcf,atcf,celcf,boxcf,t_sigma,psc)
    implicit none
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    class(atom_config),intent(inout),target::atcf
    type(cell_config),intent(in),target::celcf
    !AAAAAAAAAAAAAATTTTTTTTTTTTTTTTTEEEEEEEEEEEEEEEEEEENNNNNTTTTTIIIIOOOONNN (intent(in) a rammettre)
!    type(box_config),intent(in)::boxcf
    class(box_config)::boxcf
    type(para_space_config)::psc
    real(double),intent(out)::potistcf,sigcf(3,3)

    
    integer :: ipot
    logical,optional, intent(in)  ::t_sigma
    boxcf%lperiod=lperiod
    lcalcsigc=.false.
    test_sigma=.false.
    if (present(t_sigma))test_sigma=t_sigma
    if((test_sigma).and.(celcf%ltpcel))then
       lcalcsigc=.true.
       sigc=>celcf%sigc
    end if
    

    if(celcf%icaltabt.ne.atcf%icaltabt) then
       write (6,*)'incoherence dans icaltabt calfo',celcf%icaltabt,atcf%icaltabt
       call arret_ndm
    end if
    
    potistcalfo=0.
    potis1=0. ; potis2=0.; potis3=0. ; potcp=0.; potisP=0.
    potisTersoff=0.; potiszbl=0
    potisrep=0.; potisglue=0.; potiseam=0.


    lprteat=.false.
    lsigat=.false.
    if (test_sigma) then
       sigcalfo(:,:)=0.d0 ; if (lcalcsigc.EQV..true.) sigc=0
    end if
    select type(atcf)
    class is (atom_config_e)
       if (atcf%lsigat)then
          lsigat=.true.
          atcf%sigat=0
       end if
       if (atcf%lprteat) then
          lprteat=.true.
          atcf%eat=0
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

       call calcforce_lammps2(boxcf%at,atcf%im,atcf%imm,atcf%xp,atcf%ityp,atcf%fp,potistcalfo,sigcalfo)

    else
#endif  


       do ipot=0,npotmax
          if (lpotentiel(ipot).EQV..true.) then
             ipotentiel=ipot
             if(ipotentiel.lt.10) then
                select case (ipotentiel)
                case(0,1,3,4,5,6,7,8,9)
                   if (atcf%ltabvois) then
                      call calfo2ctabvois(atcf,celcf,boxcf) !im,imm,xp,   fp,  iwmax, ityp,indi,at,bg,volu )
!!$                      call calfo2ctabvois (atcf%im,atcf%imm,atcf%xp,   atcf%fp, atcf%iwmax, atcf%ityp,atcf%indi,&
!!$                           &boxcf%at,boxcf%bg,boxcf%volu)
                   else
                      call calfo2ccel(atcf,celcf,boxcf)
                      if (test_sigma)sigcalfo=sigcalfo+sig2p ! out of calfo2ccel to accomodate ARPS
                   endif

                   ! Potentiel total
                   potisP = potis1    !  +potis2    !+potis3
                   potistcalfo=potistcalfo+potisP

                case(2)
                   ! !!! le cas parallele n'est pas pris en compte !!!

                   if (nprocspace==1) call calfow(atcf,celcf,boxcf)

                case default
                end select

                ! !!! le cas parallele n'est pas pris en compte !!!
                if ((nprocspace==1).and.l3c) call calfo3c(atcf,celcf,boxcf)

                ! !!! le cas parallele n'est pas pris en compte !!!
                !potentiels EAM
             else
                select case (ipotentiel)
                case(12)
                   ! !!! le cas parallele n'est pas pris en compte !!!
                   if (atcf%ltabvois) then 
                      call calfojuli(atcf,celcf,boxcf)

                   else
                      call calfojulicel(atcf,celcf,boxcf)

                   end if
                case(13,14,15)
                   if (atcf%ltabvois) then
                      ! !!! le cas parallele n'est pas pris en compte !!!
                      if (nprocspace==1) call force_tersoff(atcf,celcf,boxcf)
                   else
                      call force_tersoff_cel(atcf,celcf,boxcf,psc)
                   endif
                   potistcalfo=potistcalfo+potisTersoff+potiszbl
                case (10,11,16)
                   if (atcf%ltabvois) then
                      ! !!! le cas parallele n'est pas pris en compte !!!
                      if (nprocspace==1) then
                         call calfoeamtabvois(atcf,celcf,boxcf)
                      end if
                   else
                      call calfoeamcel(atcf,celcf,boxcf,psc)
                   endif
                   potistcalfo=potistcalfo+potiseam
                                      
#ifdef ML
                case (20)
                   call mld_calfo(atcf,boxcf,potistcalfo,sigcalfo,celcf) 
#endif          
                end select
             end if
          end if
       end do
       if ((iewald.gt.0).and.(iewald.ne.3))then
          call  calfoew(atcf,celcf,boxcf)
          potistcalfo=potistcalfo+potis3
          if (any(zz.ne.0)) then
             call calfozz(atcf)
             potistcalfo=potistcalfo+potis2
          end if
       end if

#ifdef LAMMPS_VERSION
    endif
#endif  

    sigcf=sigcalfo;potistcf=potistcalfo


    return
  end subroutine calfo

end module calfo_mod
