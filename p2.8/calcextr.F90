!****************************************************************
module calcextr_mod
    USE T_kind_param_m, ONLY:  double
  USE arret_ndm_mod,only:arret_ndm
  USE gen_com_m, ONLY:lperiod
  use atomconfig,only: atom_config
  USE cellconfig,only:cell_config,caltabtc
  use boxconfig,only:box_config
  use Tpara,only:para_space_config,nprocspace
  USE caltabi_mod,only: caltabi
  USE calfo_mod,only: calfo
  ! ceci est un test du calcul des forces sur un sous-ensemble des atomes
  ! ne fonctionne que pour les pots de paires
  
  
  implicit none
contains

  subroutine calfoextr(atcomp,celndm,boxndm,psc,t_sigma)
    class(atom_config),intent(inout)::atcomp
    type(box_config),intent(in)::boxndm
    type(para_space_config),intent(in)::psc
    type(cell_config),intent(in)::celndm
    logical,optional::t_sigma

    logical::tsig=.true.
    type(atom_config)::atextr
    type(cell_config)::celextr
    type(box_config)::boxextr
    real(double)::sigextr(3,3),potisextr
    integer::i
    

    if (present(t_sigma))tsig=t_sigma
    boxextr=boxndm
    
    call celndm%print (unit=445)
    call atcomp%sort(atextr)
!    call celndm%copy(celextr,boxndm,lzeroinit=.true.)
    write(6,*)'CALCEXTR',atextr%ltabvois
    call caltabtc(celextr,atextr,lperiod,boxextr,lextr=.true.,psc=psc)
!    call celextr%print (unit=444)
    if (atextr%ltabvois) call caltabi(atextr,celextr,boxndm,lextr=.true.)

    call calfo(sigextr,potisextr,atextr,celextr,boxextr,tsig,psc)
    do i=1,atextr%im
       write(345,*)i,atextr%xp(:,i)
       write(345,*)i,atextr%fp(:,i)
    end do
  end subroutine calfoextr
    
end module calcextr_mod
