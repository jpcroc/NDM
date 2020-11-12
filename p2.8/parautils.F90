module parautils

  use paraconfig,only:para_config
#ifdef PARA
  use Tpara,only:NDM_MPI_REAL_DOUBLE
  USE mod_para,only:maj_atomes_frt_ftm
#endif
  use T_kind_param_m, ONLY:  double

  use atomconfig,only: atom_config,atom_config_d
  USE boxconfig,only:box_config
  USE cellconfig,only:cell_config,caltabtC
  use calfo_mod,only:calfo
  use constrconf_mod,only:repartition
   USE caltabi_mod,only: caltabi
  implicit none
contains
  subroutine initloc(atcomp,cellcomp,atloc,celloc,box,div,rum,lperiod)
    USE setcell,only:setcellconf
    class(atom_config_d),intent(in),target::atcomp
    type(cell_config),intent(in),target::cellcomp
    type(box_config)::box
    class(atom_config_d),pointer::atloc
    type(cell_config),pointer::celloc
    type(para_config),intent(in)::div
    real(double),intent(in)::rum
    logical::lperiod
    
    if (div%npim.gt.1) then
       call cellcomp%copy_cell(celloc)
       call repartition(atcomp,atloc,box,celloc) ! mettre les éléments de la répartition dans un type
       call setcellconf(celloc,atloc,box,atcomp%im,rum)
    else
       atloc=>atcomp
       celloc=>cellcomp
    end if
  
       call caltabtC(celloc,atloc,lperiod,box)
  end subroutine initloc
  
  subroutine pointer_caltabt_calfo(sig,potist,atcomp,cellcomp,box,atloc,celloc,div,lperiod,ltabvois,it,itetabvois,lchg)

    real(double)::sig(3,3),potist
    type(atom_config_d),target::atcomp
    type(cell_config),target::cellcomp
    type(box_config)::box
    type(para_config)::div
    class(atom_config_d),pointer::atloc
    type(cell_config),pointer::celloc
    logical,intent(in)::lperiod
    logical,optional,intent(in)::ltabvois
    integer,optional,intent(in)::itetabvois,it
    logical,optional::lchg
    
    logical::lchange=.true.
    
    if(present(lchg))lchange=lchg


#ifdef PARA
    if (lchange) then
       if (div%npim.gt.1) then
          call atcomp%atom_config%master2loc(atloc,div)
          call maj_atomes_frt_ftm(atloc,celloc)
       else
          nullify(atloc);atloc=>atcomp
          nullify(celloc);celloc=>cellcomp                 
       end if
    end if
#else
    nullify(atloc);atloc=>atcomp
    nullify(celloc);celloc=>cellcomp
       
#endif          

  
    call caltabtC(celloc,atloc,lperiod,box)
    if (present(ltabvois)) then
       if (ltabvois.and.((it==1).or.(mod(it,itetabvois)==0)))&
            &call caltabi(atloc,celloc,box)
    end if

    CALL CalFo(sig,potist,atloc,celloc,box,t_sigma=.true.) 
    
!#ifdef PARA
    if (div%npim.gt.1) then
       call atloc%vers_master(atcomp,div)
    end if
!#endif

    return
  end subroutine pointer_caltabt_calfo



    
end module parautils



