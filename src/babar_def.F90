module babar_def_mod
  USE T_kind_param_m, ONLY:  double
  use gen_com_m,only:bk
  implicit none
    type babar_config
     real(double)::energie,beta
     integer::indice
   contains
     procedure,pass:: temp=>temperature
     procedure,pass:: temp2beta=>setbeta
  end type babar_config

contains
  function temperature(bbt)
    class(babar_config)::bbt
    real(double)::temperature
    temperature=1/(bk*bbt%beta)
  end function temperature
  subroutine setbeta(bbt,temp)
    class(babar_config)::bbt
    real(double)::temp
    bbt%beta=1/(bk*temp)
  end subroutine setbeta
end module babar_def_mod
