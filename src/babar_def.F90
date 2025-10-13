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

  type id2babartype
     integer:: itbb,rank
  end type id2babartype
  
  ! rank to replica ID array
  integer, dimension(:) , ALLOCATABLE :: rank2id ! , sum_rank2id
  
  ! replica ID to rank array
  type(id2babartype), dimension(:) , ALLOCATABLE :: id2babar
  INTEGER, DIMENSION(:), ALLOCATABLE :: exchange_accepted,exchange_attempted
  real(double),allocatable:: replica_betas(:)
  integer::self_rank,self_id,replica_id
  real(double)::betamin,betamax,delta_beta
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
