module babar_def_mod
!NDM is © 2021, Jean-Paul Crocombette, CEA Saclay, SRMP
!NDM is published and distributed under the Academic Software License v1.0 (ASL).
!NDM is distributed in the hope that it will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the ASL for more details.
!You should have received a copy of the ASL along with this program; if not, write to jpcrocombette@cea.fr. It is also published at https://github.com/jpcroc/NDM/blob/ndm2025/LICENSE.md.
!You may contact the original licensor at jpcrocombette@cea.fr.



  USE T_kind_param_m, ONLY:  double
  use gen_com_m,only:bk,uwrt
  use mpi
  USE Tpara,only:NDM_MPI_REAL_DOUBLE,ierr,mpi_communicator,endmpi

  implicit none

  type babar_config
     real(double)::energie
     real(double)::beta,temp
     integer:: indice ! donne beta et temp !change correspond à l'id de rex_exch
     integer:: itbb,rank ! donne le proc (rank) et le numéro du calcul pour ce proc ! fixe  sur le calcul
     integer:: ictot ! donne le icalc correspondant du babartot 
   contains
     procedure,pass:: setbetatemp
     procedure,pass:: send2proc
     procedure,pass:: recv
     
  end type babar_config
!!$
!!$  type,extends(babar_config)::babar_conf_tot
!!$     integer:: itbb,rank,indice
!!$  end type babar_conf_tot
!!$  
!!$  type,extends(babar_config)::babar_conf_loc
!!$     integer:: ictot
!!$  end type babar_conf_loc


   ! rank to replica ID array
  INTEGER, DIMENSION(:), ALLOCATABLE :: exchange_accepted,exchange_attempted,id2calc(:)
  real(double), ALLOCATABLE :: id_beta(:),replica_energie(:),id_temp(:)
  real(double)::betamin,betamax,delta_beta,delta_temp

  
contains
  subroutine send2proc(bbc,rgcib,mpic)
    type(mpi_communicator),intent(in)::mpic
    class(babar_config)::bbc
    integer,intent(in)::rgcib
    real(double)::rbuffer(3)
    integer::ibuffer(3)
    rbuffer(1)=bbc%energie
    rbuffer(2)=bbc%beta
    rbuffer(3)=bbc%temp
    ibuffer(1)=bbc%itbb
    ibuffer(2)=bbc%rank
    ibuffer(3)=bbc%indice
    if (bbc%rank.ne.rgcib) then
       write(uwrt,*)' PB babar send',rgcib,bbc%rank
    end if
    call mpic%send(ibuffer,rgcib,411)
    call mpic%send(rbuffer,rgcib,412)
  end subroutine send2proc
  
  subroutine recv(bbc,rgem,mpic)
    type(mpi_communicator),intent(in)::mpic
    class(babar_config)::bbc
    integer,intent(in)::rgem
    real(double)::rbuffer(3)
    integer::ibuffer(3)
    call mpic%recv(ibuffer,rgem,411)
    call mpic%recv(rbuffer,rgem,412)
    bbc%energie=    rbuffer(1)
    bbc%beta=    rbuffer(2)
    bbc%temp=    rbuffer(3)
    bbc%itbb=    ibuffer(1)
    bbc%rank=    ibuffer(2)
    bbc%indice=    ibuffer(3)
  end subroutine recv
    
  subroutine setbetatemp(bbt,temp,beta)
    class(babar_config)::bbt
    real(double),optional::temp,beta
    if (present(temp).and.present(beta)) then
       write(6,*)'specification temp et beta ==> stop'
       stop
    end if
    if (.not.present(temp).and.(.not.present(beta))) then
       write(6,*)'specification ni de temp ni de beta ==> stop'
       stop
    end if
    if (present(temp)) then
       bbt%temp=temp
       bbt%beta=1/(bk*temp)
    else
       bbt%temp=1/(bk*beta)
       bbt%beta=beta
    end if
  end subroutine setbetatemp
end module babar_def_mod
