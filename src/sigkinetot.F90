module sigkinetot_mod
  USE T_kind_param_m, ONLY:  double
  USE arret_ndm_mod,only:arret_ndm

  USE atomconfig,only : atom_config_d
  USE gen_com_m, ONLY:lspaceNDM
  USE boxconfig,only:box_config_lpr,box_config
  USE var_pot, ONLY:cm
    
#ifdef PARA
  use Tpara, only:comm_space,nprocspace
#endif
  implicit none
contains
  subroutine sigkinetotMC(at_n,at_np1,box,lambda,sig,sigkine,sigtot)
    class(box_config)::box
        
    class(atom_config_d),intent(in)::at_n,at_np1
    real(double),intent(in)::lambda
    real(double),dimension(3,3)::sig,sigkine,sigtot
    integer::i,j
!    call box%print
    
    sigkine(:,:)=0.d0
    do i = 1, at_n%im
       do j = 1,3
          sigkine(1:3,j) = sigkine(1:3,j) + (1-lambda)*cm(at_n%ityp(i))*at_n%vp(1:3,i)*at_n%vp(j,i)
       enddo
    enddo
    do i = 1, at_np1%im
       do j = 1,3
          sigkine(1:3,j) = sigkine(1:3,j) + lambda*cm(at_np1%ityp(i))*at_np1%vp(1:3,i)*at_np1%vp(j,i)
       enddo
    enddo
    sigkine(1:3,1:3) =sigkine(1:3,1:3)/box%Volu
!!$#ifdef PARA
!!$
!!$    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
!!$       call comm_space%sum(sigkine)
!!$    end if
!!$#endif
    sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )       
  end subroutine sigkinetotMC

  subroutine sigkinetot(atpr,boxndm,sig,sigkine,sigtot)
      class(box_config)::boxndm
        
    class(atom_config_d),intent(in)::atpr
    real(double),dimension(3,3)::sig,sigkine,sigtot
    integer::i,j
    
    sigkine(:,:)=0.d0
    do i = 1, atpr%im
       do j = 1,3
          sigkine(1:3,j) = sigkine(1:3,j) + cm(atpr%ityp(i))*atpr%vp(1:3,i)*atpr%vp(j,i)
       enddo
    enddo
    sigkine(1:3,1:3) =sigkine(1:3,1:3)/boxndm%Volu

#ifdef PARA

    if ((nprocspace.gt.1).and.(lspaceNDM.eqv..true.)) then
       call comm_space%sum(sigkine)
    end if
#endif
    sigtot = 0.5d0*(sigkine + Transpose(sigkine) + sig + Transpose(sig) )       
  end subroutine sigkinetot
    
end module sigkinetot_mod
