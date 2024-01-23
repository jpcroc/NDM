module calfoberend_mod
  USE arret_ndm_mod,only:arret_ndm
  use atomconfig,only:atom_config_d
  USE T_kind_param_m, ONLY:  double
  use tempinstT_mod,only:tempinstT
  USE var_pot, ONLY:gamlt,cm
  USE gen_com_m, ONLY:bk,pi,text,tstep,tautcon,text,lspaceNDM
  implicit none
contains
  subroutine calfoberend(atcf)
!#ifdef PARA
!    USE Tpara,only:COMM_space,nprocspace
!#else
!    USE Tpara,only:nprocspace
!#endif
    class(atom_config_d)::atcf
    integer :: i,ic
    real(double) :: gamb,fact,tempm1


    tempm1=tempinstT(atcf)

    !      write(6,*)'jy suis'
    gamb=1./(2.*tauTcon)
    !      write(6,*)gamb,text,tempm1

    do i=1,atcf%im
       fact=cm(atcf%ityp(i))*gamb*(Text/tempm1-1.0)
       do ic=1,3
          !            write(6,*)fp(ic,i),fact*vp(ic,i)
          atcf%fp(ic,i)=atcf%fp(ic,i)+fact*atcf%vp(ic,i)
       end do

    end do

  end subroutine calfoberend

  subroutine dynlangevin(atdml,il)
    USE gen_com_m, ONLY:
    USE var_pot, ONLY:
    use atomconfig,only:atom_config_e
    class (atom_config_e)::atdml
    integer::il
    real(double)::rga
    integer :: i,ic
    real(double) :: u1,u2
    !langevin code partir du poly de Gabriel Stolz page 84, dans une version avec expoentielle comme Manuel et Cosmin
    select case (il)
    case(1)
       !     write(6,*)'rga',rga,Gl(ic,i)*sqrt(cm(ityp(1))*bk*text*(1-rga**2))/cm(ityp(1)),vp(1,1)
       do i=1,atdml%im
          rga=exp(-gamlt(atdml%ityp(i))*tstep/2)
          do ic=1,3
             !  write(6,*)'ct',cm(ityp(1)),tstep
             call random_number(u1)
             call random_number(u2)
             atdml%Glangv(ic,i)=sqrt(-2.*log(u1))*cos(2.*pi*u2)   
             atdml%vp(ic,i) = atdml%vp(ic,i)*rga+ atdml%fp(ic,i)*tstep/(cm(atdml%ityp(i))*2)&
                  &+atdml%Glangv(ic,i)*sqrt(cm(atdml%ityp(i))*bk*text*(1-rga))/cm(atdml%ityp(i))
          end do
       end do
    case(2)

       do i=1,atdml%im
          rga=exp(-gamlt(atdml%ityp(i))*tstep/2)
          do ic=1,3
             atdml%vp(ic,i) = atdml%vp(ic,i)*rga+ atdml%fp(ic,i)*tstep/(cm(atdml%ityp(i))*2)&
                  &+atdml%Glangv(ic,i)*sqrt(cm(atdml%ityp(i))*bk*text*(1-rga))/cm(atdml%ityp(i))
          end do
       end do
    case default 
       write(6,*)'check ilangevin'
       call arret_ndm
    end select



  end subroutine dynlangevin


end module calfoberend_mod

