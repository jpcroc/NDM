module one_calc_mod
  USE calfo_mod,only: calfo
  USE T_kind_param_m, ONLY:  double
  USE atomconfig,only : atom_config_d, atom_config_e
  USE cellconfig, only:cell_config ,caltabtC

  USE boxconfig,only:box_config,initbox
  use var_pot,only:ntyp
  USE gen_com_m, ONLY: erg2eV,rang,lspaceNDM,potist,sig,evA2dyn,itmax

  USE eloss, ONLY : calceloss,ibrake !, tcelec,ecelec,ibrake,elstopforce,elosselectot,elosselectot1,elosselec1,ngrdel,elosselec
  USE elec_cell, ONLY :i2t       
USE calfoberend_mod,only:calfoberend
use Tpara,only:para_space_config
use endrunT_mod,only:endrunT
  USE Mat_utils_mod,only:  MatInv
USE cryst_to_cart_mod,only: cryst_to_cart
 USE rasmolT_mod,only: rasmolT

  implicit none 
contains
  ! boucle de DM pour velocity Verlet
  ! ************************************************

  subroutine one_calc(atdml,celndm,boxndm,psc)

#ifdef PARA
    USE Tpara,only:COMM_space,nprocspace
#else
    USE Tpara,only:nprocspace
#endif
    USE recips_mod,only: recips,calcvol
    implicit none
    type(para_space_config)::psc
    type(box_config)::boxndm
    class(atom_config_d)::atdml
    type(cell_config):: celndm

    real(double)::unitP=1d-9
    integer::ic,i,ic2
    real(double), dimension(3,3) :: h, hDot,NAT
    real(double), dimension(3,3) :: trh, invh, invtrh, Gmat, invGmat, Gdot,forcebox
    real(double) :: invVolu,pre,x

    ! Vecteurs de la boîte et grandeurs associées à l'instant initial
    do i=1,itmax
       h(:,:)=boxndm%at(:,:)
       trh=Transpose(h)
       Gmat = MatMul(trh,h)
       CALL MatInv(Gmat,invGmat)
       call MatInv(h,invh)
       invtrh = Transpose(invh)
       boxndm%volu = calcvol(h(1:3,1),h(1:3,2),h(1:3,3))
       invVolu = 1.d0/boxndm%volu


       call rasmolT(atdml,boxndm,i,latcomp=.true.)
       CALL CalFo(sig,potist,atdml,celndm,boxndm,t_sigma=.true.,psc=psc)
       forcebox(:,:)=MatMul( sig(:,:) , invtrh(:,:) )*boxndm%volu

       write(6,*)
       write(6,*)'----------------------------------------------------'
       do ic = 1, 3
          write (6, '(I1,3(A,I1),A,3E18.10)') ic,' at (1,', ic, ') (2,', ic, &
               ') (3,', ic, ') =',boxndm%at(1:3,ic)
       end do
       do ic = 1, 3
          write (6, '(I1,3(A,I1),A,3E18.10)') ic,' invtrh(1,', ic, ') (2,', ic, &
               ') (3,', ic, ') =',invtrh(1:3,ic)
       end do
       write(6,*)'volu',boxndm%volu
       write(6,*)
       write(6,*)'Potist cgs ',potist
       write(6,*)
       do ic = 1, 3
          write (6, '(I1,3(A,I1),A,3G18.10)') ic,' sigma (1,', ic, ') (2,', ic, &
               ') (3,', ic, ') =',sig(1:3,ic)
       end do
       Pre=(sig(1,1)+sig(2,2)+sig(3,3))/3
       write(6,'(A,G15.7)')'pression',Pre*unitP

       !    write(6,*) invtrh(:,:) ,boxndm%volu
       write(6,*) 
       do ic = 1, 3
          write (6, '(I1,3(A,I1),A,3G18.10)') ic,' forcebox (1,', ic, ') (2,', ic, &
               ') (3,', ic, ') =',forcebox(1:3,ic)
       end do
       write(6,*)
!!$    write(6,*) 
!!$    write(6,*)'Potist eV ',potist*erg2ev
!!$    write(6,*)
!!$    do ic = 1, 3
!!$       write (6, '(I1,3(A,I1),A,3G18.10)') ic,' sigma (1,', ic, ') (2,', ic, &
!!$            ') (3,', ic, ') =',sig(1:3,ic)*unitP
!!$    end do
!!$
!!$    write(6,*) 
!!$    do ic = 1, 3
!!$       write (6, '(I1,3(A,I1),A,3G18.10)') ic,' forcebox(1,', ic, ') (2,', ic, &
!!$            ') (3,', ic, ') =',forcebox(1:3,ic)
!!$    end do
!       call atdml%print(caracT='x')
       call cryst_to_cart (atdml%im, atdml%xp, boxndm%bg, -1) 
!       call atdml%print(caracT='x')
       do ic = 1, 3
          do ic2=1,3
             NAT(ic2,ic)= boxndm%at(ic2,ic)+forcebox(ic2,ic)*1d-7
!             if(ic.ne.ic2) nat(ic2,ic)=0
          end do
       enddo
!       x=NAT(1,1)
!       NAt=boxndm%at
!       nat(1,1)=x
       write(6,*)'NAT',nat
       
       call initbox(boxndm,NAT)
!       write(6,*)boxndm%at


       call cryst_to_cart (atdml%im, atdml%xp, boxndm%at, 1)
!       call atdml%print(caracT='x')
       call caltabtC(celndm,atdml,.true.,boxndm)
       do ic = 1, 3
          write (6, '(I1,3(A,I1),A,3G18.10)') ic,' at (1,', ic, ') (2,', ic, &
               ') (3,', ic, ') =',boxndm%at(1:3,ic)
       end do
       write(6,*)'volu',boxndm%volu
       write(6,*)
       
    end do
  end subroutine one_calc
end module one_calc_mod
