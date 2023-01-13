module calfo2ctabvois_mod
  USE potrep_mod,only: potrep
  USE calfocommon
  USE atomconfig,only : atom_config,atom_config_d,atom_config_e
  USE cellconfig, only : cell_config
  use boxconfig,only: box_config
  implicit none
contains
  ! **********************************************************
  subroutine calfo2ctabvois(atcf,celcf,boxcf)
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY:lcalcjq,pi,potis1
    USE var_pot, ONLY:alpha,csive,ipo,rue_pair,ipo,pot
    use vect_dist_mod,only:vect_dist
    USE jqmod
    implicit none
    class(atom_config),intent(inout)::atcf
    type(cell_config),intent(in)::celcf
    class(box_config),intent(in)::boxcf
    !-----------------------------------------------
    integer :: iw2, iti, l, iw1, i, j, itj, k, iw
    logical::linter
    real(double) :: aux, alp,  sk, phu
    real(double) ::  dr, r,partsig,deltaepot
    REAL(double), dimension(1:3) :: dxp,  gradij

    aux = 23.06134575D-20
    alp = alpha/sqrt(pi)*aux
    iw2 = 0
    do i = 1, atcf%im

       iti = atcf%ityp(i)
!       l = ipo(iti,iti)
    end do

    ! --------------------------
    !   OUVERTURE BOUCLE SUR I
    ! --------------------------
!    call cryst_to_cart (imm, xp, bg, -1)    !cart vers cryst
    do i = 1, atcf%im-1

       iti = atcf%ityp(i)
       iw1 = iw2+1
       iw2 = atcf%iwmax(i)
       do iw = iw1, iw2
          j = atcf%indi(iw)
          itj=atcf%ityp(j);l=ipo(iti,itj)
          
          call vect_dist(atcf,celcf,boxcf,i,j,VJI=dxp,lperiod=boxcf%lperiod,rum=rue_pair(l),linter=linter,dist=r)
          
          if (.not.linter) cycle
          
          gradij(:)=dxp(:)/r
          sk = r/csive
          k = sk
          ! spline
          dr = r-float(k)*csive
          deltaepot=0.5*(pot(1,l,k)+dr*(pot(2,l,k)+dr*(pot(3,l,k)+dr*pot(4,l,k))))
          phu = -1.0*(pot(2,l,k)+dr*(2.0*pot(3,l,k)+dr*(3.0*pot(4,l,k))))
          potis1 = potis1+2*deltaepot
          atcf%fp(:,i)=atcf%fp(:,i)+phu*gradij(:)
          atcf%fp(:,j)=atcf%fp(:,j)-phu*gradij(:)
          ! A commenter qd lcalcjq=false pour ne pas perdre de temps dans le test
          !ra(3)=Force de j sur i
          if (lprteat) then
             select type (atcf)
             class is (atom_config_e)
                atcf%eat(i) = atcf%eat(i)+deltaepot
                atcf%eat(j) = atcf%eat(j)+deltaepot
             end select
             !              end if
          end if
!!$          if (lcalcjq) then
!!$             jqf=0.0
!!$             eat(i) = eat(i)+deltaepot
!!$             eat(j) = eat(j)+deltaepot
!!$
!!$             do ic=1,3
!!$                jqf=jqf-0.5*(ra(ic)*(vp(ic,i)+vp(ic,j)))
!!$             end do
!!$             do ic=1,3
!!$                jq(ic)=jq(ic)-jqf*dxp(ic)
!!$             end do
!!$          end if

          if (test_sigma) then
             partsig=phu/boxcf%volu
             ! calcul de sigma contrainte
             sig(1,1) = sig(1,1)+partsig*gradij(1)*dxp(1)
             sig(1,2) = sig(1,2)+partsig*gradij(1)*dxp(2)
             sig(1,3) = sig(1,3)+partsig*gradij(1)*dxp(3)
             sig(2,1) = sig(2,1)+partsig*gradij(2)*dxp(1)
             sig(2,2) = sig(2,2)+partsig*gradij(2)*dxp(2)
             sig(2,3) = sig(2,3)+partsig*gradij(2)*dxp(3)
             sig(3,1) = sig(3,1)+partsig*gradij(3)*dxp(1)
             sig(3,2) = sig(3,2)+partsig*gradij(3)*dxp(2)
             sig(3,3) = sig(3,3)+partsig*gradij(3)*dxp(3)
             if(lsigat)then
                select type (atcf)
                class is (atom_config_e)
                   atcf%sigat(1,1,i) = atcf%sigat(1,1,i)+0.5*partsig*gradij(1)*dxp(1)
                atcf%sigat(1,2,i) = atcf%sigat(1,2,i)+0.5*partsig*gradij(1)*dxp(2)
                atcf%sigat(1,3,i) = atcf%sigat(1,3,i)+0.5*partsig*gradij(1)*dxp(3)
                atcf%sigat(2,1,i) = atcf%sigat(2,1,i)+0.5*partsig*gradij(2)*dxp(1)
                atcf%sigat(2,2,i) = atcf%sigat(2,2,i)+0.5*partsig*gradij(2)*dxp(2)
                atcf%sigat(2,3,i) = atcf%sigat(2,3,i)+0.5*partsig*gradij(2)*dxp(3)
                atcf%sigat(3,1,i) = atcf%sigat(3,1,i)+0.5*partsig*gradij(3)*dxp(1)
                atcf%sigat(3,2,i) = atcf%sigat(3,2,i)+0.5*partsig*gradij(3)*dxp(2)
                atcf%sigat(3,3,i) = atcf%sigat(3,3,i)+0.5*partsig*gradij(3)*dxp(3)
                atcf%sigat(1,1,j) = atcf%sigat(1,1,j)+0.5*partsig*gradij(1)*dxp(1)
                atcf%sigat(1,2,j) = atcf%sigat(1,2,j)+0.5*partsig*gradij(1)*dxp(2)
                atcf%sigat(1,3,j) = atcf%sigat(1,3,j)+0.5*partsig*gradij(1)*dxp(3)
                atcf%sigat(2,1,j) = atcf%sigat(2,1,j)+0.5*partsig*gradij(2)*dxp(1)
                atcf%sigat(2,2,j) = atcf%sigat(2,2,j)+0.5*partsig*gradij(2)*dxp(2)
                atcf%sigat(2,3,j) = atcf%sigat(2,3,j)+0.5*partsig*gradij(2)*dxp(3)
                atcf%sigat(3,1,j) = atcf%sigat(3,1,j)+0.5*partsig*gradij(3)*dxp(1)
                atcf%sigat(3,2,j) = atcf%sigat(3,2,j)+0.5*partsig*gradij(3)*dxp(2)
                atcf%sigat(3,3,j) = atcf%sigat(3,3,j)+0.5*partsig*gradij(3)*dxp(3)
             end select
          end if
       end if
    end do

 end do

    !        do i=1,im,100
    !           write(6,*)i,fp(1,i),fp(2,i),fp(3,i)
    !        end do


!    call cryst_to_cart (imm, xp, at, 1)     !cryst vers cart

    return
  end subroutine calfo2ctabvois
end module calfo2ctabvois_mod
