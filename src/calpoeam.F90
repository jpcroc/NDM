module calpoeam_mod
  USE arret_ndm_mod,only:arret_ndm
  USE zieg2_mod,only: zieg2
  USE spline_mod,only: cspline
  use calpo_mod,only:coulombbuild
  implicit none
contains
  subroutine calpoeam
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m
    USE gen_com_m, only:uwrt,lwrt,rang,pi
    USE eam
    USE eamerco
    
    USE var_pot, ONLY:csive,ipotentiel,lprtpot,rue_pot,typ_and_pot,ngrid,catom,eamrep_d,ipo,npair,ntyp,roff2,typ_pot_pair,&
         &eamrep,roff1,lu_roff_pair,eamrho,eamglue,eamrho_d,eamglue_d,auxe,alpha,iewald,rhomax,rhomin

    USE SMjuli
    implicit none

    integer :: k,l,iti,lw,ngrp1
    real(double) ::xsp(ngrid),ysp(ngrid),bsp(ngrid),csp(ngrid),dsp(ngrid)
    real(double) ::ysp_d(ngrid),bsp_d(ngrid),csp_d(ngrid),dsp_d(ngrid)
    real(double):: ktor
    real(double),dimension(:),allocatable::ktorho
    real(double) :: rk,rhok,rk2,rue,minrho,maxrho



    !interface
    !
    !subroutine zieg2(pot, pot_d, csive,ngrid, ntyp,npair, catom, roff1, roff2,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)
    !  !-----------------------------------------------
    !  !   M o d u l e s
    !  !-----------------------------------------------
    !  USE T_kind_param_m, ONLY:  double
    !
    !  implicit none
    !  !-----------------------------------------------
    !  !   D u m m y   A r g u m e n t s
    !  !-----------------------------------------------
    !  integer, dimension(:,:), allocatable  :: ipo                      ! indice des paires d'atomes
    !  integer, allocatable:: typ_pot_pair(:) ! donne le type d'interaction de la paire
    !  integer , intent(in) :: ngrid,ipotentiel
    !  integer  :: ntyp
    !  integer  :: npair
    !  real(double) , intent(in) :: csive
    !  real(double)  :: auxe= 23.06134575D-20 
    !  real(double) , intent(inout) :: pot(4,npair,0:ngrid+1),pot_d(4,npair,0:ngrid+1)
    !  real(double)  :: catom(ntyp)
    !  real(double)  :: roff1(npair)
    !  real(double)  :: roff2(npair)
    !  logical :: lu_roff_pair(npair)
    !
    !end subroutine zieg2
    !end interface

    allocate(ktorho(ntyp))
    ngrp1=ngrid+1
    eamrep(:,:,:)=0.0
    eamrho(:,:,:)=0.0
    eamglue(:,:,:)=0.0
!    if (allocated(eamrep_d)) then
       eamrep_d(:,:,:)=0.0
       eamrho_d(:,:,:)=0.0
       eamglue_d(:,:,:)=0.0
 !   end if
    rue=rue_pot(ipotentiel)
    ktor=rue/ngrid
    !      write(uwrt,*) 'rue ngrid ktor ', rue,ngrid,ktor
    !repulsion
    do l=1,npair
!            write(uwrt,*)'pair pot', l,typ_pot_pair(l),lu_roff_pair(l)
       if (typ_pot_pair(l).ne.ipotentiel) cycle
       do k=1,ngrid            
          rk=(k*ktor) !; rk2=rk**2
          xsp(k)=rk
          !            write(uwrt,*)k,rk
          select case(ipotentiel)
          case(10)
             call extrapolateRep(reppair(l),SPreppair(l),rk,Erep=ysp(k))
             if (lforcetabulate) then
                call extrapolateRep(reppair_d(l),SPreppair_d(l),rk,Erep=ysp_d(k))
             end if
          case(11)
             call extrapolateReperco(rk,ysp(k))

          case(12)
             call extrapolateRepjl(reppairjl(l),rk2,Erep=ysp(k))
             !            write(uwrt,*)'REP k,repk ',l, k,ysp(k)
          case(16)
             call extrapolateRepCRG(rk,l,ysp(k))
          end select
       end do
       call cspline (ngrid,xsp,ysp,bsp,csp,dsp)
       eamrep(1,l,1:ngrid)=ysp(1:ngrid)
       eamrep(2,l,1:ngrid)=bsp(1:ngrid)
       eamrep(3,l,1:ngrid)=csp(1:ngrid)
       eamrep(4,l,1:ngrid)=dsp(1:ngrid)
       if (lforcetabulate) then
          call cspline (ngrid,xsp,ysp_d,bsp_d,csp_d,dsp_d)
          eamrep_d(1,l,1:ngrid)=ysp_d(1:ngrid)/A2cm
          eamrep_d(2,l,1:ngrid)=bsp_d(1:ngrid)/A2cm
          eamrep_d(3,l,1:ngrid)=csp_d(1:ngrid)/A2cm
          eamrep_d(4,l,1:ngrid)=dsp_d(1:ngrid)/A2cm
       end if
       !         if (l.eq.3) eamrep(:,l,:)=0.0
       if (lprtpot)then
          lw=320+l
          do k=1,ngrid 
             write(lw,*)xsp(k)*1d8,eamrep(1,l,k),eamrep(2,l,k)
          end do
       end if
       if (roff1(l).le.0) cycle
       if (lu_roff_pair(l).EQV..false.)cycle
       call zieg2(eamrep,eamrep_d,csive,ngrid,ntyp,npair,catom,roff1,roff2,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)
       !re-spline
       if (lprtpot)then
          lw=340+l
          do k=1,ngrid 
             write(lw,*)xsp(k)*1d8,eamrep(1,l,k),eamrep(2,l,k)
          end do
       end if
       
       ysp(1:ngrid)=eamrep(1,l,1:ngrid)
       call cspline (ngrid,xsp,ysp,bsp,csp,dsp)
       eamrep(1,l,1:ngrid)=ysp(1:ngrid)
       eamrep(2,l,1:ngrid)=bsp(1:ngrid)
       eamrep(3,l,1:ngrid)=csp(1:ngrid)
       eamrep(4,l,1:ngrid)=dsp(1:ngrid)

       if (lforcetabulate) then
          write(*,*) 'There is no implemantation for the Ziegler ON and ltforcetabulate TRUE '
          write(*,*) 'Switch OFF Ziegler or put lforcetabulate to FALSE'
          write(*,*) 'Hopefully you know what you are doing!'
          call arret_ndm
          ysp_d(1:ngrid)=eamrep_d(1,l,1:ngrid) 
          call cspline (ngrid,xsp,ysp_d,bsp_d,csp_d,dsp_d)
          eamrep_d(1,l,1:ngrid)=ysp_d(1:ngrid)
          eamrep_d(2,l,1:ngrid)=bsp_d(1:ngrid)
          eamrep_d(3,l,1:ngrid)=csp_d(1:ngrid)
          eamrep_d(4,l,1:ngrid)=dsp_d(1:ngrid)
       end if

    end do
    !stop
    if (ipotentiel==16) then 

       call coulombbuild(eamrep,ipotentiel,iewald,ngrp1)

    end if


    !rho
    minrho=1d30 ; maxrho=0.0
    select case (ipotentiel)
    case(10)
       do iti=1,ntyp
          if (typ_and_pot(iti,ipotentiel).eqv..false.) cycle
          !        if (rang==0) write(uwrt,*)'type ',iti
          do k=1,ngrid
             rk=(k*ktor) ; rk2=rk**2
             xsp(k)=rk
             call extrapolateRho(rhotyp(iti),SPrhotyp(iti),rk2,rho=ysp(k))
             if (lforcetabulate) then
                call extrapolateRho(rhotyp_d(iti),SPrhotyp_d(iti),rk2,rho=ysp_d(k))
             end if
          end do
          minrho=min(minrho,minval(ysp))
          maxrho=max(maxrho,maxval(ysp))
          call cspline (ngrid,xsp,ysp,bsp,csp,dsp)
          eamrho(1,iti,1:ngrid)=ysp(1:ngrid)
          eamrho(2,iti,1:ngrid)=bsp(1:ngrid)
          eamrho(3,iti,1:ngrid)=csp(1:ngrid)
          eamrho(4,iti,1:ngrid)=dsp(1:ngrid)

          if (lforcetabulate) then
             call cspline (ngrid,xsp,ysp_d,bsp_d,csp_d,dsp_d)
             eamrho_d(1,iti,1:ngrid)=ysp_d(1:ngrid)/A2cm
             eamrho_d(2,iti,1:ngrid)=bsp_d(1:ngrid)/A2cm
             eamrho_d(3,iti,1:ngrid)=csp_d(1:ngrid)/A2cm
             eamrho_d(4,iti,1:ngrid)=dsp_d(1:ngrid)/A2cm
          end if
          if (lprtpot)then
             lw=620+iti
             do k=1,ngrid 
                write(lw,*)xsp(k)*1d8,eamrho(1,iti,k),eamrho(2,iti,k)
             end do
          end if

          !        do k=1,ngrid 
          !           write(712,*)xsp(k),eamrho(1,iti,k),eamrho(2,iti,k)
          !        end do
!          if (rang==0) write(*,*) 'calpoeam RHO MIN: ', minrho
!          if (rang==0) write(*,*) 'calpoeam RHO MAX: ', maxrho
       end do
       do iti=1,ntyp
          if (rhomax(iti)==0.)then
!             if(rang==0)       write(uwrt,*)'calpoeam Minrho MaxrhoBB ',minrho,maxrho
             rhomax(iti)=maxrho*14 ; rhomin(iti)=minrho*14
!             rhomax(iti)=maxrho*44 ; rhomin(iti)=minrho*14
          end if
       end do
    case(11)
       do iti=1,ntyp
          if (typ_and_pot(iti,ipotentiel).eqv..false.) cycle
          do k=1,ngrid
             rk=(k*ktor) ; rk2=rk**2
             xsp(k)=rk
             call extrapolateRhoerco(rk2,ysp(k))
          end do
          minrho=min(minrho,minval(ysp))
          maxrho=max(maxrho,maxval(ysp))

          call cspline (ngrid,xsp,ysp,bsp,csp,dsp)
          eamrho(1,iti,1:ngrid)=ysp(1:ngrid)
          eamrho(2,iti,1:ngrid)=bsp(1:ngrid)
          eamrho(3,iti,1:ngrid)=csp(1:ngrid)
          eamrho(4,iti,1:ngrid)=dsp(1:ngrid)
       end do
       do iti=1,ntyp
          if (rhomax(iti)==0.)then
!             if(rang==0)       write(uwrt,*)'calpoeam Minrho MaxrhoBB ',minrho,maxrho
             rhomax(iti)=maxrho*14 ; rhomin(iti)=minrho*14
          end if
       end do
       if (rang==0)write(uwrt,*)'rhomm',rhomin,rhomax
    case(16)
       minrho=0
       do iti=1,ntyp
          if (typ_and_pot(iti,ipotentiel).eqv..false.) cycle
          do k=1,ngrid
             rk=(k*ktor) ; rk2=rk**2
             xsp(k)=rk
             call extrapolateRhoCRG(rk,iti,ysp(k))
          end do
!!$          minrho=min(minrho,minval(ysp))
!!$          maxrho=max(maxrho,maxval(ysp))
          call cspline (ngrid,xsp,ysp,bsp,csp,dsp)
          eamrho(1,iti,1:ngrid)=ysp(1:ngrid)
          eamrho(2,iti,1:ngrid)=bsp(1:ngrid)
          eamrho(3,iti,1:ngrid)=csp(1:ngrid)
          eamrho(4,iti,1:ngrid)=dsp(1:ngrid)
       if (lprtpot)then
          lw=620+iti
          do k=1,ngrid 
             write(lw,*)xsp(k)*1d8,eamrho(1,iti,k),eamrho(2,iti,k)
          end do
       end if
!       rhomax(iti)=maxrho*14;rhomin(iti)=0
    end do


    case(12) !Juli
       do l=1,npair
          do k=1,ngrid
             rk=(k*ktor) ; rk2=rk**2
             xsp(k)=rk
             call extrapolateRhojl(rhotypjl(l),rk2,rho=ysp(k))
             !               write(uwrt,*)'RHO k,rhok ', l,k,rk,ysp(k)
          end do
          minrho=min(minrho,minval(ysp))
          maxrho=max(maxrho,maxval(ysp))
          call cspline (ngrid,xsp,ysp,bsp,csp,dsp)
          eamrho(1,l,1:ngrid)=ysp(1:ngrid)
          eamrho(2,l,1:ngrid)=bsp(1:ngrid)
          eamrho(3,l,1:ngrid)=csp(1:ngrid)
          eamrho(4,l,1:ngrid)=dsp(1:ngrid)
       end do
       do iti=1,ntyp
          if (rhomax(iti)==0.)then
             !             if(rang==0)       write(uwrt,*)'calpoeam Minrho MaxrhoBB ',minrho,maxrho
             rhomax(iti)=maxrho*14 ; rhomin(iti)=minrho*14
          end if
       end do
    case default
       write (uwrt, *) rang,'Bienvenue dans le cote obscur de la force : pas de potentiel ?'
       call arret_ndm
    end select



    !glue
    select case(ipotentiel)
    case default
!       if(rang==0)       write(uwrt,*)'calpoeam Rhomin Rhomax ',rhomin,rhomax
       ktorho=(rhomax-rhomin)/ngrid
       do iti=1,ntyp
          if (typ_and_pot(iti,ipotentiel).eqv..false.) cycle
          !     if (rang==0) write(uwrt,*)'type ',iti
          do k=1,ngrid
             rhok=(k*ktorho(iti)+rhomin(iti))
             xsp(k)=rhok
             select case(ipotentiel)
             case(10)
                call extrapolateEam(embtyp(iti),SPembtyp(iti),rhok,Embf=ysp(k))
                if (lforcetabulate) then
                   call extrapolateEam(embtyp_d(iti),SPembtyp_d(iti),rhok,Embf=ysp_d(k))
                end if
             case(11)
                call extrapolateEamerco(rhok,ysp(k))
             case(12)
                call extrapolateEamjl(embtypjl(iti),rhok,Embf=ysp(k))
                !            write(uwrt,*)'GLUE k,gluek ', iti, k,ysp(k)
             end select

          end do

          call cspline (ngrid,xsp,ysp,bsp,csp,dsp)
          eamglue(1,iti,1:ngrid)=ysp(1:ngrid)
          eamglue(2,iti,1:ngrid)=bsp(1:ngrid)
          eamglue(3,iti,1:ngrid)=csp(1:ngrid)
          eamglue(4,iti,1:ngrid)=dsp(1:ngrid)
          if (lforcetabulate) then
             call cspline (ngrid,xsp,ysp_d,bsp_d,csp_d,dsp_d)
             eamglue_d(1,iti,1:ngrid)=ysp_d(1:ngrid)
             eamglue_d(2,iti,1:ngrid)=bsp_d(1:ngrid)
             eamglue_d(3,iti,1:ngrid)=csp_d(1:ngrid)
             eamglue_d(4,iti,1:ngrid)=dsp_d(1:ngrid)
          end if
          if (lprtpot)then
             lw=920+iti
             do k=1,ngrid 
                write(lw,*)xsp(k),eamglue(1,iti,k),eamglue(2,iti,k)
             end do
          end if

          !        do k=1,ngrid 
          !           write(812,*)xsp(k),eamglue(1,iti,k),eamglue(2,iti,k)
          !        end do
          if (rang==0) write(uwrt,*) 'calpoeam ktorho and the inverse: ',ktorho,1.d0/ktorho

       end do
    case(16)
       rhomin=0
       do iti=1,ntyp
          if (typ_and_pot(iti,ipotentiel).eqv..false.) cycle
         ! if (rang==0) write(uwrt,*)'type ',iti,rhomin(iti),rhomax(iti)
          ktorho(iti)=rhomax(iti)/ngrid
          do k=1,ngrid
             rhok=(k*ktorho(iti)) 
             xsp(k)=rhok
             ysp(k)=-ev2erg*crg%G(iti)*dsqrt(rhok)
          end do

          call cspline (ngrid,xsp,ysp,bsp,csp,dsp)
          eamglue(1,iti,1:ngrid)=ysp(1:ngrid)
          eamglue(2,iti,1:ngrid)=bsp(1:ngrid)
          eamglue(3,iti,1:ngrid)=csp(1:ngrid)
          eamglue(4,iti,1:ngrid)=dsp(1:ngrid)
          if (lprtpot)then
             lw=920+iti
             do k=1,ngrid 
                write(lw,*)xsp(k),eamglue(1,iti,k),eamglue(2,iti,k)
             end do
          end if
       end do
    end select

    return
  end subroutine calpoeam
end module calpoeam_mod
