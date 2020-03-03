module calpoeam_mod
        use zieg2_mod
        use spline_mod
        use arret_ndm_mod
        implicit none
        contains
subroutine calpoeam
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m
  use gen_com_m
  use eam
  use eamerco
  use var_pot
  use SMjuli
  implicit none

  integer :: k,l,iti,lw
  real(double) ::xsp(ngrid),ysp(ngrid),bsp(ngrid),csp(ngrid),dsp(ngrid)
  real(double) ::ysp_d(ngrid),bsp_d(ngrid),csp_d(ngrid),dsp_d(ngrid)
  real(double):: ktor,ktorho
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
!  integer, dimension(:,:), pointer  :: ipo                      ! indice des paires d'atomes
!  integer, pointer:: typ_pot_pair(:) ! donne le type d'interaction de la paire
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




   eamrep(:,:,:)=0.0
  eamrho(:,:,:)=0.0
  eamglue(:,:,:)=0.0
  eamrep_d(:,:,:)=0.0
  eamrho_d(:,:,:)=0.0
  eamglue_d(:,:,:)=0.0
  rue=rue_pot(ipotentiel)
  ktor=rue/ngrid
  !      write(6,*) 'rue ngrid ktor ', rue,ngrid,ktor
  !repulsion
  do l=1,npair
!     write(6,*)'pair pot', l,typ_pot_pair(l)
     if (typ_pot_pair(l).ne.ipotentiel) cycle
     if (rang==0) write(6,*)'paire ',l
     do k=1,ngrid            
        rk=(k*ktor) ; rk2=rk**2
        xsp(k)=rk
        !            write(6,*)k,rk
        select case(ipotentiel)
        case(10)
           call extrapolateRep(reppair(l),SPreppair(l),rk2,Erep=ysp(k))
           if (lforcetabulate) then
            call extrapolateRep(reppair_d(l),SPreppair_d(l),rk2,Erep=ysp_d(k))
           end if
        case(11)
           call extrapolateReperco(rk2,ysp(k))

        case(12)
           call extrapolateRepjl(reppairjl(l),rk2,Erep=ysp(k))
           !            write(6,*)'REP k,repk ',l, k,ysp(k)
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
           write(lw,*)xsp(k),eamrep(1,l,k),eamrep(2,l,k)
        end do
     end if
     if (roff1(l).le.0) cycle
     if (lu_roff_pair(l).EQV..false.)cycle
     call zieg2(eamrep,eamrep_d,csive,ngrid,ntyp,npair,catom,roff1,roff2,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)
     !re-spline
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
     stop
     ysp_d(1:ngrid)=eamrep_d(1,l,1:ngrid) 
     call cspline (ngrid,xsp,ysp_d,bsp_d,csp_d,dsp_d)
     eamrep_d(1,l,1:ngrid)=ysp_d(1:ngrid)
     eamrep_d(2,l,1:ngrid)=bsp_d(1:ngrid)
     eamrep_d(3,l,1:ngrid)=csp_d(1:ngrid)
     eamrep_d(4,l,1:ngrid)=dsp_d(1:ngrid)
    end if

  end do
  !stop


  !rho
  minrho=1d30 ; maxrho=0.0
  select case (ipotentiel)
  case(10)
     do iti=1,ntyp
        if (typ_and_pot(iti,ipotentiel).eqv..false.) cycle
!        if (rang==0) write(6,*)'type ',iti
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
           write(lw,*)xsp(k),eamrho(1,iti,k),eamrho(2,iti,k)
        end do
     end if

!        do k=1,ngrid 
!           write(712,*)xsp(k),eamrho(1,iti,k),eamrho(2,iti,k)
!        end do
        if (rang==0) write(*,*) 'calpoeam RHO MIN: ', minrho
        if (rang==0) write(*,*) 'calpoeam RHO MAX: ', maxrho
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

  case(12) !Juli
     do l=1,npair
        do k=1,ngrid
           rk=(k*ktor) ; rk2=rk**2
           xsp(k)=rk
           call extrapolateRhojl(rhotypjl(l),rk2,rho=ysp(k))
           !               write(6,*)'RHO k,rhok ', l,k,rk,ysp(k)
        end do
        minrho=min(minrho,minval(ysp))
        maxrho=max(maxrho,maxval(ysp))
        call cspline (ngrid,xsp,ysp,bsp,csp,dsp)
        eamrho(1,l,1:ngrid)=ysp(1:ngrid)
        eamrho(2,l,1:ngrid)=bsp(1:ngrid)
        eamrho(3,l,1:ngrid)=csp(1:ngrid)
        eamrho(4,l,1:ngrid)=dsp(1:ngrid)
     end do
  
  case default
     write (6, *) rang,'Bienvenue dans le cote obscur de la force : pas de potentiel ?'
     call arret_ndm
  end select
  if (rhomax==0.)then
     if(rang==0)       write(6,*)'calpoeam Minrho Maxrho ',minrho,maxrho
     rhomax=maxrho*14 ; rhomin=minrho*14
  end if


  !glue

  if(rang==0)       write(6,*)'calpoeam Rhomin Rhomax ',rhomin,rhomax
  ktorho=(rhomax-rhomin)/ngrid
  do iti=1,ntyp
     if (typ_and_pot(iti,ipotentiel).eqv..false.) cycle
!     if (rang==0) write(6,*)'type ',iti
     do k=1,ngrid
        rhok=(k*ktorho) +rhomin
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
           !            write(6,*)'GLUE k,gluek ', iti, k,ysp(k)
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
     if (rang==0) write(6,*) 'calpoeam ktorho and the inverse: ',ktorho,1.d0/ktorho

  end do

  !      if(ipotentiel.eq.12) then
  !         eamglue(:,3,:)=0.0
  !     end if

  return
end subroutine calpoeam
end module
