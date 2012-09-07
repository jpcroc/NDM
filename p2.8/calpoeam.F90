subroutine calpoeam
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m
  use gen_com_m
  use eam
  use var_pot
  use SMjuli
  implicit none
  
  integer :: k,l,iti
  real(double) ::xsp(ngrid),ysp(ngrid),bsp(ngrid),csp(ngrid),dsp(ngrid)
  real(double):: ktor,ktorho
  real(double) :: rk,rhok,rk2,rue

  eamrep(:,:,:)=0.0
  eamrho(:,:,:)=0.0
  eamglue(:,:,:)=0.0
  rue=rue_pair(1)
  ktor=rue/ngrid
  !      write(6,*) 'rue ngrid ktor ', rue,ngrid,ktor
  !repulsion
  do l=1,npair
     write(6,*)'pair pot', l,typ_pot_pair(l)
     if (typ_pot_pair(l).ne.ipotentiel) cycle
     if (rang==0) write(6,*)'paire ',l
     do k=1,ngrid            
        rk=(k*ktor) ; rk2=rk**2
        xsp(k)=rk
        !            write(6,*)k,rk
        select case(ipotentiel)
        case(10)
           call extrapolateRep(reppair(l),rk2,Erep=ysp(k))
          !if (k.lt.10) then 
          ! write(*, '(i6,"(",2d14.7,") )",2D15.7," )")') k,xsp(k),reppair(l)%potr(k)*ev2erg, rk,ysp(k)
          !end if
          !  if (k.gt.(ngrid-10)) then
          !   write(*,*)k, xsp(k), rk
          !  end if 
          !if (k.eq.ngrid) stop

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
     !         if (l.eq.3) eamrep(:,l,:)=0.0
     if (roff1(l).le.0) cycle
     if (lu_roff_pair(l).EQV..false.)cycle
     call zieg2(eamrep,csive,ngrid,ntyp,npair,catom,roff1,roff2,lu_roff_pair)
     !re-spline
     ysp(1:ngrid)=eamrep(1,l,1:ngrid)
     call cspline (ngrid,xsp,ysp,bsp,csp,dsp)
     eamrep(1,l,1:ngrid)=ysp(1:ngrid)
     eamrep(2,l,1:ngrid)=bsp(1:ngrid)
     eamrep(3,l,1:ngrid)=csp(1:ngrid)
     eamrep(4,l,1:ngrid)=dsp(1:ngrid)
     !         if (l==3) then
     !         do k=1,ngrid 
     !            write(6,*)xsp(k),eamrep(1,l,k),eamrep(2,l,k)
     !         end do
     !      end if

  end do
  !stop


  !rho
  rhomin=0.0 ; rhomax=0.0
  select case (ipotentiel)
  case(10)
     do iti=1,ntyp
        if (typ_and_pot(iti,ipotentiel).eqv..false.) cycle
        if (rang==0) write(6,*)'type ',iti
        do k=1,ngrid
           rk=(k*ktor) ; rk2=rk**2
           xsp(k)=rk
           call extrapolateRho(rhotyp(iti),rk2,rho=ysp(k))
        end do
        rhomin=min(rhomin,minval(ysp))
        rhomax=max(rhomax,maxval(ysp))
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
        rhomin=min(rhomin,minval(ysp))
        rhomax=max(rhomax,maxval(ysp))
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

  if (rang==0) then 
     if (rhominzero.eqv..true.) then
        if(rhomin.lt.0.0)  write(6,*) '******** RHOMIN <0 ****** ??????'
        
        if(rhomin.lt.0.0)  write(6,*) 'RHOMIN MIS A ZERO'
        if(rhomin.lt.0.0)  write(6,*) '******** RHOMIN <0 ****** ??????'
     endif
     if (rhomin.le.0)rhomin=0.
     if(rang==0)       write(6,*)'Rhomin Rhomax ',rhomin,rhomax
  end if
  rhomax=rhomax*14
  if(rang==0)       write(6,*)'Rhomin Rhomax ',rhomin,rhomax
  !glue


  ktorho=(rhomax-rhomin)/ngrid
  do iti=1,ntyp
     if (typ_and_pot(iti,ipotentiel).eqv..false.) cycle
     if (rang==0) write(6,*)'type ',iti
     do k=1,ngrid
        rhok=(k*ktorho) +rhomin
        xsp(k)=rhok
        select case(ipotentiel)
        case(10)
           call extrapolateEam(embtyp(iti),rhok,Embf=ysp(k))
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
  end do

  !      if(ipotentiel.eq.12) then
  !         eamglue(:,3,:)=0.0
  !     end if

  return
end subroutine calpoeam
