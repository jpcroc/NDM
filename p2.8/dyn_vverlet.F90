module dyn_vverlet_mod
  USE calfo_mod,only: calfo
  USE calfoberend_mod,only: calfoberend 
  use var_pot,only:ntyp
  USE gen_com_m, ONLY:ilangevin,itab,dmtype,fnemd,lcalcjq,lnemd,lperiod,lpr,eatom,ltranche,bg,&
       &l2T,llangevin,lsuivinonpbc,itesigma
  USE cellconfig,only: cell_config,ndm2cellconfig,cellconfig2ndm,caltabtC


#ifdef PARA
  USE layer_mod,only: layer
  USE atomconfig,only : atom_config_d,ndm2config, config2ndm
#endif
  implicit none
contains
  ! *************************************************************
  subroutine dyn_vverlet
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

  USE tab_imm_m,only:xp,xpp,vp,fp,iwmax,ityp,ielat,num_at_glob,ax,glangv
    USE jqmod
    USE suivinonpbc
    USE elec_cell,ONLY: dynelec
#ifdef PARA
  use mpi
  USE mod_para,only:MPI_COMM_space,ierr,NDM_MPI_REAL_DOUBLE,status,nprocs,temps_debpara,temps_para,maj_atomes_frt_ftm
#endif

    USE elec_cell, ONLY:TTlangevin
    USE Parrinello_Rahman
    use calfoberend_mod,only:dynlangevin
    use period_mod,only:period
    implicit none
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: i, ia,ic,il
    real(double), dimension(ntyp) :: aux
    real(double), save :: tmoyinst, imesureT
    !-----------------------------------------------
    !real(double), external :: tempinst
    real(double)::eatommoy
#ifdef PARA
    real(double)::jq_tot(3)
#endif
    type(box_config)::boxndm
    type(atom_config_d)::atdml
    type(cell_config):: celndm

    logical::test_sigma=.false.
    !      write (*,*) 'sub dynvverlet'

    timel = timel+tstep
    aux(:ntyp) = tstep/cm(:ntyp)/2.d0


    ! step 1 First half-step velocities update, v(t) -> v(t+dt/2)



    if (dmtype==2) then
       do i = 1, im
          do ic = 1, 3
             if (vp(ic,i)*fp(ic,i)<0) then
                vp(ic,i)=0.
             end if
          end do
       end do
    end if

    if (lLangevin) then
       il=2*(ilangevin-1)+1
       call dynlangevin(im,xp,vp,fp,ityp,il,Glangv)
    elseif (l2T) then
       il=2*(ilangevin-1)+1
       call TTlangevin(xp,vp,fp,ityp,il,Glangv,num_at_glob)
    else
       DO i=1, imd
          vp(1:3,i) = vp(1:3,i) + aux(iTyp(i))*fp(1:3,i)
       END DO
    end if


    !step 2  Coordinate update, x(t)-> x(t+dt)

    DO i=1, imd
       xpp(1:3,i)=xp(1:3,i)
       xp(1:3,i) = xp(1:3,i) + tstep*vp(1:3,i)
    END DO


    !conditions periodiques
    if (lsuivinonpbc) then
       DO i=1,imd
          tmpsuivi(1:3,i)=tmpsuivi(1:3,i)+tstep*vp(1:3,i)
       END DO
    end if
    if (lperiod)  call period  (imm,xp,ax)

    ! repartition des atomes dans la nouvelle boite


    ! cell dispatching

    if (.not.lpr) then
       if (itab/=0) then
          if (mod(it,itab)==0) then
                 call ndm2cellconfig(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
             call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,i&
                  &wmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
             call caltabtC(celndm,atdml,lperiod,bg)
             call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,&
                  &xpp=xpp)
             call cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !sans doute inutile
             
          endif
       endif
    end if


#ifdef PARA
    temps_debpara=MPI_Wtime()
    ! Mise a jour des atomes (locaux/frontieres/fantomes) sur tous les processeurs
    call maj_atomes_frt_ftm
    temps_para=temps_para+MPI_Wtime()-temps_debpara
    if (ltranche) call layer
#endif


    if (l2T) then
       call dynelec
    end if

    ! Force calculation
  call ndm2boxconfig(at,bg,zl,zls2,nzl,volu,normat,boxndm)
      call ndm2cellconfig(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
  call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
       &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
  jq=0.0
  if (itesigma>0) test_sigma=(mod(it,itesigma)==0)
  CALL CalFo(sig,potist,atdml,celndm,boxndm,t_sigma=test_sigma)

      if (l2t)then
       if (i2t==1)  call calceloss (atdml%im,atdml%fp,atdml%vp,atdml%ityp,atdml%ielat,atdml%num_at_glob)
    else
       if(ibrake.gt.0) call calceloss (atdml%im,atdml%fp,atdml%vp,atdml%ityp,atdml%ielat,atdml%num_at_glob)
    end if
    if (lTberendsen) call calfoberend(atdml%im,atdml%imm,atdml%xp,atdml%vp,atdml%fp,atdml%ityp)
!  write(6,*)'dml potist ',potist,atdml%potist
    call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,xpp=xpp)
    call cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !inutile (calfo ne change pas celndm) mais laissé par sécurite
    

    if (lnemd) then
       eatommoy=0.
       do i=1,imd
          eatommoy=eatommoy+eatom(i)/float(imd)
       end do
       do i=1,imd
          fp(1,i)=fp(1,i)+(eatom(i)-eatommoy)*Fnemd
       end do
    end if




    ! Second half-step velocities update, v(t+1/2dt) -> v(t+dt)
    if (llangevin.eqv..true.) then
       il=2*(ilangevin-1)+2
       call dynlangevin(im,xp,vp,fp,ityp,il,Glangv)
    elseif (l2T) then
       il=2*(ilangevin-1)+2
       call TTlangevin(xp,vp,fp,ityp,il,Glangv,num_at_glob)
    else
       DO i=1, imd
          vp(1:3,i) = vp(1:3,i) + aux(iTyp(i))*fp(1:3,i)
       END DO
    end if

    if (allocated(eatom))  eatom(1:im)=eatom(1:im)+0.5*cm(ityp(1:im))*(vp(1,1:im)**2+vp(2,1:im)**2+vp(3,1:im)**2)

    if (lcalcjq) then
       jqp=jq ; jqk=0.0 !; expvect(:)=0.0
       do i=1,imd
          !        eatom(i)=eatom(i)+0.5*cm(ityp(i))*(vp(1,i)**2+vp(2,i)**2+vp(3,i)**2)
          expvect(:)=expvect(:)+eatom(i)*xpp(:,i)
          do ic=1,3
             jqk(ic)=jqk(ic)+eatom(i)*vp(ic,i)
             !               jq(ic)=jq(ic)+eatom(i)*vp(ic,i)
          end do
       end do
       jq=jqp+jqk

#ifdef PARA
       call MPI_ALLREDUCE(jq,jq_tot,3,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
       jq=jq_tot
#endif



       if(rang==0)write(65,'(I8,3D16.8)')it-1,jq(1),jq(2),jq(3)
       !         write(66,'(I8,3D15.6)')it-1,expvect(1),expvect(2),expvect(3)
    end if



    return
  end subroutine dyn_vverlet
end module dyn_vverlet_mod
