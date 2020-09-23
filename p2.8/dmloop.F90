module dmloop_mod
        USE calfo_mod,only: calfo
        USE dyn_mod,only: dyn
        USE analyse_mod,only: analyse
        USE controle_mod,only: controle
        USE trempe_mod,only: trempe
        USE sauvegarde_mod,only: sauvegarde
        USE correl_mod,only: correlvp
        USE sauveforce_mod,only: sauveforce
        USE sauveposition_mod,only: sauveposition
        USE gen_com_m, ONLY:itesauvforce,itesauvposition,lcorrelvp,lfire
        USE atomconfig,only : atom_config,atom_config_d,ndm2config, config2ndm
        USE cellconfig, only:cell_config,ndm2cellconfig,cellconfig2ndm
        USE boxconfig,only:box_config,boxconfig2ndm,ndm2boxconfig
        USE eloss, ONLY : calceloss,ibrake !, tcelec,ecelec,ibrake,elstopforce,elosselectot,elosselectot1,elosselec1,ngrdel,elosselec
        USE elec_cell, ONLY :i2t       
        USE calfoberend_mod,only:calfoberend

        USE gen_com_m,only: dmtype,indi,it,itesauv,ltabvois,potist,rang,sig,nvois,&
             &nox,noy,noz,noxyz,natperc,nato,ncel,atincel,deltadist,celsize,lsigat,l2t,ltpcel,&
             &sigc,sigkine,sigat,sigtot,itesigma,ltberendsen,volu,itesigma,&
             &at,bg,zl,zls2,nzl,volu,normat

        use var_pot, only: cm! iewald,l3c,npotmax,potiseam,lpotentiel,cm,ipotentiel,potisglue,potisrep,potiseam
    
        implicit none
        contains
! ************************************************
!           Sous-programme dmloop.f
!          Version MPI du 21 fevrier 2001
! ************************************************

subroutine dmloop
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double

  USE tab_imm_m,only:xp,xpp,vp,fp,iwmax,ityp,ielat,num_at_glob,ax
  USE FireModule

  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, iti,ilocal
  REAL(double) :: fire_dt, fire_alph
  INTEGER :: fire_nstep
  !-----------------------------------------------
  !
    type(atom_config_d)::atdml
    type(cell_config):: celndm
    type(box_config)::boxndm
#ifdef PARA
    real(double), dimension(3,3) :: sig_tot,sigkine_tot

#endif
    logical:: test_sigma=.false.


    ! MPI
  if (rang==0) write (6, *) '***** PREMIERE ITERATION  VERLET STD ***'

  ! Initialization
  IF ((dmtype.EQ.2).AND.lFire) THEN
          CALL init_trempe_fire(fire_dt, fire_nstep, fire_alph)
  END IF

  !      write(6,*)'im',im
1 continue
  it = it+1
  
  write(6,*)
  write(6,*)'***** ITERATION  ****', it
  
  ! appel de la routine generale des forces
  call ndm2boxconfig(at,bg,zl,zls2,nzl,volu,normat,boxndm)
  call ndm2cellconfig(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize)
  call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
       &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
  !  call celndm%print
  if (itesigma>0)      test_sigma=(mod(it,itesigma)==0)
    if (test_sigma) then
       sig(:,:)=0.d0 ; if (ltpcel.EQV..true.) sigc=0
       if(lSigat) sigat(:,:,:)=0. ; 
    end if

  CALL CalFo(sig,potist,atdml,celndm,boxndm,t_sigma=test_sigma)
    if (l2t)then
       if (i2t==1)  call calceloss (atdml%im,atdml%fp,atdml%vp,atdml%ityp,atdml%ielat,atdml%num_at_glob)
    else
       if(ibrake.gt.0) call calceloss (atdml%im,atdml%fp,atdml%vp,atdml%ityp,atdml%ielat,atdml%num_at_glob)
    end if

       if (test_sigma) then                   
          sigkine=0.
          do ilocal = 1, atdml%im
             sigkine(1:3,1) = sigkine(1:3,1) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(1,ilocal)
             sigkine(1:3,2) = sigkine(1:3,2) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(2,ilocal)
             sigkine(1:3,3) = sigkine(1:3,3) + &
                  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(3,ilocal)
             if (lsigat) then 
                sigat(1:3,1,ilocal) = sigat(1:3,1,ilocal) +  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(1,ilocal)
                sigat(1:3,2,ilocal) = sigat(1:3,2,ilocal) +  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(2,ilocal)
                sigat(1:3,3,ilocal) = sigat(1:3,3,ilocal) +  cm(atdml%ityp(ilocal))*atdml%vp(1:3,ilocal)*atdml%vp(3,ilocal)
             end if
          end do
          sigkine(1:3,1:3) = sigkine(1:3,1:3)/volu

#ifdef PARA

          !  call MPI_ALLREDUCE(sig,sig_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
          !  sig=sig_tot
          call MPI_ALLREDUCE(sigkine,sigkine_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
          sigkine=sigkine_tot

#endif

          sigtot = sigkine+sig
       end if
    if (lTberendsen) call calfoberend(atdml%im,atdml%imm,atdml%xp,atdml%vp,atdml%fp,atdml%ityp)
  
!  write(6,*)'dml potist ',potist,atdml%potist
    call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,xpp=xpp)
    call cellconfig2ndm(celndm,noxyz,nox,noy,noz,natperc,nato,ncel,atincel,deltadist,celsize) !inutile (calfo ne change pas celndm) mais laissé par sécurite
    !    indi=indiCF

!  call calfo
  select case (dmtype)

  case(1)

     call dyn 
     if (lcorrelvp) call correlvp(xp,xpp,vp,ax,fp,ityp)

  case (2) 
          IF (lFire) THEN
                  call trempe_fire (xp, xpp, vp,  fp, ielat, iwmax, ityp, &
                        fire_dt, fire_nstep, fire_alph)
          ELSE
                  call trempe (xp, xpp, vp, fp, ielat, iwmax, ityp)
          END IF

  case default
     write (6, *) 'ne sait pas quoi faire stop'
     stop
  end select


  call analyse 
  ! MPI
  if (rang==0) then
!          write(6,*)'analyse -> sauvegarde'
     if (itesauv.GT.0) then
        if (mod(it,itesauv)==0) call sauvegarde 
     endif

!          write(6,*)'analyse -> sauveposition'
     if (itesauvposition.GT.0) then
        if (mod(it,itesauvposition)==0) call sauveposition ( it)
     endif
     if (itesauvforce.GT.0) then
        if (mod(it,itesauvforce)==0) call sauveforce ( it)
     endif
!          write(6,*)'sauvposition -> control'
  endif                                   ! fin rang=0

  call controle 


  !     write(6,*)' controle ->'


  go to 1


  return

end subroutine dmloop
end module
