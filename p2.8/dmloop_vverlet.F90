module dmloop_vverlet_mod
  USE calfo_mod,only: calfo
  USE analyse_mod,only: analyse
  USE controle_mod,only: controle
  USE dyn_vverlet_mod,only: dyn_vverlet
  USE calctemp_mod,only: calctemp
  USE sauvegarde_mod,only: sauvegarde
  USE sauveposition_mod,only: sauveposition
  USE sauveforce_mod,only: sauveforce
  USE correl_mod,only: correlvp
  USE atomconfig
  use var_pot,only:ntyp
  USE gen_com_m, ONLY: itesauvforce,itesauvposition,lcorrelvp,at,ecyl,ev2erg,im,lgc,rang,rayonc,&
       &tstep,vdc,pc,vdc,itdes,itesauv,itesigma,ldesinteg,lsigat,lsigtyp,ltpcel,sigat,sigc,sigtyptyp,sigtyp,noxyz
  implicit none 
contains
  ! boucle de DM pour velocity Verlet
  ! ************************************************

  subroutine dmloop_vverlet
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE Parrinello_Rahman
    USE tab_imm_m
    USE suivinonpbc


#ifdef PARA
    USE mod_para
#endif
    implicit none
    character :: extension*2
    integer::lenfn2,i
    integer::ilocal
    real(double) sigkine_tot(3,3)
    real(double) :: temptyp(ntyp)
#ifdef PARA
    ! declarations supplementaires pour MPI
    real(double), dimension(3,3,noxyz) :: sigc_tot

#endif

    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    ! MPI
    type(atom_config_d)::atdml

    if (rang==0) write (6, *) '***** PREMIERE ITERATION  ****'
#ifdef PARA
    temps_para=0.
#endif

    if (lsuivinonpbc) call init_suivinonpbc()
    !  write(6,*)'RG i ',rang,it
    ! Appel de la routine generale des forces
!    call calfo 
    call ndm2config(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob=num_at_glob,ltabvois=ltabvois,&
         &iwmax=iwmax,indi=indi,nvois=nvois,vp=vp,xpp=xpp)
    CALL CalFo(sig,potist,atdml) !(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
    call config2ndm(atdml,im,imm,xp,fp,ityp,ielat,num_at_glob,ltabvois,iwmax=iwmax,indi=indi,vp=vp,xpp=xpp)



    !  call analyse
    call calctemp (temptyp) 
1   continue
    it = it+1


    if (ldesinteg)itdes=itdes+1


    call dyn_vverlet
    ! les positions et les vitesses sont synchrones en ce point ; les atomes sont bien r�partis en cellules
    !  write(6,*)'RG i ',rang,it
    ! calcul de sigtot
    !if (lpr==.false.) then



    sigkine=0.
    do ilocal = 1, im
       sigkine(1:3,1) = sigkine(1:3,1) + &
            cm(ityp(ilocal))*vp(1:3,ilocal)*vp(1,ilocal)
       sigkine(1:3,2) = sigkine(1:3,2) + &
            cm(ityp(ilocal))*vp(1:3,ilocal)*vp(2,ilocal)
       sigkine(1:3,3) = sigkine(1:3,3) + &
            cm(ityp(ilocal))*vp(1:3,ilocal)*vp(3,ilocal)
       if (lsigat) then 
          sigat(1:3,1,ilocal) = sigat(1:3,1,ilocal) +  cm(ityp(ilocal))*vp(1:3,ilocal)*vp(1,ilocal)
          sigat(1:3,2,ilocal) = sigat(1:3,2,ilocal) +  cm(ityp(ilocal))*vp(1:3,ilocal)*vp(2,ilocal)
          sigat(1:3,3,ilocal) = sigat(1:3,3,ilocal) +  cm(ityp(ilocal))*vp(1:3,ilocal)*vp(3,ilocal)
       end if
       if ((mod(it,itesigma)==0).and.(lTPcel.EQV..true.)) then
          sigc(1:3,1,ielat(ilocal)) = sigc(1:3,1,ielat(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(1,ilocal)*noxyz/volu
          sigc(1:3,2,ielat(ilocal)) = sigc(1:3,2,ielat(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(2,ilocal)*noxyz/volu
          sigc(1:3,3,ielat(ilocal)) = sigc(1:3,3,ielat(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(3,ilocal)*noxyz/volu
       end if
       if ((mod(it,itesigma)==0).and.(lsigtyp.EQV..true.)) then
          sigtyp(1:3,1,ityp(ilocal)) = sigtyp(1:3,1,ityp(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(1,ilocal)
          sigtyp(1:3,2,ityp(ilocal)) = sigtyp(1:3,2,ityp(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(2,ilocal)
          sigtyp(1:3,3,ityp(ilocal)) = sigtyp(1:3,3,ityp(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(3,ilocal)

          sigtyptyp(1:3,1,ityp(ilocal),ityp(ilocal)) = sigtyptyp(1:3,1,ityp(ilocal),ityp(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(1,ilocal)
          sigtyptyp(1:3,2,ityp(ilocal),ityp(ilocal)) = sigtyptyp(1:3,2,ityp(ilocal),ityp(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(2,ilocal)
          sigtyptyp(1:3,3,ityp(ilocal),ityp(ilocal)) = sigtyptyp(1:3,3,ityp(ilocal),ityp(ilocal)) + &
               cm(ityp(ilocal))*vp(1:3,ilocal)*vp(3,ilocal)
       end if
    end do
    sigkine(1:3,1:3) = sigkine(1:3,1:3)/volu
#ifdef PARA

    !  call MPI_ALLREDUCE(sig,sig_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    !  sig=sig_tot
    call MPI_ALLREDUCE(sigkine,sigkine_tot,9,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
    sigkine=sigkine_tot
    if (lsigtyp.and.mod(it,itesigma) == 0) then
       call MPI_ALLREDUCE(sigtyp,sigtyp_loc,9*ntyp,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
       sigtyp=sigtyp_loc
       call MPI_ALLREDUCE(sigtyptyp,sigtyptyp_loc,9*ntyp*ntyp,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
       sigtyptyp=sigtyptyp_loc

    end if
    if (allocated(sigc)) then
       call MPI_ALLREDUCE(sigc,      sigc_tot,      9*noxyz,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_WORLD,ierr)
       sigc=sigc_tot
    end if
#endif
    sigtot = sigkine+sig

    call analyse 
    if (lcorrelvp) call correlvp(xp,xpp,vp,ax,fp,ityp)
    ! MPI
    !     write(6,*)'analyse -> sauvegarde'XS
    if (itesauv.GT.0) then
       if (mod(it,itesauv)==0) call sauvegarde
    endif

    !     write(6,*)'analyse -> sauveposition'
    if (itesauvposition.GT.0) then
       if (mod(it,itesauvposition)==0) then
          call sauveposition (it)
          if (lsuivinonpbc) then
             call reset_suivinonpbc
             call sauvepositionnonpbc (it)
          end if
       end if
    endif
    if (itesauvforce.GT.0) then
       if (mod(it,itesauvforce)==0) call sauveforce (it)
    end if
    !     write(6,*)'sauvposition -> control'

    call controle
    !     write(6,*)' controle ->'
    go to 1

    return
  end subroutine dmloop_vverlet
end module dmloop_vverlet_mod
