module ForceMatrix_mod
  USE arret_ndm_mod,only:arret_ndm
  USE gen_com_m,only:  lperiod,lenfnam,lspaceNDM,rang,firsttime_lammps,erg2ev,fnam,fnamcout,imm_glob,fnam,lenfnam,iteration
  USE atomconfig,only:atom_config
  USE cellconfig, only:cell_config, caltabtC
  USE var_pot,only:ntyp,cm,gamlt
  USE calfo_mod,only: calfo
  USE T_kind_param_m, ONLY:  double
  USE cryst_to_cart_mod, ONLY: cryst_to_cart
  USE caltabi_mod,only: caltabi
  USE boxconfig,only:box_config,periodbox
  use paraconfig,only:para_config,commconstr,initparapuresp
#ifdef PARA
  use Tpara,only:grp_world,nprocs,myidsp,MPI_COMM_space,nprocspace,ierr,mpi_comm_world,&
       &NDM_MPI_REAL_DOUBLE,para_space_config,status,comm_space,mpi_world
  use mod_para,only:maj_atomes_frt_ftm
  USE init_vois_mod,only: init_voisinage
#else
  use Tpara,only:myidsp,nprocspace,para_space_config,nprocs
#endif
  use read_val,only:rvois,ltabvois
  use var_pot,only:ipotentiel,rumax
  USE parautils,only:initloc,pointer_caltabt_calfo
#ifdef LAMMPS_VERSION
  use vars_lammps
  use lammps_util_mod,only:init_lammps
#endif  
  use config2data_mod,only:config2data
  USE constrconf_mod,only:read_cin
  USE parautils,only:driver_caltabt_DM

  
  implicit none

  type(para_config),target::paraFM
  type(para_space_config)::pscFM
  logical::lbigmaster,lmaster,lparaFM
  logical::lwritefreq
  integer::nparaFM !nombre calculs de forces en parallele
  integer::ndecal
  real(double)::decal
  class(atom_config),pointer::atfmloc
  type(cell_config),pointer::celfmloc
  type(cell_config),target:: cellcible ! ne sert qu'à faire pointer cellnebloc sur quelquechose
  type(atom_config),target::atcible
  
contains 


  subroutine calcFM (atfm,celfm,boxfm)
    type(atom_config)::atfm
    type(box_config)::boxfm
    type(cell_config)::celfm

    integer::ideb,ifin,natl,nq,i,i1,i2,ic2,ic1,im,it1,it2,im3,ic,idecal
    real(double)::eig
    real(double),allocatable::FMat(:,:)
    real(double),allocatable::Fpzero(:,:)
    real(double),allocatable::FpMat(:,:,:,:,:)
    real(double)::fmi,potist,sig(3,3)
    real(double),allocatable::eigval(:),work(:)
    integer::info,nwork,nskip
    character*80::fnamfreqout

#ifdef PARA    
    celfmloc=>cellcible
    atfmloc=>atcible
#endif    
    
    
    im3=3*atfm%im
    im=atfm%im
    allocate(FMat(im3,im3))
    allocate(Fpzero(3,im))
    allocate(FpMat(-ndecal:ndecal,3,im,3,im))
    FMat=0.
    FpMat=0.
    Fpzero=0.0
    nq=int(atfm%im/nparaFM)
    ideb=paraFM%image*nq+1
    if (paraFM%image==paraFM%nimage-1) then 
       ifin=atfm%im
    else
       ifin=(1+paraFM%image)*nq
    end if
   
    call initloc(atfm,celfm,atfmloc,celfmloc,boxfm,paraFM,rumax,lperiod,ldistrib=.false.,psc=pscfm) !initloc contient caltabtc sur atloc
    call pointer_caltabt_calfo(sig,potist,atfm,celfm,boxfm,atfmloc,celfmloc,parafm,&
         &lperiod,lchg=.false.,psc=pscfm,lcalcvois=.false.)
    if (lmaster) then 
       Fpzero(1:3,1:im)=atfm%fp(1:3,1:im)
    end if
    iteration=1 !(empeche le recalcul de la table des voisins dans driver_caltabt_DM)
    do idecal=-ndecal,ndecal
       if (idecal==0) cycle ! pas de calcul pour décalage=0
       write(6,*)'decal rang ideb ifin',idecal,rang,ideb,ifin
          write(6,*)lmaster,rang,i,idecal,ideb,ifin
       do i=ideb,ifin
          do ic=1,3
             atfm%xp(ic,i)= atfm%xp(ic,i)+idecal*decal
             call driver_caltabt_DM(atfm,celfm,boxfm,pscfm,lperiod)
             call pointer_caltabt_calfo(sig,potist,atfm,celfm,boxfm,atfmloc,celfmloc,parafm,&
                  &lperiod,lchg=.true.,psc=pscfm,lcalcvois=.false.)
             if (lmaster)then
                Fpmat(idecal,ic,i,1:3,1:im)=atfm%fp(1:3,1:im)
             end if
             atfm%xp(ic,i)= atfm%xp(ic,i)-idecal*decal
          end do
       end do
    end do
    if (lmaster)then
       call paraFM%mpi_master%sum(Fpmat)
       do idecal=-ndecal,ndecal
          do i=1,im
             do ic=1,3
                Fpmat(idecal,ic,i,1:3,1:im)=Fpmat(idecal,ic,i,1:3,1:im)-Fpzero(1:3,1:im)
             end do
          end do
       end do
    end if
    if (rang==0) then
       do i1=1,im
          do i2=1,im
             do ic1=1,3
                do ic2=1,3
                   fmi=0.
                   do idecal=-ndecal,ndecal
                      if (idecal==0) cycle
                      fmi=fmi+Fpmat(idecal,ic1,i1,ic2,i2)/(idecal*decal)
                   end do
                   fmi=fmi/(2*ndecal)
                   it1=ic1+(i1-1)*3
                   it2=ic2+(i2-1)*3
                   fmi=-fmi/sqrt(cm(atfm%ityp(i1))*cm(atfm%ityp(i2)))
                   Fmat(it1,it2)=fmi
                end do
             end do
          end do
       end do
    end if
    call paraFM%mpi_orig%bcast(0,Fmat)

        if (rang==0) then
#ifdef MKL

       nwork=3*im3-1
       allocate(work(nwork))
       allocate(eigval(im3))
       write (6,*)'PRE MKL'
       call DSYEV('N','L',im3,Fmat,im3,eigval,work,nwork,info)

       nskip=0
       do i=1,im3
          eig=eigval(i)
          if ((dabs(eig) < 1.d10 ).or.(eig<0))           nskip=nskip+1
       end do
       if (nskip==3) then
          write(6,*)'3 non positive frequencies OK '
       else
          write(6,*)nskip, 'non positive frequencies: saddle point ? '
          write(6,*)nskip, 'THERMO IS DUBIOUS'
       end if
       if (lwritefreq) then
          fnamfreqout = fnam(1:lenfnam)//'.freq.dat'
          open(unit=122,file= fnamfreqout, form='formatted', status='unknown')
          do i=1,im3
             eig=eigval(i)
             if ((dabs(eig) < 1.d10 ).or.(eig<0)) then
                eig=sqrt(-eig)
                write(122,*)i,'I',eig
             else
                eig=sqrt(eig)
                write(122,*)i,eig
             end if
          end do
       end if

       close(122, status='keep')             
       call thermocalc(eigval,im3)

       
#else
       write(6,*)"diagonalization works with lapack or MKL"
       write(6,*)"these libraries are NOT linked by default"
       write(6,*)"link them in Makefile.ndm_your_makefile"
       write(6,*)"and recompile with make MKL=1 ndm_your_makefile"
       call arret_ndm
       
#endif
    end if
  end subroutine calcFM

  subroutine thermocalc (eigval,nmat)
    real(double),intent(in)::eigval(:)
    integer,intent(in)::nmat

    integer::ntemp,nT
    real(double)::temperature,dtemp,hz_to_ev,temperature_to_ev,kbt,kbt2,xp2,eig

    real(double),allocatable,dimension(:)::Fmin,Smin,Fcla,Scla
    integer::jmat,nskip
    logical::lskip
    character*80::fnamthout

    ntemp=300

    allocate(Fmin(ntemp))
    allocate(Fcla(ntemp))
    allocate(Smin(ntemp))
    allocate(Scla(ntemp))
    Fmin=0;Smin=0;Fcla=0;Scla=0
    hz_to_ev=0.004135665538536d-15
    temperature_to_ev=8.6173303E-05
    

    dtemp=10.


    nskip=0
    do jmat=1,nmat
       eig=eigval(jmat)
       lskip=.false.
       if ( (eig < 0).or.(dabs(eig)<1.d10)) then
          nskip=nskip+1
          cycle
       else
          eig=sqrt(eig)

          do nT=1,ntemp
             temperature=nt*dtemp
             kBT=temperature*temperature_to_ev
             kBT2=2.0d0*kBT
             xp2=eig*hz_to_ev/kBT2
             Fmin(nT)=Fmin(nT)+kBT*DLOG(2.0d0*DSINH(xp2) )
             Smin(nT)=Smin(nT)+(xp2/DTANH(xp2)-DLOG(2.0d0*DSINH(xp2)))
             Fcla(nT)=Fcla(nT)+kBT*DLOG(2.0d0*xp2 )
             Scla(nT)=Scla(nT)-(DLOG(2.d0*xp2)+1.0d0)
!             write(6,'(I5,f12.5, 5e25.15)') jmat,temperature, Fmin(nT), Fcla(nT), Smin(nT), Scla(nT),xp2
                   
          end do
       end if  !jmat
    end do
    write(6,*)'NSKIP',nskip,nmat-nskip
    fnamthout = fnam(1:lenfnam)//'.thermo.dat'
    open(unit=123,file= fnamthout, form='formatted', status='unknown')
    do nT=1,ntemp
      temperature=nT*dtemp
      write(123,'(f12.5, 4e25.15)') temperature, Fmin(nT), Fcla(nT), Smin(nT), Scla(nT)
   end do
   return
!    close(123)
  end subroutine thermocalc

!************************************     
  subroutine init_mpi_FM

#ifdef PARA
    write(6,*)'INITMPIFM'
    if (lparaFM) then 
       if (mod(nprocs,nparaFM).ne.0) then
          write(6,*)'nprocs/nparaFM <>0 STOP'
          call MPI_FINALIZE(ierr)
          call arret_ndm
       end if
       paraFM%mpi_orig%nproc=nprocs
       paraFM%mpi_orig%rank=rang
       paraFM%nimage=nparaFM
       call MPI_COMM_DUP(MPI_COMM_WORLD,paraFM%mpi_orig%comm,ierr)
       call MPI_COMM_GROUP(paraFM%mpi_orig%comm,paraFM%mpi_orig%group,ierr)
       call commconstr(paraFM)
       myidsp=paraFM%mpi_image%rank
       call MPI_COMM_free(mpi_comm_space,ierr)
       MPI_COMM_space=paraFM%mpi_image%comm
       nprocspace=paraFM%mpi_image%nproc
       call comm_space%init(MPI_COMM_SPACE)

    else
       call initparapuresp(paraFM,rang,mpi_WORLD)
    end if
    if (rang==0) lbigmaster=.true.
    lmaster=paraFM%lmaster
    


#else

    paraFM%mpi_orig%nproc=1
    paraFM%mpi_orig%rank=0
    paraFM%mpi_image%nproc=1
    paraFM%lmaster=.true.
    lmaster=.true.
    lbigmaster=.true.
    lparaFM=.false.
    nparaFM=1

#endif
    if (rang==0) then
       write(6,*)'nparafM,nprocspace, nprocs',nparafM,nprocspace, nprocs
    end if
  end subroutine init_mpi_FM

end module ForceMatrix_mod
