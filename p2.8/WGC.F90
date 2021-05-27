module WGC_mod

  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:  inv_angst, lperiod, rang,itmax,leev,sig, &
       it, itesauv, itesauvposition, itesauvforce,itmax, fnam,lenfnam,fnamcout,&
       inv_angst, erg2ev, angst,fpstop,fsumstop,itetabvois, iterasmol,&
       dmtype, potist,mdcg_noise,formatsauv,lspaceNDM,latcomp,sigstop,sigext,ihbox0,unitP,lprahman
  USE sauvegardeT_mod,only: sauvegardeT
  USE endrunT_mod,only: endrunT
  USE arret_ndm_mod,only: arret_ndm
  USE initspeed_mod,only: bruit_xp
#ifdef PARA
  use Tpara,only:COMM_space,myidsp,nprocspace,para_space_config
  use mod_para,only:maj_atomes_frt_ftm
#else
  use Tpara,only:nprocspace,para_space_config
#endif
  USE atomconfig,only : atom_config
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
  use paraconfig,only:para_config
  USE parautils,only:initcomp,depeche_mode
  USE Mat_utils_mod,only:  MatInv
  USE scalebox_mod,only: scalebox
  USE boxconfig,only:box_config,periodbox,initbox
  USE recips_mod,only: recips ,calcvol
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE rasmolT_mod,only: rasmolT
  
  implicit none

  type(box_config)::boxcgmin
  type(box_config),target::boxcg
  type(atom_config),target::atcgcomp
  type(atom_config)::atcgmin
  type(cell_config),target::cellcgcomp
  type(para_space_config),target::pscCG
  class(atom_config),pointer::atcgloc
  type(cell_config),pointer::cellcgloc
  type(cell_config),target:: cellcible ! ne sert qu'à faire pointer cellnebloc sur quelquechose
  type(atom_config),target::atcible
  type(para_config),target::gcpara
  real(double), dimension(3,3) :: trh, invh, invtrh, forcebox,h,sigsym
  real(double),allocatable,dimension (:)::R,F,Rmin
  integer::N,ndir,nstep,ityprel
  real(double)::betaguess,V,betaV,betaP,beta
  integer::ncalls,nextsauv,nextmol
  logical::lvm
  real(double)::fpstop0,fpstopsig
  logical,target:: lchg
contains

  subroutine initsteep
    real(double), allocatable :: bruitmd(:,:)
    integer::i,i1,i2,ip
    allocate(bruitmd(3,1:atcgcomp%im))
    select case(ityprel)

    case(1)
       fpstop=fpstop0
       N=atcgcomp%im*3
       if (allocated (R).or.allocated(F)) then
          deallocate(R,F,Rmin)
       end if
       allocate(R(N))
       allocate(Rmin(N))
       allocate(F(N))
       if (mdcg_noise /= 0 ) then
          call bruit_xp (atcgcomp%xp,bruitmd,atcgcomp%im)
          bruitmd=bruitmd*1d-8
       else
          bruitmd=0
       end if

       R(:)=0;V=0;F(:)=0
       atcgmin=atcgcomp
       do i=1,atcgcomp%im
          i1=atcgcomp%num_at_glob(i)
          ! Variables = cartesian coordinates (in cm)
          atcgcomp%xp(1:3,i)= atcgcomp%xp(1:3,i)+bruitmd(1:3,i)
          R(3*i1-2:3*i1) = atcgcomp%xp(1:3,i)
       end do
       Rmin=R
    case(2)
       if (allocated (R).or.allocated(F)) then
          deallocate(R,F,Rmin)
       end if



       fpstop=(sigstop/unitP)*(boxcg%volu**0.6666666)
	fpstopsig=fpstop
       call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%bg, -1) 
       N=9
       allocate(R(N))
       allocate(Rmin(N))
       allocate(F(N))
       R(:)=0;V=0;F(:)=0
       boxcgmin=boxcg
       ip=0
       do i1=1,3
          do i2=1,3
             ip=ip+1
             R(ip)=boxcg%at(i1,i2)
          end do
       end do
       Rmin=R
    end select
  end subroutine initsteep

  subroutine test_conv(Nt,Ft,lover,Vt,Rt,lvm )
    integer,intent(in)::Nt
    real(double),intent(in)::Ft(Nt)
    real(double),intent(in),optional::Rt(Nt)
    real(double),optional::Vt
    logical::lover
    logical,optional::lvm
    real(double)::Vminabs=1d16

    real(double)::forctot,formax,deltaV,sigmax,fsigmax
    integer::i,i1,i2,ip
    if (present(lvm))lvm=.false.
    forctot=sqrt( SUM(Ft(:)**2))
    formax=0

    lover=.false.

    if (present(Vt)) then
       if (Vt.lt.Vminabs) then
          lvm=.true.
          deltaV=Vminabs-Vt
          if (Vminabs==1d16) deltaV=0
          Vminabs=Vt
          select case(ityprel)
          case(1)
             do i=1,atcgmin%im
                i1=atcgmin%num_at_glob(i)
                ! Variables = cartesian coordinates (in cm)
                atcgmin%xp(1:3,i)=Rt(3*i1-2:3*i1) 
             end do

          case(2)
             ip=0
             do i1=1,3
                do i2=1,3
                   ip=ip+1
                   h(i1,i2)=Rt(ip)
                end do
             end do
             call initbox(boxcgmin,h)
          end select
       end if
    end if
    select case (ityprel)
    case(2)
       Fsigmax=0;sigmax=0
       do ip=1,9
          Fsigmax=max(Fsigmax,abs(Ft(ip)))
       end do
       do i1=1,3
          do i2=1,3
             sigmax=max(sigmax,abs(unitP*sigsym(i1,i2)))
          end do
       end do
       if (fsigmax.le.fpstop) lover=.true.
       !       write(6,*)sigmax,sigstop
       !       if (sigmax.le.sigstop) lover=.true.
    case(1)
       forctot=sqrt( SUM(Ft(:)**2))
       formax=0
       do i=1,N
          formax = Max( formax,Abs(Ft(i)))
       end do

       forctot = forctot*erg2eV/angst
       formax  = formax*erg2eV/angst
!       write(6,'(A,4E20.11)')'TEST', forctot,formax
       if (fpstop>0) then   
          if (formax.le.fpstop) then
             lover=.true.
          end if
       end if
       if (fsumstop>0) then   
          if (forctot.le.fsumstop) then
             lover=.true.
          end if
       end if
    end select
    if (present(Vt)) then
       select case(ityprel)
       case(1)
          if(lvm) then
             write(6,'(I4,4E20.11,A, 2E20.11)')ncalls, Vt,Vt*erg2eV ,forctot,formax,' ****', deltaV, deltaV*erg2eV
          else
             write(6,'(I4,4E20.11)')ncalls, Vt,Vt*erg2eV ,forctot,formax
          end if
       case(2)
          if(lvm) then
             write(6,'(I4,4E20.11,A, 2E20.11)')ncalls, Vt,Vt*erg2eV ,Fsigmax, sigmax, ' ****', deltaV, deltaV*erg2eV
          else
             write(6,'(I4,4E20.11)')ncalls, Vt,Vt*erg2eV ,Fsigmax,sigmax
          end if

       end select
       if (lvm)then
          if (NCALLS.ge.nextsauv) then
             formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'.cout'
             call sauvegardeT(atcgcomp,cellcgcomp,boxcg,formatsauv,fnamcout,latcomp=.true.)
             nextsauv=NCALLS+itesauv
          end if
          if (NCALLS.ge.nextmol) then
             call rasmolT(atcgcomp,boxcg,it,latcomp=.true.)
             nextmol=NCALLS+iterasmol
          end if
       end if
    end if
    return
  end subroutine test_conv

  subroutine final_tconv(lover)
    logical,intent(out)::lover
    real(double):: forctot,sigmax,formax
    integer::i1,i2,i
    forctot=sqrt( SUM(atcgcomp%fp(:,:)**2))
    formax=0
    sigmax=0
    do i=1,atcgcomp%im
       do i2=1,3
          formax = Max( formax,Abs(atcgcomp%fp(i2,i)))
       end do
    end do
    do i=1,3
       do i2=1,3
          sigmax = Max( sigmax,Abs(unitp*sigsym(i,i2)))
       end do
    end do

    forctot = forctot*erg2eV/angst
    formax  = formax*erg2eV/angst

    write(6,*)'  FORCE MAX            FORCETOT            SIGMAX'
    write(6,'(3E20.11)')formax,forctot,sigmax
        write(6,'(A,3E20.11)')'seuils',fpstop0,fsumstop,sigstop
    write(6,*)
    lover=.false.
    if (lprahman) then
       if (fpstop.gT.0) then
          if ((formax.le.fpstop0).and.(sigmax.le.sigstop))lover=.true.
       end if
       if (fsumstop.gT.0) then
          if ((forctot.le.fsumstop).and.(sigmax.le.sigstop))lover=.true.
       end if
    else
       if (fpstop.gT.0) then
          if (formax.le.fpstop)lover=.true.
       end if
       if (fsumstop.gT.0) then
          if (forctot.le.fsumstop)lover=.true.
       end if
    end if
    return
  end subroutine final_tconv



  subroutine back2NDM( N,R,V,F,lover)
    integer,intent(in)::N
    real(double),intent(in)::F(N),R(N),V
    logical,intent(out) :: lover

    integer::i,i1,i2,ip

    select case(ityprel)
    case(1)
       if (lover) then
          potist=V
          do i=1,atcgcomp%im
             i1=atcgcomp%num_at_glob(i)             
             atcgcomp%fp(1:3,i)=F(3*i1-2:3*i1)
             atcgcomp%xp(1:3,i)=R(3*i1-2:3*i1)
          end do
       else
          atcgcomp=atcgmin
       end if
    case(2)
       if (lover) then
          ip=0
          do i1=1,3
             do i2=1,3
                ip=ip+1
                h(i1,i2)=R(ip)
             end do
          end do
          call initbox(boxcg,h)
          call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%at, 1) 
          potist=V
       else
          boxcg=boxcgmin
          call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%at, 1) 
       end if
    end select
    return
  end subroutine back2NDM

  subroutine setV_F(N,R,V,F,lover,lvm)

    real(double),intent(in):: R(N)
    real(double),intent(out)::V
    real(double),intent(out)::F(N)
    logical,intent(out)::lover,lvm
    integer,intent(in) ::N

    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    !-----------------------------------------------
    !-----------------------------------------------
    integer::iproc,proc_source,cellx,celly,cellz
    real(double)::aux,auy,auz
    character :: extension*2
    integer::lenfn2,ko,i1,i,i2,ip
    real(double) :: fpmax,fpn,forctot,formax,fpmax_glob
    real(double)::volu

    real(double) :: invVolu,pre,x
    logical :: lchgbox


    select case(ityprel)
    case(1)

       do i=1,atcgcomp%im
          i1=atcgcomp%num_at_glob(i)
          atcgcomp%xp(1:3,i) = R(3*i1-2:3*i1)
       end do
       lchgbox=.false.
    case(2)
       ip=0
       do i1=1,3
          do i2=1,3
             ip=ip+1
             h(i1,i2)=R(ip)
          end do
       end do
       call initbox(boxcg,h)
!       write(6,*)'R',R
       call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%at, 1) 
       lchgbox=.true.
    end select

    if (lperiod)  call periodbox (boxcg,atcgcomp)
    call depeche_mode (gcpara,'xft',lchgbox)
    V=potist
    NCALLS=NCALLS+1
    sigsym = 0.5d0*(sig + Transpose(sig) )
    select case (ityprel)
    case(1)
       do i=1,atcgcomp%im
          i1=atcgcomp%num_at_glob(i)
          F(3*i1-2:3*i1)=atcgcomp%fp(1:3,i)
          !       write(6,*)'FORCES',i,i1,F(3*i1-2:3*i1)
       end do
    case(2)

      !       do i1=1,3
       !       write(6,*)'SIG ',sigsym(:,i1)*unitP
       !       end do
       !       Pre=(sig(1,1)+sig(2,2)+sig(3,3))/3
       !       write(6,'(A,G15.7)')'pression',Pre*unitP
       call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%bg, -1) 
       h(:,:)=boxcg%at(:,:)
       trh=Transpose(h)
       call MatInv(h,invh)
       invtrh = Transpose(invh)
       volu = calcvol(h(1:3,1),h(1:3,2),h(1:3,3))
       invVolu = 1.d0/volu

       forcebox(:,:)=MatMul( sigsym(:,:) , invtrh(:,:) )*volu
       ip=0
       do i1=1,3
          do i2=1,3
             ip=ip+1
             F(ip)=forcebox(i1,i2)
          end do
       end do
    end select

    call test_conv(N,F,lover,V,R,lvm)    

    !       do i=1,N
    !       write(6,*)'R_F',i,R(i),F(i)
    !    end do
    return
  end subroutine SETV_F


  subroutine set_pointers_GC
    use parautils,only: psc_p,sig_p,potist_p,atcomp_p,cellcomp_p,box_p,div_p,atloc_p,celloc_p,ltabvois_p,&
         &itetabvois_p,it_p,lperiod_p,lchg_p
    character,target::carac(3)
    carac='xft'
    psc_p=>pscCG
    sig_p=>sig
    potist_p=>potist
    atcomp_p=>atcgcomp
    cellcomp_p=>cellcgcomp
    box_p=>boxcg
    div_p=>gcpara
    atloc_p=>atcgloc
    celloc_p=>cellcgloc
    ltabvois_p=>atcgcomp%ltabvois
    itetabvois_p=>itetabvois
    it_p=>it
    lperiod_p=>lperiod
    lchg_p=>lchg
  end subroutine set_pointers_GC

  
end module WGC_mod
