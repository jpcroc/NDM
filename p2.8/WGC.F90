module WGC_mod

  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:  inv_angst, lperiod, rang,leev,sig,cunitp, &
       iteration, itesauv, itesauvposition, itesauvforce, fnam,lenfnam,fnamcout,&
       inv_angst, erg2ev, angst,fpstop,fsumstop,itetabvois, iterasmol,&
       dmtype, potist,mdcg_noise,lspaceNDM,sigstop,sigext,ihbox0,unitP,lprahman
  USE sauvegardeT_mod,only: sauvegardeT
  USE endrunT_mod,only: endrunT
  USE arret_ndm_mod,only: arret_ndm
 
  USE initspeed_mod,only: bruit_xp
#ifdef PARA
  use Tpara,only:COMM_space,myidsp,nprocspace,para_space_config
#else
  use Tpara,only:nprocspace,para_space_config
#endif
  USE atomconfig,only : atom_config
  USE cellconfig, only:cell_config
  USE boxconfig,only:box_config
  use paraconfig,only:para_config
  USE parautils,only:initcomp,depeche_mode
  USE Mat_utils_mod,only:  MatInv
  USE boxconfig,only:box_config,periodbox,updatebox
  USE recips_mod,only: recips ,calcvol
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE rasmolT_mod,only: rasmolT
  
  implicit none
  integer,parameter::unitgc=333
  type(box_config)::boxcgmin
  type(box_config),target::boxcg
  class(atom_config),allocatable,target::atcgcomp
  class(atom_config),allocatable::atcgmin
  type(cell_config),target::cellcgcomp
  type(para_space_config),target::pscCG
  class(atom_config),pointer::atcgloc
  type(cell_config),pointer::cellcgloc
  type(cell_config),target:: cellcible ! ne sert qu'à faire pointer cellnebloc sur quelquechose
  class(atom_config),allocatable,target::atcible
  type(para_config),target::gcpara
  real(double), dimension(3,3) :: trh, invh, invtrh, forcebox,h!,sigsym
  real(double),allocatable,dimension (:)::R,F,Rmin,Fmin
  integer::Nvar,ndir,nstep,ityprel,ncgtry
  real(double)::betaguess,V,betaV,betaP,beta,betaV0,betaP0
  integer::ncalls,nextsauv,nextmol,formatsauv
  logical::lvm
  real(double)::fpstop0,fpstopsig,fstpdecr
  logical,target:: lchg,lcalcvois,lvarstop
  real(double)::  ft2,fm2,fs2,beta35
  real(double)::gammas,gammav
contains

  subroutine initsteep
    real(double), allocatable :: bruitmd(:,:)
    integer::i,i1,i2,ip
    allocate(bruitmd(3,1:atcgcomp%im))
    select case(ityprel)

    case(1)
       if (lprahman)        fpstopsig=(sigstop/unitP)*(boxcg%volu**0.6666666)
!       fpstop=fpstop0
       Nvar=atcgcomp%im*3
       if (allocated (R).or.allocated(F)) then
          deallocate(R,F,Rmin,Fmin)
       end if
       allocate(R(Nvar))
       allocate(Rmin(Nvar))
       allocate(F(Nvar))
       allocate(Fmin(Nvar))
       if (mdcg_noise /= 0 ) then
          call bruit_xp (bruitmd,atcgcomp%im)

       else
          bruitmd=0
       end if
       mdcg_noise=0
       R(:)=0;V=0;F(:)=0
       atcgmin=atcgcomp
       do i=1,atcgcomp%im
          i1=atcgcomp%num_at_glob(i)
          ! Variables = cartesian coordinates (in cm)
          atcgcomp%xp(1:3,i)= atcgcomp%xp(1:3,i)+bruitmd(1:3,i)
          R(3*i1-2:3*i1) = atcgcomp%xp(1:3,i)
       end do
       Rmin=R
       Fmin=F
    case(2)
       if (allocated (R).or.allocated(F)) then
          deallocate(R,F,Rmin,Fmin)
       end if
       fpstopsig=(sigstop/unitP)*(boxcg%volu**0.6666666)
       write(unitgc,*)'FPSTOPSIG',fpstopsig
!	fpstopsig=fpstop
       call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%bg, -1) 
       !       Nvar=9
       if (any(ihbox0.ne.1)) then
          Nvar=0
          do i1=1,3
             do i2=1,3
                if (ihbox0(i1,i2)==1) then
                   Nvar=Nvar+1
                end if
             end do
          end do
          write(unitgc,*)'Nb of cell variables',Nvar
       else
          Nvar=9
       end if

       allocate(R(Nvar))
       allocate(Rmin(Nvar))
       allocate(F(Nvar))
       allocate(Fmin(Nvar))
       R(:)=0;V=0;F(:)=0
       boxcgmin=boxcg
!!$       ip=0
!!$       do i1=1,3
!!$          do i2=1,3
!!$             ip=ip+1
!!$             R(ip)=boxcg%at(i1,i2)
!!$          end do
!!$       end do
       ip=0
       do i1=1,3
          do i2=1,3
!             hold(i1,i2)=boxcg%at(i1,i2)
             if (ihbox0(i1,i2)==1) then
                ip=ip+1
                R(ip)=boxcg%at(i1,i2)
             end if
          end do
       end do

       Rmin=R
       Fmin=F
    end select
  end subroutine initsteep

  subroutine test_conv(Nt,Ft,lover,Vt,Rt,lvm )
    integer,intent(in)::Nt
    real(double),intent(in)::Ft(Nt)
    real(double),intent(in),optional::Rt(Nt)
    real(double),optional::Vt
    logical::lover
    logical,optional::lvm
    real(double),save::Vminabs=1d16
    real(double),dimension (3,3)::invh,invtrh,forcebx
    real(double)::formax,forctot,fsigmax,deltaV,sigmax,sigm2,ppot
    integer::i,i1,i2,ip,ic
    if (present(lvm))lvm=.false.

    ft2=0;fm2=0;fs2=0
    lover=.false.
    write(unitGC,*)
    if (present(Vt)) then
       write(unitGC,*)'TESTCOMPLET'
    else
       write(unitGC,*)'TESTDIR'
    end if
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
                   if (ihbox0(i1,i2)==1) then
                      ip=ip+1
                      boxcgmin%at(i1,i2)=Rt(ip)
                   end if
                end do
             end do
             call updatebox(boxcgmin,boxcgmin%at)
          end select
       end if
    end if
    if (present(Vt)) then
       write (unitgc, *)
       write (unitgc, *) '************ STRESS in ', cunitP
       do ic = 1, 3
          write (unitgc, '(I1,3(A,I1),A,3G18.10)') ic,' sigma potentiel (1,', ic, ') (2,', ic, &
               ') (3,', ic, ') =',sig(1:3,ic)*unitP
          ppot = ppot+1.0/3.0*sig(ic,ic)
       end do
       write(unitgc,'(A,G18.10)')'PRESSURE',ppot*unitP
    end if
    select case (ityprel)
    case(2)
       forctot=sqrt( SUM(atcgcomp%fp(:,:)**2))
       formax=0
       do i=1,atcgcomp%im
          do ic=1,3
             formax=max(formax,abs(atcgcomp%fp(ic,i)))
          end do
       enddo
       forctot = forctot*erg2eV/angst
       formax  = formax*erg2eV/angst
!!$       do i=1,Nvar
!!$          formax = Max( formax,Abs(Ft(i)))
!!$       end do

       Fsigmax=0;sigmax=0
       do ip=1,Nvar
          Fsigmax=max(Fsigmax,abs(Ft(ip)))
       end do
       do i1=1,3
          do i2=1,3
             if (ihbox0(i1,i2)==1) then
                sigmax=max(sigmax,abs(unitP*sig(i1,i2)))
             end if
          end do
       end do
       sigm2=sigmax
       if (present(Vt)) then 
          write(unitgc,*)'COMP fsigmax fpstopsig',fsigmax, fpstopsig
       else
          write(unitgc,*)'DIR fsigmax fpstopsig',fsigmax, fpstopsig
       end if
       if (fsigmax.le.fpstopsig) lover=.true.

    case(1)

       call MatInv(boxcg%at(:,:),invh)
       invtrh = Transpose(invh)
       forcebx(:,:)=MatMul( sig(:,:) , invtrh(:,:) )*boxcg%volu
       fsigmax=0
       do i1=1,3
          do i2=1,3
             fsigmax=max(fsigmax,abs(forcebx(i1,i2)))
          end do
       end do

       sigm2=0
       do i1=1,3
          do i2=1,3
             sigm2=max(sigm2,abs(unitP*sig(i1,i2)))
          end do
       end do
       forctot=sqrt( SUM(Ft(:)**2))
       formax=0
       do i=1,Nvar
          formax = Max( formax,Abs(Ft(i)))
       end do

       forctot = forctot*erg2eV/angst
       formax  = formax*erg2eV/angst
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
       ft2=forctot ; fm2=formax;fs2=fsigmax
    if (present(Vt)) then
!       ft2=forctot ; fm2=formax;fs2=fsigmax
          if(lvm) then
             write(unitgc,'(I4,5E20.11,A, 1E20.11)')ncalls, Vt ,ft2,fm2,Fs2, sigm2, ' ****', deltaV
             write(6,'(I4,5E20.11,A, 1E20.11)')ncalls, Vt*erg2eV ,ft2,fm2,Fs2, sigm2, ' ****', deltaV*erg2eV
          else
             write(unitgc,'(I4,5E20.11)')ncalls, Vt ,ft2,fm2,Fs2,sigm2
             write(6,'(I4,5E20.11)')ncalls, Vt*erg2eV ,ft2,fm2,Fs2,sigm2
          end if

       if (lvm)then
          select case(ityprel)
          case(2)
             call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%at, 1) 
          end select
          if (NCALLS.ge.nextsauv) then
             formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'.cout'
             call sauvegardeT(atcgcomp,cellcgcomp,boxcg,formatsauv,fnamcout,latcomp=.true.)
             nextsauv=NCALLS+itesauv
          end if
          if (NCALLS.ge.nextmol) then
             call rasmolT(atcgcomp,boxcg,iteration,latcomp=.true.)
             nextmol=NCALLS+iterasmol
          end if
          select case(ityprel)
          case(2)
             call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%bg, -1) 
          end select
       end if
    else
       select case(ityprel)
       case(1)
          write(unitgc,'(A,2E20.11)')'test direction fm2 ft2 fs2 sigmax ',fm2,ft2,Fs2,sigmax
       case(2)
          write(unitgc,'(A,4E20.11)')'test direction fm2 ft2 fs2 sigmax ',fm2,ft2,Fs2,sigmax
       end select
    end if
    return
  end subroutine test_conv

  subroutine final_tconv(lover)
    logical,intent(out)::lover
    real(double):: forctot,sigmax,formax,fsigmax
    integer::i1,i2,i
    real(double),dimension (3,3)::invh,invtrh,forcebx
    
    forctot=sqrt( SUM(atcgcomp%fp(:,:)**2))
    formax=0
    sigmax=0
    fsigmax=0
    do i=1,atcgcomp%im
       do i2=1,3
          formax = Max( formax,Abs(atcgcomp%fp(i2,i)))
       end do
    end do
    do i=1,3
       do i2=1,3
          if (ihbox0(i,i2)==1) then
             sigmax = Max( sigmax,Abs(unitp*sig(i,i2)))
          end if
       end do
    end do
       call MatInv(boxcg%at(:,:),invh)
       invtrh = Transpose(invh)
       forcebx(:,:)=MatMul( sig(:,:) , invtrh(:,:) )*boxcg%volu
       fsigmax=0
       do i1=1,3
          do i2=1,3
             fsigmax=max(fsigmax,abs(forcebx(i1,i2)))
          end do
       end do

    forctot = forctot*erg2eV/angst
    formax  = formax*erg2eV/angst
    write(unitgc,*)
    write(unitgc,*)'  FORCE MAX            FORCETOT            FSIGMAX'
    write(unitgc,'(3E20.11)')formax,forctot,Fsigmax
    write(unitgc,'(A,3E20.11)')'thresholds',fpstop,fsumstop,fpstopsig
    write(unitgc,*)
    lover=.false.
    if (lprahman) then
       if (fpstop.gT.0) then
          if ((formax.le.fpstop).and.(fsigmax.le.fpstopsig))lover=.true.
       end if
       if (fsumstop.gT.0) then
          if ((forctot.le.fsumstop).and.(fsigmax.le.fpstopsig))lover=.true.
       end if
    else
       if (fpstop.gT.0) then
          if (formax.le.fpstop)lover=.true.
       end if
       if (fsumstop.gT.0) then
          if (forctot.le.fsumstop)lover=.true.
       end if
    end if
    ft2=forctot ; fm2=formax;fs2=fsigmax
    return
  end subroutine final_tconv



  subroutine back2NDM( N,R,V,F,lover)
    integer,intent(in)::N
    real(double),intent(in)::F(N),R(N),V
    logical,intent(in) :: lover

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
!       call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%bg, -1) 
       if (lover) then
          ip=0
          do i1=1,3
             do i2=1,3
                if (ihbox0(i1,i2)==1) then 
                   ip=ip+1
                   boxcg%at(i1,i2)=R(ip)
                end if
             end do
          end do
          call updatebox(boxcg,boxcg%at)
          call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%at, 1) 
          potist=V
       else
          boxcg=boxcgmin
          call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%at, 1) 
       end if
    end select
    return
  end subroutine back2NDM

  subroutine setV_F(N,R,V,F,lover,lvm,idesc)

    real(double),intent(in):: R(N)
    real(double),intent(out)::V
    real(double),intent(out)::F(N)
    logical,intent(out)::lover,lvm
    integer,intent(in) ::N
    integer,intent(inout)::idesc
    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    !-----------------------------------------------
    !-----------------------------------------------



    integer::i1,i,i2,ip
    real(double)::volu,Press

    real(double) :: invVolu
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
             if (ihbox0(i1,i2)==1) then
                ip=ip+1
                boxcg%at(i1,i2)=R(ip)
             end if
          end do
       end do
       call updatebox(boxcg,boxcg%at)

       call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%at, 1) 
       lchgbox=.true.
    end select
    call periodbox (boxcg,atcgcomp)
    call depeche_mode (gcpara,lchgbox)
    V=potist
    NCALLS=NCALLS+1
    select case (ityprel)
    case(1)
       do i=1,atcgcomp%im
          i1=atcgcomp%num_at_glob(i)
          F(3*i1-2:3*i1)=atcgcomp%fp(1:3,i)

       end do
    case(2)
       call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%bg, -1)
!       h(:,:)=boxcg%at(:,:)
       trh=Transpose(boxcg%at(:,:))
       call MatInv(boxcg%at(:,:),invh)
       invtrh = Transpose(invh)
       volu = calcvol(boxcg%at(1:3,1),boxcg%at(1:3,2),boxcg%at(1:3,3))
       invVolu = 1.d0/volu

       forcebox(:,:)=MatMul( sig(:,:) , invtrh(:,:) )*volu
       Press=(sig(1,1)+sig(2,2)+sig(3,3))/3
       ip=0
       do i1=1,3
          do i2=1,3
             if (ihbox0(i1,i2)==1) then
                ip=ip+1
                F(ip)=forcebox(i1,i2)
             end if
          end do
       end do
    end select
    call test_conv(N,F,lover,V,R,lvm)
    if (lvm) idesc=idesc+1
    !       do i=1,N
    !       write(unitgc,*)'R_F',i,R(i),F(i)
    !    end do
    return
  end subroutine SETV_F


  subroutine set_pointers_GC
    use parautils,only: psc_p,sig_p,potist_p,atcomp_p,cellcomp_p,box_p,div_p,atloc_p,celloc_p,lcv_p,&
         lperiod_p,lchg_p,lperiod_p,it_p
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
    lcv_p=>lcalcvois
    it_p=>iteration
    lperiod_p=>lperiod
    lchg_p=>lchg
  end subroutine set_pointers_GC

  
end module WGC_mod
