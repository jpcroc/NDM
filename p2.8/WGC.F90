module WGC_mod

  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:  inv_angst, lperiod, rang,itmax,leev,sig, &
       it, itesauv, itesauvposition, itesauvforce,itmax, fnam,lenfnam,fnamcout,&
       inv_angst, erg2ev, angst,fpstop,fsumstop,itetabvois, &
       dmtype, potist,mdcg_noise,formatsauv,lspaceNDM,latcomp,sigstop,sigext,ihbox0,unitP
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
  use paraconfig,only:para_config,initparapuresp
  USE parautils,only:initcomp,pointer_caltabt_calfo
  USE Mat_utils_mod,only:  MatInv
  USE scalebox_mod,only: scalebox
  USE boxconfig,only:box_config,periodbox,initbox
 USE recips_mod,only: recips ,calcvol
USE cryst_to_cart_mod,only: cryst_to_cart

  implicit none

  type(box_config)::boxcg,boxcgmin
  type(atom_config)::atcgcomp,atcgmin
  type(cell_config)::cellcgcomp
  type(para_space_config)::pscCG
  class(atom_config),pointer::atcgloc
  type(cell_config),pointer::cellcgloc
  type(cell_config),target:: cellcible ! ne sert qu'à faire pointer cellnebloc sur quelquechose
  type(atom_config),target::atcible
  type(para_config)::gcpara
  real(double), dimension(3,3) :: trh, invh, invtrh, forcebox,h,sigsym
    real(double),allocatable,dimension (:)::R,F
  integer::N,ndir,nstep,ityprel
  real(double)::betaguess,V
  integer::ncalls
  logical::lvm
  real(double)::fpstop0
  
contains

  subroutine initsteep
    real(double), allocatable :: bruitmd(:,:)
    integer::i,i1,i2,ip
    allocate(bruitmd(3,1:atcgcomp%im))
    select case(ityprel)

    case(1)
       N=atcgcomp%im*3
       allocate(R(N))
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

    case(2)
       unitP=1d-9
       fpstop0=fpstop
       fpstop=(sigstop/unitP)*(boxcg%volu**0.6666666)
       write(6,*)'SIGSTOP==FPSTOP=',fpstop
       call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%bg, -1) 
       N=9
       allocate(R(N))
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

    end select
  end subroutine initsteep

  subroutine test_conv(Nt,Ft,lover,Vt,Rt )
    integer,intent(in)::Nt
    real(double),intent(in)::Ft(Nt)
    real(double),intent(in),optional::Rt(Nt)
    real(double),optional::Vt
    logical::lover
    real(double)::Vminabs=1d16

    real(double)::forctot,formax,deltaV,sigmax,fsigmax
    integer::i,i1,i2,ip
    lvm=.false.
    !    write(6,*)'FT',ft
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
    end if

  end subroutine test_conv

  subroutine final_tconv
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

    write(6,*)'END CG FORCE MAX       FORCETOT      SIGMAX'
    write(6,'(3E20.11)')formax,forctot,sigmax
    return
  end subroutine final_tconv
    

  
  subroutine back2NDM( N,R,V,F,lover)
    integer,intent(in)::N
    real(double),intent(in)::F(N),R(N),V
    logical :: lover

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
    
  subroutine setV_F(N,R,V,F,lover)

    real(double),intent(in):: R(N)
    real(double),intent(out)::V
    real(double),intent(out)::F(N)
    logical,intent(out)::lover
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
    logical:: lchg


    select case(ityprel)
    case(1)

       do i=1,atcgcomp%im
          i1=atcgcomp%num_at_glob(i)
          atcgcomp%xp(1:3,i) = R(3*i1-2:3*i1)
       end do
    case(2)
       ip=0
       do i1=1,3
          do i2=1,3
             ip=ip+1
             h(i1,i2)=R(ip)
          end do
       end do
       call initbox(boxcg,h)
       call cryst_to_cart (atcgcomp%im, atcgcomp%xp, boxcg%at, 1) 

    end select

    if (lperiod)  call periodbox (boxcg,atcgcomp)
    lchg=.true.
    call pointer_caltabt_calfo(sig,potist,atcgcomp,cellcgcomp,boxcg,atcgloc,cellcgloc,gcpara,lperiod,&
         &atcgcomp%ltabvois,it,itetabvois,lchg,pscCG,'xft') 
    V=potist
    NCALLS=NCALLS+1
    select case (ityprel)
    case(1)
       do i=1,atcgcomp%im
          i1=atcgcomp%num_at_glob(i)
          F(3*i1-2:3*i1)=atcgcomp%fp(1:3,i)
          !       write(6,*)'FORCES',i,i1,F(3*i1-2:3*i1)
       end do
    case(2)
       sigsym = 0.5d0*(sig + Transpose(sig) )
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

    call test_conv(N,F,lover,V,R)    

       !       do i=1,N
    !       write(6,*)'R_F',i,R(i),F(i)
    !    end do 
    return
  end subroutine SETV_F
  
end module WGC_mod
