module WGC_mod

  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:  inv_angst, lperiod, rang,itmax,leev,sig, &
       it, itesauv, itesauvposition, itesauvforce,itmax, fnam,lenfnam,fnamcout,&
       inv_angst, erg2ev, angst,fpstop,fsumstop,itetabvois, &
       dmtype, potist,mdcg_noise,formatsauv,lspaceNDM,latcomp,sigstop,sigext,ihbox0
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
  USE boxconfig,only:box_config,periodbox

  implicit none

  type(box_config)::boxcg
  type(atom_config)::atcgcomp,atcgmin
  type(cell_config)::cellcgcomp
  type(para_space_config)::pscCG
  class(atom_config),pointer::atcgloc
  type(cell_config),pointer::cellcgloc
  type(cell_config),target:: cellcible ! ne sert qu'à faire pointer cellnebloc sur quelquechose
  type(atom_config),target::atcible
  type(para_config)::gcpara
  real(double), dimension(3,3)  ::trh0,invh0,invtrh0
  real(double)  ::volu0, invVolu0
  real(double),allocatable,dimension (:)::R,F
  integer::N,ndir,nstep,typrel
  real(double)::betaguess,Vminabs=1d16,V
  integer::ncalls
  logical::lvm
  
contains

  subroutine initsteep
    real(double), allocatable :: bruitmd(:,:)
    integer::i,i1
    allocate(bruitmd(3,1:atcgcomp%im))
    if (typrel==1) then
      
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

    else ! LPR
    end if

  end subroutine initsteep

  subroutine test_conv(Nt,Ft,lover,Vt,Rt )
    integer,intent(in)::Nt
    real(double),intent(in)::Ft(Nt)
    real(double),intent(in),optional::Rt(Nt)
    real(double),optional::Vt
    logical::lover


    real(double)::forctot,formax,deltaV
    integer::i,i1
    lvm=.false.
    !    write(6,*)'FT',ft
    forctot=sqrt( SUM(Ft(:)**2))
    formax=0
    do i=1,N
       formax = Max( formax,Abs(Ft(i)))
    end do
    lover=.false.

    if (typrel==1) then
       forctot = forctot*erg2eV/angst
       formax  = formax*erg2eV/angst
       if (present(Vt)) then

          if (Vt.lt.Vminabs) then
             deltaV=Vminabs-Vt
             do i=1,atcgmin%im
                i1=atcgmin%num_at_glob(i)
                ! Variables = cartesian coordinates (in cm)
                atcgmin%xp(1:3,i)=Rt(3*i1-2:3*i1) 
                Vminabs=Vt
             end do
             write(6,'(I4,4E20.11,A, 2E20.11)')ncalls, Vt,Vt*erg2eV ,forctot,formax,' ****', deltaV, deltaV*erg2eV
          else
             write(6,'(I4,4E20.11)')ncalls, Vt,Vt*erg2eV ,forctot,formax

          end if
       end if
       if (fpstop>0) then   
          if (formax.le.fpstop) then
!!$             if (rang==0) write(6,*)'force par atome  max  ev/Ang ', formax
!!$             if (rang==0) write (6, *) 'energie ', potist*erg2eV
!!$             !                 if (it.le.1) xp(:,:)=ax(:,:)
             lover=.true.
          end if
       end if
       if (fsumstop>0) then   
          if (forctot.le.fsumstop) then
!!$             if (rang==0) write(6,*)'  sqrt ( sum_f F_i^2 ):   ev/Ang ', forctot
!!$             if (rang==0) write(6, *) 'energie ', potist*erg2eV
             !                 if (it.le.1) xp(:,:)=ax(:,:)
             lover=.true.
             !call endrunT(atcgcomp,cellcgcomp,boxcg,latcomp)
          end if
       end if
    end if

  end subroutine test_conv


  subroutine back2NDM( N,R,V,F,lover)
    integer,intent(in)::N
    real(double),intent(in)::F(N),R(N),V
    logical :: lover

    integer::i,i1

    
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
    integer::lenfn2,ko,i1,i
    real(double) :: fpmax,fpn,forctot,formax,fpmax_glob


    logical:: lchg

    do i=1,atcgcomp%im
       i1=atcgcomp%num_at_glob(i)
       atcgcomp%xp(1:3,i) = R(3*i1-2:3*i1)
    end do
    if (lperiod)  call periodbox (boxcg,atcgcomp)
    lchg=.true.
    call pointer_caltabt_calfo(sig,potist,atcgcomp,cellcgcomp,boxcg,atcgloc,cellcgloc,gcpara,lperiod,&
         &atcgcomp%ltabvois,it,itetabvois,lchg,pscCG,'xft') 
    V=potist
    NCALLS=NCALLS+1
    do i=1,atcgcomp%im
       i1=atcgcomp%num_at_glob(i)
       F(3*i1-2:3*i1)=atcgcomp%fp(1:3,i)
!       write(6,*)'FORCES',i,i1,F(3*i1-2:3*i1)
    end do
    call test_conv(N,F,lover,V,R)    

!       do i=1,N
!       write(6,*)'R_F',i,R(i),F(i)
!    end do 
    return
  end subroutine SETV_F

end module WGC_mod
