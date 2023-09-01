module calcdepla_mod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE var_pot, ONLY:ntyp,ty
  USE gen_com_m, ONLY:tdepla,iteration,timel,iko,lcasca,lfilm,rang&
       &,ivisu,lspacendm,lcasca
  USE cellconfig,only:cell_config
  USE atomconfig,only:atom_config_e,atom_config
  use boxconfig,only: box_config
  use rasmolT_mod,only:rasmolT
  USE Tpara,only:COMM_space,nprocspace,myidsp,nprocs
  USE arret_ndm_mod,only: arret_ndm
  use vect_dist_mod,only:distat



  implicit none 
contains
  ! *******************************************************************
  subroutine calcdepla(atcf,celcf,boxcf,tdep,iteration,C1)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    implicit none

    class(atom_config_e)::atcf
    type(cell_config)::celcf
    type(box_config)::boxcf
    real(double)::tdep
    character*1,optional::C1
    integer,intent(in)::iteration


    type(atom_config)::atdep
    integer :: ndeplatot,i,iatdep,im,iti,lufilm,est_present
    integer , dimension(ntyp) :: ndepla
    integer,allocatable:: deplat(:)
    real(double)::ddep,deptot
    character*80::namefilm
    integer::lufilmpaf
    real(double)::xp_iko(4),dispiko
    integer,save::icall=0
    integer::source
!    integer::ityp_iko
    
!!$  integer :: i, iti
!!$  integer , dimension(im) :: indic
!!$  real(double), dimension(im) :: distdepl
!!$  integer :: lufilm, lufilmpaf, lufilmext, lutampon
!!$  real(double), dimension(ntyp) :: dr2
!!$  real(double) :: dri2, a1, a2, a3, c1, c2, c3, racdri2,depiko
!!$  real(double), dimension(1,3) :: cv
!!$  character :: fnamtampon*10, fnamfilmext*20, extension*10
!!$  real(double), dimension(3) :: xp_iko
!!$  integer :: ityp_iko
!!$#ifdef PARA
!!$  integer,      allocatable :: ityp_depla(:)
!!$  integer,      allocatable :: indic_depla(:)
!!$  real(double), allocatable :: xp_depla(:,:)
!!$  real(double), allocatable :: dist_depla(:)
!!$  real(double), dimension(ntyp) :: dr2_glob
!!$  integer , dimension(ntyp) :: ndepla_glob
!!$  integer :: ndeplatot_glob
!!$  integer :: est_present
!!$  integer :: ndeplatot_tmp
!!$  integer :: proc_source
!!$#endif

    icall=icall+1
    lufilmpaf = 79                             ! index fichier film du paf pour toutes les iterations
    if((myidsp==0).and.(icall==1))    open(unit=lufilmpaf, file='filmpaf', status='unknown')

    !      write(6,*)'entree dans calcdepla'
    if (rang==0) then
       write (6, *)
       write (6, *) '----------- Displacements TDEPLA ----------------',TDEPLA*1d8
    end if
    !       write(6,*)'tdepla',tdepla
    ndeplatot = 0
    ndepla(:ntyp) = 0
    if (tdep.le.0)then
       allocate(deplat(atcf%im))
    end if

    im=atcf%im
    atcf%lgul(:)=.false.
    do i=1,im
       call distat(atcf%xp(:,i),atcf%ax(:,i),boxcf,ddep)
       if (lcasca) then
          if (atcf%num_at_glob(i)==iko) dispiko=ddep
       end if
       if (ddep.ge.tdep) then
          ndeplatot=ndeplatot+1
          ndepla(atcf%ityp(i))=ndepla(atcf%ityp(i))+1
          atcf%lgul(i)=.true.
       end if
    end do
    if (tdepla.gt.0) allocate(deplat(ndeplatot))
    iatdep=0
    do i=1,im
       if (atcf%lgul(i)) then
          call distat(atcf%xp(:,i),atcf%ax(:,i),boxcf,ddep)
          iatdep=iatdep+1
          deplat(iatdep)=ddep
       end if
    end do
    if(iatdep.ne.ndeplatot) then
       write(6,*)'PB calcdepla'
       call arret_ndm
    end if

    call atcf%fab(atdep,lback=.false.)
    if (present(C1))then
       namefilm='FILM'//C1
    else
       namefilm='film'
    end if
    if (lfilm) then
       if (lspacendm.and.(nprocspace.gt.1)) then
          call rasmolT(atdep,boxcf,iteration,namefilm,latcomp=.false.,ivisumol=ivisu)
       else
          call rasmolT(atdep,boxcf,iteration,namefilm,latcomp=.true.,ivisumol=ivisu)
       end if
    end if
    
    deptot=0
    do iatdep=1,ndeplatot
       deptot=deptot+deplat(iatdep)
    end do

#ifdef PARA
    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call comm_space%sum(ndeplatot)
       call comm_space%sum(ndepla)
       call comm_space%sum(deptot)
    end if
#endif

    if (myidsp==0) then
       write(6,*)'TOTAL DISPLACED ATOMS',ndeplatot,iteration
       do iti=1,ntyp
          if (ndepla(iti).gt.0) then
             write(6,*)'TYPE DISPLACED ATOMS',iti, ndepla(iti)
          end if
       end do
       write(6,*)'TOTAL DISPLACEMENT ',deptot*1d8,iteration
    end if

    est_present=0
     if (lcasca.and.lfilm) then
        est_present=0
        do i=1, im
           if (atcf%num_at_glob(i)==iko) then
              est_present=1
              xp_iko(1:3)=atcf%xp(:,i)
              xp_iko(4)=dispiko
              exit
           endif
        enddo
        
#ifdef PARA
        
        ! Recherche du proc possedant iko
        if((lspacendm).and.(nprocspace.gt.1)) then
        ! Emission/reception des infos vers le proc 0
           if (est_present==1.and.myidsp==0) then
              !rien xp_iko deja trouvé
              !              xp_iko(:) = atcf%xp(:,i2iko)
           else if (est_present==1) then
              call comm_space%send(xp_iko,0,11001)
           else if (myidsp==0) then
              call comm_space%probe(11001,source)
              call comm_space%recv(xp_iko,11001,source)
           end if
        end if
#endif
        if(myidsp==0)then
           write (lufilmpaf, *) ' IT', iteration, ' time ', timel
           write (lufilmpaf, *)  xp_iko(1)*1D+8, xp_iko(2)*1D+8, &
                xp_iko(3)*1D+8, xp_iko(4)*1d8
        end if
     endif


    
  end subroutine calcdepla
end module calcdepla_mod
