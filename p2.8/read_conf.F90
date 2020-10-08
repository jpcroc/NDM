module read_conf
#ifdef PARA
#endif
  USE gen_com_m, ONLY:imm,rang,it,itmax,nitmax,tmean,pmean,timel,tstep,two,usdh,dilat,lvpread,&
       &oldtstep !at,bg,zls2,tstep,oldtstep,tmean,timel,nox,noy,noz,im,imm,&
!         &it,itmax,ldesinteg,lperiod,pmean,zl,xpspr,nzl,normat,rang,dilat,fmt_cin,lsuivinonpbc,&
!         &lvpread,nitmax,usdh,two
    USE arret_ndm_mod,only: arret_ndm

    use atomconfig,only:atom_config_d,atom_config,atom_config_e
    use boxconfig,only:box_config,initbox
    implicit none
    
contains
  subroutine read_cin(boxcin,itread,atcinr,immr,fnamcin,lres,fmtcin,icible,imic)
    USE T_kind_param_m, ONLY:  double
    !    USE suivinonpbc
#ifdef PARA
    !    use mpi
    !    USE mod_para,only:MPI_COMM_space,ierr,NDM_MPI_REAL_DOUBLE,status,nprocs,proc_cell
#endif
    implicit none
    character,intent(in) :: fnamcin*80
    integer,intent(in)::itread
    type(box_config)::boxcin
    class(atom_config),optional::atcinr
    integer,intent(in),optional::immr,fmtcin
    logical,intent(in),optional::lres
    integer,allocatable,optional,intent(in)::icible(:)
    integer,optional,intent(in)::imic
    
    logical::lrestart=.false.
    integer :: i, ic, icintype, icintypemod , lucin,fmt_cin=1
    integer, dimension(:),allocatable     :: ibuffer
    real(double), dimension(:,:),allocatable    :: buffer
    real(double)::at(3,3)
#ifdef PARA
    integer,      dimension(ntyp)         :: na_loc

    integer  :: pointeur_loc
    integer  :: numcell, numproc
    integer  :: i_glob
    integer  :: cellx,celly,cellz
#endif
    !      integer , dimension(imm,ntyp) :: fv
    !         integer , dimension(6000,10) :: fv    !Truc_bizarre_jmd
    integer::im_gr,i_loc
    if (present(immr))then
       allocate (ibuffer(immr))
       allocate (buffer(3,immr))
    endif
    if (present(lres))lrestart=lres
        if (present(fmtcin))fmt_cin=fmtcin
    if (rang==0) then
       write(6,*)
       write(6,*)' *-*-*-*-*-*LECTURE DE CIN*-*-*-*-*-*-'
       write(6,*)
    endif

    lucin = 93
    open(unit=lucin, file=fnamcin, form='unformatted', status='unknown', err=456)

    read (lucin, err=456) icintype

    !          if (rang==0) write (6, *) 'config type de fichier .cin : ', icintype
    if (icintype>3.or.icintype<0) then
       write (6, *) rang, 'wrong icintype'
       call arret_ndm
    endif

    icintypemod = mod(icintype,2)

    if (lrestart.and.icintypemod==0) then
       write (6, *) rang, 'not possible to restart from this file'
       call arret_ndm
    endif
    !at(vect123,xyz)
    !        if (icintype>=2) then
    read (lucin, err=456) at
    if(dilat(1).ne.0.0)then
       do i=1,3
          at(i,:)=at(i,:)*dilat(i)
       end do
    end if

    call initbox(boxcin,at)

    select case(itread)
    case(0)
       return
    case(1)
       if (.not.present(atcinr))then
          write(6,*)'atcinr pas present et itread=1'
          stop
       end if

       read (lucin, err=456) im_gr                         !number of atoms in the box
       if (im_gr>immr) then
          if(rang==0)                    write (6, *) 'im > imM', im_gr, imm
          call arret_ndm
       endif
       atcinr%im=im_gr

       !                    write(6,*)im
       read (lucin, err=456) ibuffer   !ityp
       atcinr%ityp(1:im_gr)=ibuffer(1:im_gr)

       if (rang==0) write (6, *) 'types'
       read (lucin, err=456) buffer    ! xp
       atcinr%xp(1:3,1:im_gr)=buffer(1:3,1:im_gr)
       if (rang==0) write (6, *) 'xp'
       formcin:select case (fmt_cin)
       case (0) formcin
          do i=1,im_gr
             atcinr%num_at_glob(1:im_gr) = i
          enddo
       case(1) formcin
          read (lucin, err=456) ibuffer
          atcinr%num_at_glob(1:im_gr)=ibuffer(1:im_gr)
          if (rang==0) write (6, *) 'num_at_glob'
       case default  formcin
          if (rang.eq.0) write(6,*) 'precisez le format fmt_cin'
          call arret_ndm
       end select formcin


       select type(atcinr)
       type is (atom_config)
          read (lucin, err=456) buffer                     !xpp
          read (lucin, err=456) buffer                     !vp
          lvpread=.false.
       type is (atom_config_d)
          if (icintypemod==1) then
             read (lucin, err=456) buffer                     !xpp
             atcinr%xpp(:,1:im_gr)=buffer(:,1:im_gr)
             read (lucin, err=456) buffer                     !vp
             atcinr%vp(:,1:im_gr)=buffer(:,1:im_gr)
             !             read (lucin, err=456) buffer                     !former positions
             !             atcinr%ax(:,1:im_gr)=buffer(:,1:im_gr)
             read (lucin, err=456) buffer                     !ax inutile
          end if
       type is (atom_config_e)
          if (icintypemod==1) then
             read (lucin, err=456) buffer                     !xpp
             atcinr%xpp(:,1:im_gr)=buffer(:,1:im_gr)
             read (lucin, err=456) buffer                     !vp
             atcinr%vp(:,1:im_gr)=buffer(:,1:im_gr)
             !             read (lucin, err=456) buffer                     !former positions
             !             atcinr%ax(:,1:im_gr)=buffer(:,1:im_gr)
             read (lucin, err=456) buffer                     !ax utile peut-être
             if (lrestart) then 
                if (atcinr%lax) then
                   atcinr%ax(:,1:im_gr)=buffer(:,1:im_gr)
                end if
             end if
          else
             if (atcinr%lax) then
                atcinr%ax(:,1:im_gr)=atcinr%xp(:,1:im_gr)
                lvpread=.false.
             end if
          end if
       end select
       read (lucin, err=456) oldtstep
       !  Si l'option de redemarrage (lrestart) n'est pas activee
       !  alors les positions d'origine ax deviennent les xp du fichier .cin
       !       if (.not.lrestart) then
       !          atcinr%ax(:,:im) = atcinr%xp(:,:im)
       !          if (lsuivinonpbc) axnonpbc(:,:im)=ax(:,:im)
       !       endif


       if (lrestart) then
          read (lucin, err=456) tmean, pmean, it, timel
          if (nitmax.ge.0) itmax=it+nitmax
          tstep = oldtstep

          if (rang==0) then

             write (6, *) 'restart parameters'
             write (6, *) 'it =', it, ' time =', timel
             write (6, *) 'pmean', pmean, ' tmean =', tmean
             write (6, *) 'tstep', tstep
          endif                                ! fin rang=0
          usdh = 1.0/(two*tstep)
       endif
    case(2) ! at xp et num_at_glob
       if (.not.present(atcinr))then
          write(6,*)'atcinr pas present et itread=2'
          stop
       end if

       read (lucin, err=456) im_gr                         !number of atoms in the box
       if (im_gr>immr) then
          if(rang==0)                    write (6, *) 'im > imM', im_gr, imm
          call arret_ndm
       endif
       atcinr%im=im_gr

       !                    write(6,*)im
       read (lucin, err=456) ibuffer   !ityp muet
       !       atcinr%ityp(1:im_gr)=ibuffer(1:im_gr)

       !       if (rang==0) write (6, *) 'types'
       read (lucin, err=456) buffer    ! xp
       atcinr%xp(1:3,1:im_gr)=buffer(1:3,1:im_gr)

       !       formcin:select case (fmt_cin)
       !       case (0) formcin
       !          do i=1,im_gr
       !             atcinr%num_at_glob(1:im_gr) = i
       !          enddo
       !       case(1) formcin
       read (lucin, err=456) ibuffer
       atcinr%num_at_glob(1:im_gr)=ibuffer(1:im_gr)
       !       case default  formcin
       !          if (rang.eq.0) write(6,*) 'precisez le format fmt_cin'
       !          call arret_ndm
       !       end select formcin

    case(3)
        if (.not.present(atcinr))then
          write(6,*)'atcin pas present et itread=3'
          stop
       end if

       read (lucin, err=456) im_gr                         !number of atoms in the box
       if (im_gr>immr) then
          if(rang==0)                    write (6, *) 'im > imM', im_gr, imm
          call arret_ndm
       endif
       atcinr%im=imic

       !                    write(6,*)im
       read (lucin, err=456) ibuffer   !ityp
       do i_loc=1,imic
          atcinr%ityp(i_loc)=ibuffer(icible(i_loc))
       enddo


       if (rang==0) write (6, *) 'types'
       read (lucin, err=456) buffer    ! xp muet

       select case (fmt_cin) !num_at_glob muet
       case (0) 
       !   do i=1,im_gr
       !      atcinr%num_at_glob(1:im_gr) = i
       !   enddo
       case(1) 
          read (lucin, err=456) ibuffer
       !   atcinr%num_at_glob(1:im_gr)=ibuffer(1:im_gr)
       !   if (rang==0) write (6, *) 'num_at_glob'
       end select


       select type(atcinr)
       type is (atom_config)
          read (lucin, err=456) buffer                     !xpp
          read (lucin, err=456) buffer                     !vp
          lvpread=.false.
       type is (atom_config_d)
          if (icintypemod==1) then
             read (lucin, err=456) buffer                     !xpp
             do i_loc=1,imic
                atcinr%xpp(:,i_loc)=buffer(:,icible(i_loc))
             enddo
             read (lucin, err=456) buffer                     !vp
             do i_loc=1,imic
                atcinr%vp(:,i_loc)=buffer(:,icible(i_loc))
             enddo

             !             read (lucin, err=456) buffer                     !former positions
             !             atcinr%ax(:,1:im_gr)=buffer(:,1:im_gr)
             read (lucin, err=456) buffer                     !ax inutile
          end if
       type is (atom_config_e)
          if (icintypemod==1) then
             read (lucin, err=456) buffer                     !xpp
             do i_loc=1,imic
                atcinr%xpp(:,i_loc)=buffer(:,icible(i_loc))
             enddo
             read (lucin, err=456) buffer                     !vp
             do i_loc=1,imic
                atcinr%vp(:,i_loc)=buffer(:,icible(i_loc))
             enddo
             read (lucin, err=456) buffer                     !ax utile peut-être
             if (lrestart) then 
                if (atcinr%lax) then
                   do i_loc=1,imic
                      atcinr%ax(:,i_loc)=buffer(:,icible(i_loc))
                   enddo
                end if
             end if
          else
             if (atcinr%lax) then
                atcinr%ax(:,1:imic)=atcinr%xp(:,1:imic)
                lvpread=.false.
             end if
          end if
       end select
       read (lucin, err=456) oldtstep
       !  Si l'option de redemarrage (lrestart) n'est pas activee
       !  alors les positions d'origine ax deviennent les xp du fichier .cin
       !       if (.not.lrestart) then
       !          atcinr%ax(:,:im) = atcinr%xp(:,:im)
       !          if (lsuivinonpbc) axnonpbc(:,:im)=ax(:,:im)
       !       endif


       if (lrestart) then
          read (lucin, err=456) tmean, pmean, it, timel
          if (nitmax.ge.0) itmax=it+nitmax
          tstep = oldtstep

          if (rang==0) then

             write (6, *) 'restart parameters'
             write (6, *) 'it =', it, ' time =', timel
             write (6, *) 'pmean', pmean, ' tmean =', tmean
             write (6, *) 'tstep', tstep
          endif                                ! fin rang=0
          usdh = 1.0/(two*tstep)
       endif
      
    end select
    close (lucin)
    return
456 print *,'Erreur dans la lecture du fichier .cin, verifier son format&
         & et fmt_cin ATTENTION A BIG_ENDIAN !! SI COMMPILE BIG_ENDIAN NE LIT PLUS QUE CA'


  end subroutine read_cin
  !******************************************************************************************************
   !******************************************************************************************************
  !******************************************************************************************************
  
  subroutine read_gin (boxrg,atrg,fnamgin,latr)
    USE T_kind_param_m, ONLY:  double


    character,intent(in) :: fnamgin*80
    type(atom_config)::atrg
    type(box_config)::boxrg
    integer::latr(3)

!    real(double)::rumax_init,alpha_init
    integer ,     dimension(:),   allocatable :: itypc
    real(double), dimension(:,:), allocatable :: xc
    real(double),dimension(:,:),allocatable :: tmpxc

    real(double)::at(3,3)

    integer ::  lugin, imcell, la, lb, lc, icell,ic,i,ia,ib
    
    !  si coordonnees reduites

    if (rang==0) write (6, *) '**********construction du reseau************'
    lugin=92
    !                                                !number of cells in 3 directions
     open(unit=lugin, file=fnamgin, status='unknown')
    read (lugin, *) latr(1), latr(2), latr(3)
    if (rang==0) write (6, *) 'repetition de mailles', latr


    !     **** coordonnes des vecteurs de maille en A dans une base orthonormee ****
    !                                                !a
    read (lugin, *) at(1,1), at(2,1), at(3,1)
    !b
    read (lugin, *) at(1,2), at(2,2), at(3,2)
    !                                                !c
    read (lugin, *) at(1,3), at(2,3), at(3,3)
    at=at*1d-8
    call initbox(boxrg,at)
    read (lugin, *) imcell               !number of atoms in UC
    if (imcell>imm) then
       if(rang==0)               write (6, *) 'trop d_atomes dans la cel. unite'
       call arret_ndm
    endif
    call atrg%init(imcell)
    do i = 1, imcell
       read (lugin, *) atrg%xp(1,i), atrg%xp(2,i), atrg%xp(3,i),atrg%ityp(i)
    end do
    do i=1,imcell
       WHERE ( (atrg%xp(:,i).LT.0.d0).OR.(atrg%xp(:,i).GE.1.d0) )
          atrg%xp(:,i)  = atrg%xp(:,i)  - Dble(Floor(atrg%xp(:,i)))
       END WHERE
    end do

    close(lugin)
  end subroutine read_gin

end module read_conf
