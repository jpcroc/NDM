module transf_mod
  USE gen_com_m, only:uwrt,lwrt,rang,fnam,lenfnam,fnamcout,latcomp,lwgin

  USE arret_ndm_mod,only:arret_ndm
  USE atomconfig,only:atom_config
  use cellconfig,only:cell_config
  use boxconfig,only:box_config
  use var_pot,only:ntyp
  USE sauvegardeT_mod,only: sauvegardeT!,cin2gin
  use rasmolT_mod,only:rasmolT
  use newunit_mod,only:newunit
  USE Tpara,only:COMM_space
  implicit none
  logical:: ltransf
  logical,allocatable:: lrmtype(:)
  integer,allocatable::iswitchtype2(:)

contains
  ! **************************************************************

  subroutine transfconf(atcf,celcf,boxcf)
    class(atom_config),intent(inout),target::atcf
    type(cell_config),intent(in),target::celcf
    class(box_config)::boxcf

    integer ::unitlp,i,iti,formatsauv,j
    character*80::nameo
    class(atom_config),allocatable::attrf
    namelist /inputtr/lrmtype,iswitchtype2

    !    call atcf%deftype(attrf)

    if (ltransf) then 
       allocate (attrf,mold=atcf)
       attrf=atcf
       allocate(lrmtype(ntyp))
       lrmtype(:)=.false.
       allocate(iswitchtype2(ntyp))
       do iti=1,ntyp
          iswitchtype2(iti)=iti
       end do


       call newunit(unitlp)
       open(unit=unitlp,file='transfin')
       read(unitlp,nml=inputtr)
       write(uwrt,*)lrmtype
       write(uwrt,*)iswitchtype2

       do iti=1,ntyp
          if(lrmtype(iti).eqv..true.)   then
             if (lwrt) then
                write(uwrt,*)'removal of atoms of type ',iti 
             end if
          end if
       end do
       j=0
       do i=1,atcf%im
          if (lrmtype(atcf%ityp(i)).eqv..false.) then
             j=j+1
             call atcf%copy_atom(i,attrf,j)
          end if
       end do
       attrf%im=j
       write(uwrt,*)'im',atcf%im,attrf%im
       attrf%im_glob=attrf%im
       call comm_space%sum(attrf%im_glob)
       write(uwrt,*)'im',attrf%im,attrf%im_glob


       do iti=1,ntyp
          if (iswitchtype2(iti)==iti) cycle
          if (lwrt) write(uwrt,*)'atoms of type ',iti , 'transformed to type ',iswitchtype2(iti)
          do i=1,attrf%im
             if (attrf%ityp(i)==iti) then
                attrf%ityp(i)=iswitchtype2(iti)
             end if
          end do
       end do
       formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'.trf.cout'
       call sauvegardeT(attrf,celcf,boxcf,formatsauv,fnamcout,latcomp=latcomp)
       nameo=fnam(1:lenfnam)
       nameo=nameo//'.trf'
       nameo='trf'
       call rasmolT (attrf,boxcf,namefr=nameo,latcomp=latcomp)
    else
       formatsauv = 2 ; fnamcout= fnam(1:lenfnam)//'.cout.'
       call sauvegardeT(atcf,celcf,boxcf,formatsauv,fnamcout,latcomp=latcomp)
       call rasmolT (atcf,boxcf,latcomp=latcomp)
       if (lwgin)       call rasmolT (atcf,boxcf,-1,latcomp=latcomp,ivisumol=5)
       if (rang==0) write (uwrt, *) 'generation terminee'

    end if
    !     if (lwgin) call rasmolT (attrf,boxcf,namefr=nameo,latcomp=latcomp,ivisumol=5)
    if (lwrt) write (uwrt, *) 'generation terminee'
    call arret_ndm

  end subroutine transfconf

  subroutine poscheck(xp,ityp,im,imm)
    character*1 rep
    integer im,imm,ityp(imm),iat
    real*8 xp(3,imm)

1   continue      
    write(uwrt,*)'numero de l_atome ?'
    read(5,*)iat
    write(uwrt,*)'atome no ',iat,' de type',ityp(iat),' situe en'
    write(uwrt,*)xp(1,iat)*1d+08,xp(2,iat)*1d+08,xp(3,iat)*1d+08

    write(uwrt,*)

    write(uwrt,*)'autre atome ?'
    read(5,*)rep
    if (rep.eq.'o') goto 1

    return
  end subroutine poscheck
  subroutine transf (atcf)
    !routine de transformation de la boite : ajouter, enlever, transformer des atomes,
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m
    class(atom_config)::atcf




    !-----------------------------------------------
    !   D u m m y   A r g u m e n t s
    !-----------------------------------------------
    ! Variables locales
    integer :: ichg,iv,itii,i,ic
    real(double) :: xii(3)
    character*1 rep

    write(uwrt,*)
    write(uwrt,*)'MODIFICATION DE LA BOITE '
1   continue
    write(uwrt,*)
    write(uwrt,*)'modif. de type :1; enlever un at. :2 ; ajouter un at. :3'
    read(5,*) ichg

    select case (ichg)
    case (1)
       write(uwrt,*)'verifications des positions atomiques ? o/n'
       read (5,*) rep
       if (rep.eq.'o') then
          call poscheck(atcf%xp,atcf%ityp,atcf%im,atcf%imm)
       endif

       write(uwrt,*)'numero de l_atome ?'
       read(5,*) iv
       write(uwrt,*)'atome no ',iv,' de type',atcf%ityp(iv),'situe en'
       write(uwrt,*) atcf%xp(1,iv)*1d+08,atcf%xp(2,iv)*1d+08,atcf%xp(3,iv)*1d+08

       write(uwrt,*)'nouveau type de l_atome ?'
       read(5,*)itii
       atcf%ityp(iv)=itii

    case (2)
2      continue       
       write(uwrt,*)'verifications des positions atomiques ? o/n'
       read (5,*) rep
       if (rep.eq.'o') then
          call poscheck(atcf%xp,atcf%ityp,atcf%im,atcf%imm)
       endif

       write(uwrt,*)'indiquer l indice de l_atome a supprimer '
       read(5,*)iv
       write(uwrt,*)'atome no ',iv,' de type',atcf%ityp(iv),'situe en'
       write(uwrt,*)atcf%xp(1,iv)*1d+08,atcf%xp(2,iv)*1d+08,atcf%xp(3,iv)*1d+08

       do i=iv+1,atcf%im
          atcf%ityp(i-1)=atcf%ityp(i)
          do ic=1,3
             atcf%xp(ic,i-1)=atcf%xp(ic,i)
          enddo
       enddo
       atcf%im=atcf%im-1

    case (3)
12     continue       
       write(uwrt,*)'verifications des positions atomiques ? o/n'
       read (5,*) rep
       if (rep.eq.'o') then
          call poscheck(atcf%xp,atcf%ityp,atcf%im,atcf%imm)
       endif

       write(uwrt,*)'type de l_atome ?'
       read(5,*)itii
       write(uwrt,*)'position (en A) ?'
       read(5,*)xii(1),xii(2),xii(3)
       do ic=1,3
          xii(ic)=xii(ic)*1d-08
       enddo
       atcf%im=atcf%im+1
       atcf%xp(:,atcf%im)=xii(:)
       atcf%ityp(atcf%im)=itii
       if(atcf%im.gt.atcf%imm) then
          write(uwrt,*)'im> imm' 
          call arret_ndm
       end if

    case default
       write (uwrt, *) 'mauvais type de chnagement'
       call arret_ndm
    end select

    write(uwrt,*)'AUTRE MODIFICATION ?'
    read(5,*)rep
    if (rep.eq.'o') goto 1

    return   
  end subroutine transf
end module transf_mod
