module transf_mod
  USE gen_com_m, only:uwrt,lwrt
  USE arret_ndm_mod,only:arret_ndm
  USE atomconfig,only:atom_config
  implicit none 
contains
  ! **************************************************************
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
end module transf_mod
