module layer_mod
  USE cryst_to_cart_mod
  USE gen_com_m, ONLY:at,bg,im,imd,imfree,imgi,imgs,imm,it,rang,rulayer,nzl,frozen,free
  implicit none
contains

  subroutine layer
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE tab_imm_m

#ifdef PARA
    USE mod_para
#endif 


    implicit none
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer :: j1, j2, j3, i, i1, i2, i3, j
    integer , dimension(imm) :: itypbis, itypter, itypcar
    integer :: na1, na2, na3, na4,nfr
    real(double), pointer  :: xpdyn(:,:),xppdyn(:,:),vpdyn(:,:)
    real(double), pointer :: lay(:),typtemp(:)
    real(double) :: csup,cinf
    !-----------------------------------------------
    !
    !
    !-------------------------------------------------------------------

    csup = 1-rulayer/nzl(1)
    cinf = rulayer/nzl(1)
    if ((rang==0).and.(it.le.20)) then  
       write (6, *) 'sub layer',rulayer*1d8
       !  csup = 1-rulayer/nzl(1)
       !  cinf = rulayer/nzl(1)
       write (*, *) 'surface superieur', csup
       write (*, *) 'surface inferieur', cinf
    end if


#ifdef PARA
    free(:)=.true.
    frozen(:,:)=.false.
    call cryst_to_cart (imm, xp, bg, -1) !cart vers cryst
    where (xp(1,:im)>csup) 
       free(:im) = .false.
       frozen(1,:im)=.true.
       frozen(2,:im)=.true.
       frozen(3,:im)=.true.
    elsewhere (xp(1,:im)<cinf)
       free(:im) = .false.
       frozen(1,:im)=.true.
       frozen(2,:im)=.true.
       frozen(3,:im)=.true.
    end where
    call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart
    imd=im
    nfr=0
    do i=1,im
       if (.not.free(i)) then
          vp(:,i)=0.d0
          nfr=nfr+1
          !        write(6,'(A,I3,A,I3,A,I7)')'it ',it,' Rg ',rang,' ng ',num_at_glob(i)
       endif
    end do
    if (it.le.20) write(6,*)'rang nfr',rang,nfr
    call MPI_ALLREDUCE(nfr,imFree,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    if ((rang==0).and.(it.le.20)) write(6,*)'nb d atomes libres IMFREE ',imFree
#else
    allocate(xpdyn(3,imm))
    allocate(xppdyn(3,imm))
    allocate(vpdyn(3,imm))
    allocate(typtemp(imm))
    allocate(lay(imm))
    lay(:im) = 0

    call cryst_to_cart (imm, xp, bg, -1) !cart vers cryst
    where (xp(1,:im)>csup) lay(:im) = 1
    where (xp(1,:im)<cinf) lay(:im) = 2
    call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart

    i1 = 0
    i2 = 0
    i3 = 0
    do i = 1, im
       if (lay(i)==0) then
          i1 = i1+1
          xpdyn(1,i1) = xp(1,i)
          xpdyn(2,i1) = xp(2,i)
          xpdyn(3,i1) = xp(3,i)
          typtemp(i1)=ityp(i)
          xppdyn(:,i1)=xpp(:,i) ;  vpdyn(:,i1)=vp(:,i) ; ax(:,i1)=0.
          free(i1)=.true.
          frozen(:,i1)=.false.
       endif
    end do
    imd=i1
    write(6,*)'imd',imd
    do i = 1, im
       if (lay(i)==1) then
          i2 = i2+1
          i1=i1+1
          free(i1)=.false.
          frozen(:,i1)=.true.
          xpdyn(1,i1) = xp(1,i)
          xpdyn(2,i1) = xp(2,i)
          xpdyn(3,i1) = xp(3,i)
          typtemp(i1)=ityp(i)
          xppdyn(:,i1)=xp(:,i) ;  vpdyn(:,i1)=0. ; ax(:,i1)=xpp(:,i)
       endif
    end do
    do i = 1, im
       if (lay(i)==2) then
          i3 = i3+1
          i1=i1+1
          xpdyn(1,i1) = xp(1,i)
          xpdyn(2,i1) = xp(2,i)
          xpdyn(3,i1) = xp(3,i)
          typtemp(i1)=ityp(i)
          xppdyn(:,i1)=xp(:,i) ;  vpdyn(:,i1)=0. ; ax(:,i1)=xpp(:,i)
          free(i1)=.false.
          frozen(:,i1)=.true.
       endif
    end do


    imgs = i2
    imgi = i3
    imfree=imd
    !      if (imd+imgs+imgi/=im) then
    write (6, *) ' layer', imd, imgs, imgi, im
    !         stop
    !      endif

    ityp(:)=typtemp(:)
    xp(:,:)=xpdyn(:,:)
    xpp(:,:)=xppdyn(:,:)
    vp(:,:)=vpdyn(:,:)


    !  xpp(:,imd+1:im)=xp(:,imd+1:im)
    !  free(imd:im) = .false.
    !  do i=1,im
    !     write(6,*)i,lay(i),free(i)
    !  end do
    deallocate(typtemp) ; deallocate(xpdyn)
    deallocate(xppdyn); deallocate(vpdyn)
    deallocate (lay)
    !write(6,*)'RRRRRRRRRRRRRR',xp(1,38603)

#endif
    return
  end subroutine layer
end module layer_mod
