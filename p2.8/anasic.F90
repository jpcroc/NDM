module sic
    USE extension_mod,only: Pextension    


contains

  subroutine anasic(im,xp,ityp,nvi,ivois,at,it,maxvois,eatom,itmax)
    USE gen_com_m, ONLY : ev2erg,erg2ev
    implicit none
    integer::im,it,maxvois,itmax
    real*8::xp(3,im),at(3,3),eatom(im)
    integer::ivois(maxvois,im),nvi(im),ityp(im)


    integer,save ::ncall=0
    integer::lenfn
    character*9::extension
    real*8::xp1,xp2,xp3

    integer :: i,luxlf,plottyp (im),lurasmol,j,nvhomo
    logical ::lop

    character*3:: ty(2)
    integer::nviT(2),itypT(4,2)
    real*8:: eatomT(2),eatomdiff(im)

!    write(6,*)'entree anasic',it
    ncall=ncall+1
    nvhomo=0
    nviT(:)=4 ; itypT(:,1)=2; itypT(:,2)=1
    open(unit=547,file='eatomin')
    read(547,*)eatomT(1)
!    write(6,*)eatomT(1)
    read(547,*)eatomT(2)
    close (547)
!    write(6,*)'eatomin lu'
 !    eatomT(1)=-6.0102327 ;    eatomT(2)=-6.7642684
    ty(1)='Si' ; ty(2)='C '

    do i=1,im
       eatomdiff(i)=eatom(i)*erg2ev-eatomT(ityp(i))
    end do
 !   write(6,*)'post eatD'
    luxlf=839
    inquire(unit=luxlf,OPENED=lop)
    if (.not.lop) open(luxlf, file='sic.axsf', form='formatted', &
         status='unknown')

    if ((ncall==1).and.(itmax.ne.0))write(luxlf,*)'ANIMSTEPS ',itmax
    do i=1,im
 !      write(6,*)i
       if (ityp(i)==1) then
          plottyp(i)=14
       else
          plottyp(i)=6
       end if
       if (nvi(i).lt.nviT(ityp(i))) plottyp(i)=plottyp(i)-1
       if (nvi(i).gt.nviT(ityp(i))) plottyp(i)=plottyp(i)+1

       do j=1,nvi(i)
          if (ityp(ivois(j,i)).eq.ityp(i)) then 
             
             plottyp(i)=plottyp(i)+20
!             cycle
          end if
       end do
       if (plottyp(i).ge.20) nvhomo=nvhomo+1
!       write(6,*)i,nvi(i),(ityp(ivois(j,i)),j=1,nvi(i)),plottyp(i)
    end do
    if(nvhomo.gt.0)write(6,*)'Nb d''atomes avec liaison homopolaire = ', nvhomo
    write(luxlf,*)'ATOMS',ncall
    do i=1,im
       xp1 = xp(1,i)*1D+08
       xp2 = xp(2,i)*1D+08
       xp3 = xp(3,i)*1D+08
       write(luxlf, 138) plottyp(i),xp1, xp2, xp3
    end do
138 format(I2,3f10.4,I7)

    call Pextension(it,extension,lenfn)
    lurasmol=837
    open(lurasmol, file='SiC.'//extension(1:lenfn)//'.mol', form='formatted',status='unknown')




    write (lurasmol, '(I9,A,I7,A,F12.6)') im, ' IT =', it
    at=at*1.d8
    write (lurasmol,'(9F12.6)')at(1,1),at(2,1),at(3,1),at(1,2),at(2,2),at(3,2),at(1,3),at(2,3),at(3,3)
    at=at/1.d8 

    do i=1,im
       xp1 = xp(1,i)*1D+08
       xp2 = xp(2,i)*1D+08
       xp3 = xp(3,i)*1D+08
       write (lurasmol, 136) ty(ityp(i)),xp1, xp2, xp3,eatomdiff(i)   
    end do
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    close(lurasmol)

136 format(A,3f10.4,D14.5,I7)
  end subroutine anasic
end module sic
