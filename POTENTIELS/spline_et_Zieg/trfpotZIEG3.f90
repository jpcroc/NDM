program trfpot
  implicit none
  real*8, parameter :: ev2erg=1.6021764631580d-12, erg2eV=1.d0/eV2erg   
  integer::nptr,npt,i,j,ic,iew,ndec
  real*8,pointer::xr(:),potr(:),xp(:),potp(:),br(:),cr(:),dr(:),bp(:),cp(:),dp(:)
  real*8::stp,max,dec,catom(2),roff1,roff2
  real*8,pointer::potz(:,:),potz_d(:,:),potzt(:),potzt2(:)

  write(6,*)'step= 1 ou npt=2'
  read(5,*)ic
  write(6,*)ic
  write(6,*)'Ziegler=1 ou pas=2'
  read(5,*)iew
  write(6,*)iew
  select case(ic)
  case(1)
     write(6,*)'stp,max'
     read(5,*)stp,max
     write(6,*)stp,max
     npt=int(max/stp)
     write(6,*)'==> npt',npt
  case(2)
     write(6,*)'npt max'
     read(5,*)npt,max
     write(6,*)npt, max
     stp=max/npt
     write(6,*)'==>stp',stp
  end select
  allocate(xp(npt))
  allocate(potp(npt))
  allocate(bp(npt))
  allocate(cp(npt))
  allocate(dp(npt))

  write(6,*)'nptr'
  read(5,*)nptr
  write(6,*)nptr
  allocate(xr(nptr))
  allocate(potr(nptr))
  allocate(br(nptr))
  allocate(cr(nptr))
  allocate(dr(nptr))
  open(unit=10,file='potrin')
  open(unit=11,file='potrout')
  if (iew.eq.1) then 
     open(unit=12,file='Zpotrout')
     open(unit=14,file='Zpotdecrawrout')
     open(unit=15,file='Zpotrawrout')
  end if
  do i=1,nptr
     read(10,*)xr(i),potr(i)
  end do

  call cspline (nptr,xr,potr,br,cr,dr)

  do i=1,npt
     xp(i)=float(i)*stp
     if(xp(i).lt.xr(1))then 
        potp(i)=potr(1)
     else
        do j=1,nptr-1
           if(xr(j+1).gt.xp(i))exit
        end do
!        write(6,'(A,2I5,2G15.7)')'i j xp xr ', i,j, xp(i),xr(j)
        dec=-(xr(j)-xp(i))
!        write(6,*)dec,br(j)
        potp(i)=potr(j)+dec*(br(j)+dec*(cr(j)+dec*dr(j)))
     end if
     write(11,*)xp(i),potp(i)
!     write(6,*)xp(i),potp(i)
  end do

  if (iew==2) stop

  call cspline (npt,xp,potp,bp,cp,dp)

  allocate (potz(4,0:npt+1))
  allocate (potzt(0:npt+1))
  allocate (potzt2(0:npt+1))
  allocate (potz_d(4,0:npt+1))
  potz(1,1:npt)=potp(1:npt)*ev2erg
  potz(2,1:npt)=bp(1:npt)*ev2erg/1d-8
  potz(3,1:npt)=cp(1:npt)*ev2erg/(1d-8*1d-8)
  potz(4,1:npt)=dp(1:npt)*ev2erg/(1d-8*1d-8*1d-8)
  write(6,*)'catom 1 2'
  read(5,*)catom(1),catom(2)
  write(6,*)catom(1),catom(2)
  write(6,*)'roff1 2'
  read(5,*)roff1,roff2
  write(6,*)roff1,roff2
  read(5,*)ndec
  write(6,*)ndec

  do i=1,npt
!     xp(i)=float(i)*stp
!     if(xp(i).lt.xr(1))then 
!        potp(i)=potr(1)
!     else
!        do j=1,nptr-1
!           if(xr(j+1).gt.xp(i))exit
!        end do
!        dec=xr(j)-xp(i)
!        potp(i)=potr(j)+dec*(br(j)+dec*(cr(j)+dec*dr(j)))
!     end if
     write(13,*)xp(i)*1d-8,potz(1,i)
!     write(6,*)xp(i),potp(i)
  end do




  roff1=roff1*1d-8
  roff2=roff2*1d-8
  stp=stp*1d-8
  call zieg2(potz, potz_d, stp,npt, catom, roff1, roff2,potzt,potzt2,ndec)
  potz=potz*erg2ev
  potzt=potzt*erg2ev
  potzt2=potzt2*erg2ev

  do i=1,npt
     write(12,'(2F25.8)')xp(i),potz(1,i)
     write(14,'(2F25.8)')xp(i),potzt(i)
     write(15,'(2F25.8)')xp(i),potzt2(i)

  end do



end program trfpot
