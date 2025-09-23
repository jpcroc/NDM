program depla
  implicit none
  integer::nat,nat2,i,j,ic
  real*8, allocatable::xp(:,:),xref(:,:)
  real*8::at(3,3),bg(3,3),atref(3,3),bgref(3,3),dec(1,3),rcut,distab(100),dist
  character*2,allocatable::ty(:)
  character*2::tyr
  character (len=256):: fileref,filecur,line
  character*21::h1,h2,h3,h4
  write(6,*)'filecur'
  read(5,*) filecur
  open (unit=11,file=filecur,form='formatted')
  read(11,*)h1,h2,h3,h4,nat
  write(6,*)h1,h2,h3,h4,nat
  write(6,*)'fileref'
  read(5,*) fileref
  open (unit=10,file=fileref,form='formatted')
  read(10,*)h1,h2,h3,h4,nat2
  write(6,*)h1,h2,h3,h4,nat2
  if (nat.ne.nat2) then
     write(6,*)'stop nat nat2'
     stop
  end if
  allocate(xp(3,nat))
  allocate(xref(3,nat))
  allocate(ty(nat))
  read(10,*) !1
  read(11,*)
  read(10,*)!UC1
  read(11,*)

  write(6,*)'CUR'
  read(11,'(A)')line
  read(line(11:),*)at(1,1)
  write(6,*)line(1:10),at(1,1)
  read(11,'(A)')line
  read(line(11:),*)at(2,1)
  write(6,*)line(1:10),at(2,1)
  read(11,'(A)')line
  read(line(11:),*)at(3,1)
  write(6,*)line(1:10),at(3,1)
  read(11,*)

    read(11,'(A)')line
  read(line(11:),*)at(1,2)
  write(6,*)line(1:10),at(1,2)
  read(11,'(A)')line
  read(line(11:),*)at(2,2)
  write(6,*)line(1:10),at(2,2)
  read(11,'(A)')line
  read(line(11:),*)at(3,2)
  write(6,*)line(1:10),at(3,2)
  read(11,*)

    read(11,'(A)')line
  read(line(11:),*)at(1,3)
  write(6,*)line(1:10),at(1,3)
  read(11,'(A)')line
  read(line(11:),*)at(2,3)
  write(6,*)line(1:10),at(2,3)
  read(11,'(A)')line
  read(line(11:),*)at(3,3)
  write(6,*)line(1:10),at(3,3)


  write(6,*)
  write(6,*)'REF'
  read(10,'(A)')line
  read(line(10:),*)atref(1,1)
  write(6,*)line(1:10),atref(1,1)
  read(10,'(A)')line
  read(line(10:),*)atref(2,1)
  write(6,*)line(1:10),atref(2,1)
  read(10,'(A)')line
  read(line(10:),*)atref(3,1)
  write(6,*)line(1:10),atref(3,1)
  read(10,*)

    read(10,'(A)')line
  read(line(10:),*)atref(1,2)
  write(6,*)line(1:10),atref(1,2)
  read(10,'(A)')line
  read(line(10:),*)atref(2,2)
  write(6,*)line(1:10),atref(2,2)
  read(10,'(A)')line
  read(line(10:),*)atref(3,2)
  write(6,*)line(1:10),atref(3,2)
  read(10,*)

    read(10,'(A)')line
  read(line(10:),*)atref(1,3)
  write(6,*)line(1:10),atref(1,3)
  read(10,'(A)')line
  read(line(10:),*)atref(2,3)
  write(6,*)line(1:10),atref(2,3)
  read(10,'(A)')line
  read(line(10:),*)atref(3,3)
  write(6,*)line(1:10),atref(3,3)


  
  call recips (atref(1:3,1), atref(1:3,2), atref(1:3,3), bgref(1:3,1), bgref(1:3,2), bgref(1:3,3))
  call recips (at(1:3,1), at(1:3,2), at(1:3,3), bg(1:3,1), bg(1:3,2), bg(1:3,3))
     
  read(10,*);  read(10,*)
  read(11,*);  read(11,*)

  do i=1,nat
     read(11,*)
     read(11,'(A2,A5,I15)')tyr,h1,j
 !    write(6,*)tyr,h1,j
     read(11,*)xp(1,j),xp(2,j),xp(3,j)
     ty(j)=tyr
!     write(6,*)j,ty(j),xp(:,j)
  end do
  write(6,*)'REF'
  do i=1,nat
     read(10,*)
     read(10,'(A2,A5,I15)')tyr,h1,j
!     write(6,*)tyr,h1,j
     read(10,*)xref(1,j),xref(2,j),xref(3,j)
     if (ty(j).ne.tyr) then
        write(6,*) 'stop types'
        stop
     end if
!     write(6,*)j,ty(j),xp(:,j)
  end do
  distab=0
  do i=1,nat
     dec(1,:)=xp(:,i)-xref(:,i)
     do ic=1,3
        if ( (dec(1,ic).GT.0.5d0).OR.(dec(1,ic).LT.-0.5d0) )then
           dec(1,ic) = dec(1,ic) - Dble(Nint(dec(1,ic)))
        end if
     end do
     call cryst_to_cart (1, dec, at, 1) !cryst vers cart sur cv
     dist=sqrt(dec(1,1)**2+dec(1,2)**2+dec(1,3)**2)
     j=1+int(dist*10)
     distab(j)=distab(j)+1
     
  end do

  do j=100,1,-1
     if(distab(j).ne.0) then
        write(6,*)'distab',float(j)/10,distab(j)
     end if
  end do
end program depla

  
subroutine cryst_to_cart(nvec, vec, trmat, iflag)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------

  !             Version du 01 septembre 2000
  !-----------------------------------------------------------------------
  !
  !     This routine transforms the atomic positions or the k-point
  !     components from crystallographic to carthesian coordinates ( iflag=1)
  !     and viceversa ( iflag=-1 ).
  !     Output carth. coordinates are stored in the input ('vec') array.
  !
  !
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer , intent(in) :: nvec
  integer , intent(in) :: iflag
  real*8 , intent(inout) :: vec(3,nvec)
  real*8 , intent(in) :: trmat(3,3)
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: nv
  real*8, dimension(3) :: vau
  !-----------------------------------------------
  do nv = 1, nvec
     if (iflag==1) then
        vau = trmat(:,1)*vec(1,nv)+trmat(:,2)*vec(2,nv)+trmat(:,3)*vec(3,nv&
             )
     else
        vau = trmat(1,:)*vec(1,nv)+trmat(2,:)*vec(2,nv)+trmat(3,:)*vec(3,nv&
             )
     endif
     vec(:,nv) = vau

  end do
  !
  return
end subroutine cryst_to_cart


subroutine recips(a1, a2, a3, b1, b2, b3)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------

  !---------------------------------------------------------------------
  !
  !   This routine generates the reciprocal lattice vectors b1,b2,b3
  !   given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.
  !
  !
  !     first the input variables
  !
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real*8 , intent(in) :: a1(3)
  real*8 , intent(in) :: a2(3)
  real*8 , intent(in) :: a3(3)
  real*8 , intent(out) :: b1(3)
  real*8 , intent(out) :: b2(3)
  real*8 , intent(out) :: b3(3)
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: iperm, i, j, k, l, ipol
  real*8 :: den, s
  !-----------------------------------------------
  !
  !   then the local variables
  !
  !
  !    first we compute the denominator
  !
  den = 0
  i = 1
  j = 2
  k = 3
  s = 1.D0
  do iperm = 1, 3
     den = den+s*a1(i)*a2(j)*a3(k)
     l = i
     i = j
     j = k
     k = l
  end do
  i = 2
  j = 1
  k = 3
  s = -s
  do while(s<0.D0)
     do iperm = 1, 3
        den = den+s*a1(i)*a2(j)*a3(k)
        l = i
        i = j
        j = k
        k = l
     end do
     i = 2
     j = 1
     k = 3
     s = -s
  end do
  !
  !    here we compute the reciprocal vectors
  !
  i = 1
  j = 2
  k = 3
  do ipol = 1, 3
     b1(ipol) = (a2(j)*a3(k)-a2(k)*a3(j))/den
     b2(ipol) = (a3(j)*a1(k)-a3(k)*a1(j))/den
     b3(ipol) = (a1(j)*a2(k)-a1(k)*a2(j))/den
     l = i
     i = j
     j = k
     k = l
  end do
  return
end subroutine recips
