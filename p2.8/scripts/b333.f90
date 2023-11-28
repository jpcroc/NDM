program b333
  implicit none
  real*8,allocatable::depla(:),x5d(:,:),x3d(:,:),x5de(:,:),x5f(:,:),x3f(:,:),x5fe(:,:)
  integer,allocatable::ind(:),ityp(:),ityp5(:)
  integer::nat,i,ic,ibd(3),ibf(3),j,ic1,ic2,ic3,nat3e,nat3f,nat3eO,nat3fO
  real*8::a5,a3,dec

  integer::cubv(5,5,5),indmin(3),indmax(3),delta(3),indminm1(3),indmaxm1(3)
  
  delta=0; cubv=0;indmin=100;indmax=0; cubv=0
  nat3e=0;nat3f=0;nat3eO=0;nat3fO=0
  open (unit=11,file='deb_I1UO2.gin',form='formatted')
  open (unit=12,file='fin_I1UO2.gin',form='formatted')
  open (unit=21,file='deb_I1UO2.333.gin',form='formatted')
  open (unit=22,file='fin_I1UO2.333.gin',form='formatted')
  open (unit=10,file='distimages',form='formatted')

  read(11,*);read(12,*)
  read(11,*)a5;read(12,*)
  read(11,*);read(12,*)
  read(11,*);read(12,*)
  read(11,*)nat;read(12,*)
  allocate (x5d(3,nat));  allocate (x3d(3,nat));  allocate (x5de(3,nat))
  allocate (x5f(3,nat));  allocate (x3f(3,nat));  allocate (x5fe(3,nat))
  allocate (ind(nat));   allocate (ityp(nat)); allocate (ityp5(nat));   allocate (depla(nat))
  a3=3./5.*a5
  write(6,*)a5,a3
  write(21,*)'1 1 1 ';  write(22,*)'1 1 1 '
  write(21,*) a3, ' 0 0 '
  write(21,*) ' 0 ',a3, ' 0 '
  write(21,*) ' 0 0 ',a3 
  write(22,*) a3, ' 0 0 '
  write(22,*) ' 0 ',a3, ' 0 '
  write(22,*) ' 0 0 ',a3
  write(21,*)' 330 ';  write(22,*)' 330 '
  do i=1,nat
     read(11,*)x5d(1,i),x5d(2,i),x5d(3,i),ityp5(i)
     read(12,*)x5f(1,i),x5f(2,i),x5f(3,i),ityp(i)
     if (ityp(i).ne.ityp5(i)) then
        write(6,*)'BOUH!' ; stop
     end if
     x5d(:,i)=x5d(:,i)+0.01
     x5f(:,i)=x5f(:,i)+0.01
  end do
  where(x5d.ge.1) x5d=x5d-1.0
  where(x5d.lt.0) x5d=x5d+1.0
  where(x5f.ge.1) x5f=x5f-1.0
  where(x5f.lt.0) x5f=x5f+1.0
  x5fe=5*x5f
  x5de=5*x5d


  do i=1,nat
     read(10,*)ind(i),depla(i)
  end do
  call tri(nat,depla,ind)
  do  j=1,nat
     i=ind(j)
     do ic=1,3
 !       write(6,*)i,x5de(ic,i),x5fe(ic,i)
        ibd(ic)=int(x5de(ic,i))+1
        ibf(ic)=int(x5fe(ic,i))+1
        indminm1=indmin;indmaxm1=indmax
        if (indmin(ic).gt.ibd(ic))indmin(ic)=ibd(ic)
        if (indmin(ic).gt.ibf(ic))indmin(ic)=ibf(ic)
        if (indmax(ic).lt.ibd(ic))indmax(ic)=ibd(ic)
        if (indmax(ic).lt.ibf(ic))indmax(ic)=ibf(ic)
        delta(ic)=indmax(ic)-indmin(ic)
!        write(6,*)i,x5de(ic,i),x5fe(ic,i),ibd(ic),ibf(ic)
     end do
     write(6,'(5I6,G17.7,6I4)')i,j,delta(:),depla(j),indmax,indmin

     if (any(delta==3)) then
        indmin=indminm1
        indmax=indmaxm1
        exit
     end if
  end do
  write(6,*)indmin
  write(6,*)indmax
  do ic1=indmin(1),indmax(1)
     do ic2=indmin(2),indmax(2)
        do ic3=indmin(3),indmax(3)
!           write(6,*)'TT',ic1,ic,ic3
           cubv(ic1,ic2,ic3)=1
        end do
     end do
  end do
!  write(6,*)cubv
 ! stop
  do i=1,nat
     do ic=1,3
        ibd(ic)=int(x5de(ic,i))+1
        ibf(ic)=int(x5fe(ic,i))+1
!        write(6,*)i,x5de(ic,i),x5fe(ic,i),ibd(ic),ibf(ic)
     end do

     if ((cubv(ibd(1),ibd(2),ibd(3))==1).or.(cubv(ibf(1),ibf(2),ibf(3))==1)) then
        nat3e=nat3e+1 ; if (ityp(i)==1) nat3eO=nat3eO+1
        nat3f=nat3f+1; if (ityp(i)==1) nat3fO=nat3fO+1
     end if
  end do
  write(6,*)nat3e,nat3f,nat3EO,nat3fO
  if ((nat3e.ne.330).or. (nat3f.ne.330).or.(nat3eO.ne.220).or. (nat3fO.ne.220)) then
     write(6,*)'POUET'
     stop
  end if

  do i=1,nat
     do ic=1,3
        ibd(ic)=int(x5de(ic,i))+1
        ibf(ic)=int(x5fe(ic,i))+1
!        write(6,*)i,x5de(ic,i),x5fe(ic,i),ibd(ic),ibf(ic)
     end do

     if ((cubv(ibd(1),ibd(2),ibd(3))==1).or.(cubv(ibf(1),ibf(2),ibf(3))==1)) then
        do ic=1,3
           dec=(indmin(ic)-1)/5.0
           x3d(ic,i)=(x5d(ic,i)-dec)*5./3.
           x3f(ic,i)=(x5f(ic,i)-dec)*5./3.
        end do
        write(21,'(3G19.10,I4)')x3d(:,i),ityp(i)
        write(22,'(3G19.10,I4)')x3f(:,i),ityp(i)

     end if
  end do

  
end program b333
  
subroutine tri(nat,depla,ind)
  integer,intent(in)::nat
  integer::ind(nat)
  real*8::depla(nat)

  integer::i,j,itmp
  real*8::tmp

  do i=1,nat
     do j=i+1,nat
        if(depla(i).lt.depla(j))then
           tmp=depla(i)
           itmp=ind(i)
           depla(i)=depla(j)
           depla(j)=tmp
           ind(i)=ind(j)
           ind(j)=itmp
        end if
     end do
  end do
!  do i=1,nat
!     write(6,*)i,ind(i),depla(i)
!  end do
end subroutine tri
