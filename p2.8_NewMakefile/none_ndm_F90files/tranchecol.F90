program tranchecol
  !calcule les positions successives le long d'un col pour input drag ou NEB
  implicit none
  integer :: i,j,ic,im,ntr,itr,rep(3)
  integer,pointer::ityp(:)

  real*8::dr,at(3,3),deltr


  real*8,pointer::xi(:,:),xf(:,:),dx(:,:),xtr(:,:)
  character :: fnam*15
  character(len=2) :: extension
  !  write(6,*)'nb d atomes'

 
  open(unit=14, file='init_final', status='old')
  open(unit=15, file='contrainte', status='unknown')
  read(14,*)rep(:)
  read(14,*)at(:,1)
  read(14,*)at(:,2)
  read(14,*)at(:,3)
  read(14,*)im
  allocate (xi(3,im));  allocate (dx(3,im));  allocate (xf(3,im));  allocate (xtr(3,im))
  allocate(ityp(im))

  do i=1,im
     read(14,*)xi(:,i),ityp(i)
  end do
  read(14,*)
  read(14,*)
  read(14,*)
  read(14,*)
  read(14,*)
  do i=1,im
     read(14,*)xf(:,i)
  end do
  dx(:,:)=xf(:,:)-xi(:,:)
  dr=sqrt(SUM(dx*dx))

  do i=1,im
     write(15,'(3f14.8,i5,f14.8)') dx(1:3,i),  i, dr
  end do

  write(6,*)'ntr'
  read(5,*)ntr
  ntr=ntr+1
  do itr=1,ntr
     deltr=float(itr)/float(ntr)
     write(extension,'(i2.2)') itr
     fnam='tranche'//extension
     open(unit=16, file=fnam)

     xtr(:,:)=xi(:,:)+dx(:,:)*deltr
     write(16,*)rep(:)
     write(16,*)at(:,1)
     write(16,*)at(:,2)
     write(16,*)at(:,3)
     write(16,*)im
     do i=1,im
        write(16,'(3F14.8,I3)')xtr(:,i),ityp(i)
     end do
     close(16)
  end do
end program tranchecol
