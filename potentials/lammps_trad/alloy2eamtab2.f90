program alloy2eamtab
  implicit none
  integer::i,j,k,iti,ipr,ntyp,npair,nrho,nr,ip,ip2,ios,idum,nlr,nlrho
  real*8, allocatable:: rho(:,:,:),emb(:,:,:),rep(:,:,:)
  real*8::rcut,dr,drho
  integer,allocatable:: ipolmp(:,:),ipondm(:,:),tndm(:,:),tlmp(:,:)
  open(10,file='eam.alloy',form='formatted')
  open(11,file='eamtab.end',form='formatted')
  read(10,*);  read(10,*);  read(10,*)
  read(10,*) ntyp
  npair=ntyp*(ntyp+1)/2
  read(10,*) nrho,drho,nr,dr,rcut

  allocate (emb(ntyp,2,nrho))
  allocate (rho(ntyp,2,nr))
  allocate (rep(npair,2,nr))
  allocate (ipondm(ntyp,ntyp))
  allocate (ipolmp(ntyp,ntyp))
  allocate (tndm(2,npair))
  allocate (tlmp(2,npair))
  k=0
  do i=1,ntyp
     do j=i,ntyp
        k=k+1
        ipondm(i,j)=k
        ipondm(j,i)=k
        tndm(1,k)=i
        tndm(2,k)=j
     end do
  end do
  k=0
  do i=1,ntyp
     do j=1,i
        k=k+1
        ipolmp(i,j)=k
        ipolmp(j,i)=k
        tlmp(1,k)=i
        tlmp(2,k)=j
     end do
  end do

  do i=1,nr
     rho(:,1,i)=i*dr
     rep(:,1,i)=i*dr
  end do
  do i=1,nrho
     emb(:,1,i)=i*drho
  end do
  nlr=nr/5
  nlrho=nrho/5
  do iti=1,ntyp
     read(10,*) idum
     write(6,*)idum
     do j=1,nlrho
        k=1+5*(j-1)
        read(10,*,iostat=ios)emb(iti,2,k),emb(iti,2,k+1),emb(iti,2,k+2),emb(iti,2,k+3),emb(iti,2,k+4)
!        write(6,*)emb(iti,2,k),emb(iti,2,k+4)
     end do
     do j=1,nlr
        k=1+5*(j-1)
        read(10,*,iostat=ios)rho(iti,2,k),rho(iti,2,k+1),rho(iti,2,k+2),rho(iti,2,k+3),rho(iti,2,k+4)
!        write(6,*)rho(iti,2,k),rho(iti,2,k+4)
     end do
  end do
  do ip=1,npair
     write(6,*)ip
     do j=1,nlr
        k=1+5*(j-1)
        read(10,*,iostat=ios)rep(iti,2,k),rep(iti,2,k+1),rep(iti,2,k+2),rep(iti,2,k+3),rep(iti,2,k+4)
!        write(6,*)rep(iti,2,k),rep(iti,2,k+4)
     end do
  end do

  do iti=1,ntyp
     write(11,*)iti
     write(11,*) nrho,drho
     do j=1,nrho
        write(11,*)emb(iti,1,j),emb(iti,2,j)
     end do
     write(11,*)iti
     write(11,*) nr,dr
     do j=1,nr
        write(11,*)rho(iti,1,j),rho(iti,2,j)
     end do
  end do
  ip2=0
  do i=1,ntyp
     do j=i,ntyp
        ip=ipolmp(i,j)
        ip2=ip2+1
        write(11,*)ip2, ip
        write(11,*)nr,dr
        do k=1,nr
           write(11,*)rep(ip,1,k),rep(ip,2,k)
        end do
     end do
  end do
end program alloy2eamtab
           
  
