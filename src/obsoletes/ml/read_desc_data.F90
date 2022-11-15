module read_desc_data_mod
        implicit none
        contains

subroutine read_descriptor(path,lpath,l,descriptor_type)

use ml_in_ndm_module, only : rangml

implicit none

!interface ! read_g2ml
!   subroutine read_g2ml(path,lpath,l,g2_out,nb,n_eta,n_rs)
!   
!   implicit none
!   
!   integer,intent(in) :: lpath
!   character(len=lpath),intent(in) :: path
!   integer,intent(in) :: l
!   double precision,dimension(:,:,:),pointer :: g2_out
!   integer,intent(out) :: nb,n_eta,n_rs
!   
!   end subroutine read_g2ml
!!end interface read_g2ml
!
!!interface read_g3ml
!   subroutine read_g3ml(path,lpath,l,g3_out,nb,n_eta,n_lambda,n_zeta)
!   
!   implicit none
!   
!   integer,intent(in) :: lpath
!   character(len=lpath),intent(in) :: path
!   integer,intent(in) :: l
!   double precision,dimension(:,:,:,:),pointer :: g3_out
!   integer,intent(out) :: nb,n_eta,n_lambda,n_zeta
!   
!   end subroutine read_g3ml
!end interface  ! read_g3ml


integer,intent(in) :: lpath
character(len=lpath),intent(in) :: path
integer,intent(in) :: l
integer,intent(in) :: descriptor_type
!real(double), allocatable, dimension(:,:,:) :: g2
!real(double), allocatable, dimension(:,:,:,:) :: g3
double precision,dimension(:,:,:),pointer :: g2
double precision,dimension(:,:,:,:),pointer :: g3
integer :: im,n_eta,n_rs,n_lambda,n_zeta
integer,parameter ::  descriptor_g2=1,descriptor_g3=2,descriptor_behler=3


nullify(g2,g3)

select case(descriptor_type)
   case(descriptor_g2)
      call read_g2ml(path,lpath,l,g2,im,n_eta,n_rs)
   case(descriptor_g3)
      call read_g3ml(path,lpath,l,g3,im,n_eta,n_lambda,n_zeta)
   case(descriptor_behler)
      call read_g2ml(path,lpath,l,g2,im,n_eta,n_rs)
      call read_g3ml(path,lpath,l,g3,im,n_eta,n_lambda,n_zeta)
   case default
     if (rangml==0) write(*,*) "No descriptor type has been selected",descriptor_type
     if (rangml==0) write(*,*) "File : read_desc_data.F90"
end select
   

!deallocate(g2,g3)

return
end subroutine read_descriptor


subroutine read_g2ml(path,lpath,l,g2_out,nb,n_eta,n_rs)

implicit none

integer,intent(in) :: lpath
character(len=lpath),intent(in) :: path
integer,intent(in) :: l
!double precision,allocatable,dimension(:,:,:),intent(out) :: g2
double precision,dimension(:,:,:),pointer :: g2_out
integer,intent(out) :: nb,n_eta,n_rs

double precision,dimension(:),allocatable :: eta,rs
character(len=90) :: command
integer :: i,p1,p2,poscarunit,g2unit,na,ni
double precision :: rcut
character(len=60) :: poscar,g2ml
character(len=8) :: t1,t2
character(len=7) :: cnumber
integer :: fnumber, number

number=1e+6
fnumber=number+l
write(cnumber,'(i7)') fnumber
poscar=trim(adjustl(path))//cnumber(2:7)//".poscar"
g2ml=trim(adjustl(path))//cnumber(2:7)//".g2ml"
poscarunit=30
g2unit=31
!$$   if ( l.le.9 ) then
!$$      write(poscar,'(A<lpath>,A5,I1,A7)')path,"00000",l,".poscar"
!$$      write(g2ml,'(A<lpath>,A5,I1,A5)')path,"00000",l,".g2ml"
!$$   elseif ( l.le.99 ) then
!$$      write(poscar,'(A<lpath>,A4,I2,A7)')path,"0000",l,".poscar"
!$$      write(g2ml,'(A<lpath>,A4,I2,A5)')path,"0000",l,".g2ml"
!$$   elseif ( l.le.999 ) then
!$$      write(poscar,'(A<lpath>,A3,I3,A7)')path,"000",l,".poscar"
!$$      write(g2ml,'(A<lpath>,A3,I3,A5)')path,"000",l,".g2ml"
!$$   elseif ( l.le.9999 ) then
!$$      write(poscar,'(A<lpath>,A2,I4,A7)')path,"00",l,".poscar"
!$$      write(g2ml,'(A<lpath>,A2,I4,A5)')path,"00",l,".g2ml"
!$$   elseif ( l.le.99999 ) then
!$$      write(poscar,'(A<lpath>,A1,I5,A7)')path,"0",l,".poscar"
!$$      write(g2ml,'(A<lpath>,A1,I5,A5)')path,"0",l,".g2ml"
!$$   else
!$$      write(poscar,'(A<lpath>,I6,A7)')path,l,".poscar"
!$$      write(g2ml,'(A<lpath>,I6,A5)')path,l,".g2ml"
!$$   endif



   open(poscarunit,file=poscar,status='unknown')
   read(poscarunit,*)
   read(poscarunit,*)
   read(poscarunit,*)
   read(poscarunit,*)
   read(poscarunit,*)
   read(poscarunit,'(2(I3,1x))')na,ni
   nb=na+ni
   close(poscarunit)

   write(command,*)"grep -c eta ",g2ml," > neta.txt" 
   call system(command)
   write(command,*)"grep -c rs ",g2ml," >> neta.txt" 
   call system(command)
   open(40,file='neta.txt',action="read")
   read(40,*)n_eta
   read(40,*)n_rs
   if (n_rs==0) then
       n_rs=1
   endif
   close(40)
   call system('rm neta.txt')
   nullify(g2_out)
   allocate(g2_out(n_eta,n_rs,nb),eta(n_eta),rs(n_rs))

   open(g2unit,file=g2ml,status='unknown')

   do p1=1,n_rs
      do p2=1,n_eta
         read(g2unit,'(A3,1X,F5.3,1X,A4,1X,F6.3)')t1,eta(p1),t2,rcut
         read(g2unit,*)
         do i=1,na
            read(g2unit,*)g2_out(p2,p1,i)
         enddo
         read(g2unit,*)
         do i=1,ni
            read(g2unit,*)g2_out(p2,p1,na+i)
         enddo
      enddo
   enddo

   close(g2unit)
   deallocate(eta,rs) 

return
end subroutine read_g2ml


subroutine read_g3ml(path,lpath,l,g3_out,nb,n_eta,n_lambda,n_zeta)

implicit none

integer,intent(in) :: lpath
character(len=lpath),intent(in) :: path
integer,intent(in) :: l
!double precision,allocatable,dimension(:,:,:,:),intent(out) :: g3
double precision,dimension(:,:,:,:),pointer :: g3_out
integer,intent(out) :: nb,n_eta,n_lambda,n_zeta

double precision,dimension(:),allocatable :: eta,zeta,lambda
character(len=150) :: command
integer :: i,p1,p2,p3,poscarunit,g3unit,na,ni
double precision :: rcut
character(len=60) :: poscar,g3ml
character(len=8) :: t1,t2,t3,t4
integer :: fnumber, number
character(len=7) :: cnumber

number=1e+6
fnumber=number+l
write(cnumber,'(i7)') fnumber
poscar=trim(adjustl(path))//cnumber(2:7)//".poscar"
g3ml=trim(adjustl(path))//cnumber(2:7)//".g3ml"

poscarunit=30
g3unit=31

!$$  if ( l.le.9 ) then
!$$      write(poscar,'(A<lpath>,A5,I1,A7)')path,"00000",l,".poscar"
!$$      write(g3ml,'(A<lpath>,A5,I1,A5)')path,"00000",l,".g3ml"
!$$   elseif ( l.le.99 ) then
!$$      write(poscar,'(A<lpath>,A4,I2,A7)')path,"0000",l,".poscar"
!$$      write(g3ml,'(A<lpath>,A4,I2,A5)')path,"0000",l,".g3ml"
!$$   elseif ( l.le.999 ) then
!$$      write(poscar,'(A<lpath>,A3,I3,A7)')path,"000",l,".poscar"
!$$      write(g3ml,'(A<lpath>,A3,I3,A5)')path,"000",l,".g3ml"
!$$   elseif ( l.le.9999 ) then
!$$      write(poscar,'(A<lpath>,A2,I4,A7)')path,"00",l,".poscar"
!$$      write(g3ml,'(A<lpath>,A2,I4,A5)')path,"00",l,".g3ml"
!$$   elseif ( l.le.99999 ) then
!$$      write(poscar,'(A<lpath>,A1,I5,A7)')path,"0",l,".poscar"
!$$      write(g3ml,'(A<lpath>,A1,I5,A5)')path,"0",l,".g3ml"
!$$   else
!$$      write(poscar,'(A<lpath>,I6,A7)')path,l,".poscar"
!$$      write(g3ml,'(A<lpath>,I6,A5)')path,l,".g3ml"
!$$   endif
   open(poscarunit,file=poscar,status='unknown')
   read(poscarunit,*)
   read(poscarunit,*)
   read(poscarunit,*)
   read(poscarunit,*)
   read(poscarunit,*)
   read(poscarunit,'(2(I3,1x))')na,ni
   nb=na+ni
   close(poscarunit)

   write(command,*)"grep eta ",g3ml," | awk '{print $2}' | uniq | wc -l > neta.txt" 
   call system(command)
   write(command,*)"grep zeta ",g3ml," | awk '{print $6}' | sort -n | uniq | wc -l >> neta.txt" 
   call system(command)
   write(command,*)"grep lambda ",g3ml," | awk '{print $8}' | sort -n | uniq | wc -l >> neta.txt" 
   call system(command)
   open(40,file='neta.txt',action="read")
   read(40,*)n_eta
   read(40,*)n_zeta
   read(40,*)n_lambda
   close(40)
   call system('rm neta.txt')

   nullify(g3_out)
   allocate(g3_out(n_eta,n_lambda,n_zeta,nb),eta(n_eta),zeta(n_zeta),lambda(n_lambda))

   open(g3unit,file=g3ml,status='unknown')

   do p3=1,n_zeta
      do p2=1,n_lambda
         do p1=1,n_eta
            read(g3unit,'(A3,1X,F6.3,1X,A4,1X,F6.3,1X,A4,1X,F6.3,1X,A6,1X,F6.3)')t1,eta(p1),t2,rcut,t3,zeta(p3),t4,lambda(p2)
            read(g3unit,*)
            do i=1,na
               read(g3unit,*)g3_out(p1,p2,p3,i)
            enddo
            read(g3unit,*)
            do i=1,ni
               read(g3unit,*)g3_out(p1,p2,p3,na+i)
            enddo
         enddo
      enddo
   enddo

   close(g3unit)
   deallocate(eta,zeta,lambda) 

return
end subroutine read_g3ml

end module
