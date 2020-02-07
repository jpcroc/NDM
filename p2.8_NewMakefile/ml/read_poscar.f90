


subroutine read_poscar(pref,l,ef)

use gen_com_m, ONLY: im,at,bg
use tab_imm_m, only : ityp,xp,fp
use var_pot, only : ntyp
use ml_in_ndm_module, only : rangml,path,lpath,build_subdata,natm

implicit none

character(len=2),intent(in) :: pref
integer,intent(in) :: l ! N° of poscar
double precision,intent(out) :: ef

integer :: i,j,fileunit,nitems
double precision :: e,ecor
character(len=60) :: fullfilenamein
character(len=13) :: filenamein
!character(len=10) :: filename
character(len=100) :: command,string
integer :: number,fnumber
character(len=7)::cnumber

lpath=len_trim(path)

fileunit=30


number=1e+6
fnumber=number+l

write(cnumber,'(i7)') fnumber

fullfilenamein=trim(adjustl(path))//cnumber(2:7)//".poscar"

!$$  if ( l.le.9 ) then
!$$      write(fullfilenamein,'(A<lpath>,A5,I1,A7)')path,"00000",l,".poscar"
!$$      write(filenamein,'(A5,I1,A7)')"00000",l,".poscar"
!$$   elseif ( l.le.99 ) then
!$$      write(fullfilenamein,'(A<lpath>,A4,I2,A7)')path,"0000",l,".poscar"
!$$      write(filenamein,'(A4,I2,A7)'),"0000",l,".poscar"
!$$  elseif ( l.le.999 ) then
!$$     write(fullfilenamein,'(A<lpath>,A3,I3,A7)')path,"000",l,".poscar"
!$$     write(filenamein,'(A3,I3,A7)')"000",l,".poscar"
!$$  elseif ( l.le.9999 ) then
!$$     write(fullfilenamein,'(A<lpath>,A2,I4,A7)')path,"00",l,".poscar"
!$$     write(filenamein,'(A2,I4,A7)')"00",l,".poscar"
!$$  elseif ( l.le.99999 ) then
!$$     write(fullfilenamein,'(A<lpath>,A1,I5,A7)')path,"0",l,".poscar"
!$$     write(filenamein,'(A1,I5,A7)')"0",l,".poscar"
!$$  else
!$$     write(fullfilenamein,'(A<lpath>,I6,A7)')path,l,".poscar"
!$$     write(filenamein,'(I6,A7)')l,".poscar"
!$$  endif
   if (rangml==0) write(*,*) 'poscar ....', fullfilenamein
   open(fileunit,file=fullfilenamein,status='unknown')
   !old read(fileunit,'(G14.9,1x,G14.9,1x,G14.9)')e,ef,ecor
   !old read(fileunit,'(3(F13.10,1x))')at(1,1),at(1,2),at(1,3)
   !old read(fileunit,'(3(F13.10,1x))')at(2,1),at(2,2),at(2,3)
   !old read(fileunit,'(3(F13.10,1x))')at(3,1),at(3,2),at(3,3)
   read(fileunit,*) e,ef,ecor
   read(fileunit,*) at(1,1),at(1,2),at(1,3)
   read(fileunit,*) at(2,1),at(2,2),at(2,3)
   read(fileunit,*) at(3,1),at(3,2),at(3,3)
!   read(fileunit,*)at(1,1),at(1,2),at(1,3)
!   read(fileunit,*)at(2,1),at(2,2),at(2,3)
!   read(fileunit,*)at(3,1),at(3,2),at(3,3)
   call recips(at(:,1), at(:,2), at(:,3), bg(:,1), bg(:,2), bg(:,3))
   read(fileunit,'(A)')string
   ntyp=nitems(string)
   if (allocated(natm)) deallocate(natm)
   allocate(natm(ntyp))
   read(fileunit,*)natm(:)
   im=sum(natm)
   read(fileunit,*)

   if (associated(xp)) deallocate(xp); allocate(xp(3,im))
   if (associated(fp)) deallocate(fp); allocate(fp(3,im))
   if (associated(ityp)) deallocate(ityp); allocate(ityp(im))

   do i=1,im
      read(fileunit,*)(xp(j,i),j=1,3)
      !read(fileunit,'(3(F9.5,1x))')(xp(j,i),j=1,3)
   enddo
!   xp(:,1:im)=matmul(at(:,:),xp(:,1:im)) ! frac to angstrom
   !read(fileunit,*)
   !do i=1,im
  !    read(fileunit,'(3(F9.6,1x))')(fp(j,i),j=1,3)
   !enddo
   close(fileunit)

   do i=1,ntyp
      do j=1,natm(i)
         if (i==1) then
            ityp(j)=i
         else
            ityp(natm(i-1)+j)=i
         endif
      enddo
   enddo

   if (build_subdata) then
      write(command,*)"cp ",fullfilenamein," ",pref,"_",filenamein
      call system(command)
   endif

return
end subroutine read_poscar

subroutine convert_A2cm(i,xa,xpos)
! i= 1 convert positions from A  -> cm
! i=-1 convert positions from cm -> A
use gen_com_m, only : A2cm,at,bg
use tab_imm_m, only : xp
use ml_in_ndm_module, ONLY: rangml

implicit none

integer,intent(in) :: i
logical,intent(in) :: xa,xpos

if (i==1) then ! A2cm
   if (xa) then
         at(:,:)=at(:,:)*A2cm
         bg(:,:)=bg(:,:)/A2cm
   end if
   if (xpos) xp(:,:)=xp(:,:)*A2cm
elseif (i==-1) then ! cm2A
   if (xa) then
         at(:,:)=at(:,:)/A2cm
         bg(:,:)=bg(:,:)*A2cm
   end if
   if (xpos) xp(:,:)=xp(:,:)/A2cm
else
   if (rangml==0) write(*,*) "Wrong input : input = 1 or -1"
   stop "fatal in convert_A2cm"
endif


return
end subroutine convert_A2cm
