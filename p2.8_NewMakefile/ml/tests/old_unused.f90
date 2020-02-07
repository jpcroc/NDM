

subroutine write_bispectrum_so4(bispectrum_so4,bispectrum_so4f)

  USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im,imm
  use ml_in_ndm_module, ONLY: n2_typ,j_max,path,l => i_poscar

  implicit none

  double complex,dimension(0:int(2*j_max),0:int(2*j_max),0:int(2*j_max),n2_typ,imm),intent(in) :: bispectrum_so4
  double complex,dimension(0:int(2*j_max),0:int(2*j_max),0:int(2*j_max),3,n2_typ,imm),intent(in) :: bispectrum_so4f
  character*70 :: bispectrum_so4ml, bispectrum_so4fml
  integer :: j,p1,p2,p3,ctyp,lpath,bispectrum_so4unit,bispectrum_so4funit

  namelist /input_ml/ path
  lpath=len_trim(path)

  call generate_name_for_write_bispectrum_so4(l, bispectrum_so4ml, bispectrum_so4fml)


  bispectrum_so4unit=43
  bispectrum_so4funit=44

   open(bispectrum_so4unit,file=bispectrum_so4ml,status='unknown')
   open(bispectrum_so4funit,file=bispectrum_so4fml,status='unknown')

   do ctyp=1,n2_typ
      write(bispectrum_so4unit,'(A3,1X,I2)')"typ",ctyp
      write(bispectrum_so4funit,'(A3,1X,I2)')"typ",ctyp
      do p3=0,int(2*j_max)
!        do p2=0,l_max
         p2=p3
         do p1=abs(p3-p2),min(int(2*j_max),p3+p2)
            write(bispectrum_so4unit,'(A8,2X,I2,2X,I2,2X,I2)') "l  l2 l1",p1,p2,p3
            write(bispectrum_so4funit,'(A8,2X,I2,2X,I2,2X,I2)')"l  l2 l1",p1,p2,p3
            write(bispectrum_so4funit,'(A12,1X,2(A15,1X))')"bispectrum_so4f_x","bispectrum_so4f_y","bispectrum_so4f_z"
            do j=1,im

               write(bispectrum_so4unit,'(G16.8)')bispectrum_so4(p1,p2,p3,ctyp,j)
               write(bispectrum_so4funit,'(3(G16.8,1x))')bispectrum_so4f(p1,p2,p3,1,ctyp,j),bispectrum_so4f(p1,p2,p3,2,ctyp,j),bispectrum_so4f(p1,p2,p3,3,ctyp,j)

            enddo
         enddo
!        enddo
      enddo
   enddo

   close(bispectrum_so4unit)
   close(bispectrum_so4funit)
return
end subroutine write_bispectrum_so4

subroutine generate_name_for_write_bispectrum_so4(l, bispectrum_so4ml, bispectrum_so4fml)
! input l
! output pow_so4ml, pow_so4fml


  USE T_kind_param_m, ONLY:  double
  use ml_in_ndm_module, ONLY: path
  integer, intent(in) :: l
  character*70, intent(out) :: bispectrum_so4ml, bispectrum_so4fml
  character*80 :: chfmt, chfmtf
  integer :: lpath


  lpath=len_trim(path)

   if ( l.le.9 ) then
      write(chfmt,*)'(A',lpath,',A5,I1,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I1,A10)'

      write(bispectrum_so4ml,FMT=chfmt)path,"00000",l,".bisso4ml"
      write(bispectrum_so4fml,FMT=chfmtf)path,"00000",l,".bisso4fml"
   elseif ( l.le.99 ) then
      write(chfmt,*)'(A',lpath,',A5,I2,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I2,A10)'
      write(bispectrum_so4ml,FMT=chfmt)path,"0000",l,".biso4ml"
      write(bispectrum_so4fml,FMT=chfmtf)path,"0000",l,".biso4fml"
   elseif ( l.le.999 ) then
      write(chfmt,*)'(A',lpath,',A5,I3,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I3,A10)'
      write(bispectrum_so4ml,FMT=chfmt)path,"000",l,".biso4ml"
      write(bispectrum_so4fml,FMT=chfmtf)path,"000",l,".biso4fml"
   elseif ( l.le.9999 ) then
      write(chfmt,*)'(A',lpath,',A5,I4,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I4,A10)'
      write(bispectrum_so4ml,FMT=chfmt)path,"00",l,".biso4ml"
      write(bispectrum_so4fml,FMT=chfmtf)path,"00",l,".biso4fml"
   elseif ( l.le.99999 ) then
      write(chfmt,*)'(A',lpath,',A5,I5,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I5,A10)'
      write(bispectrum_so4ml,FMT=chfmt)path,"0",l,".biso4ml"
      write(bispectrum_so4fml,FMT=chfmtf)path,"0",l,".biso4fml"
   else
      write(chfmt,*)'(A',lpath,',A5,I6,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I6,A10)'
      write(bispectrum_so4ml,FMT=chfmt)path,l,".biso4ml"
      write(bispectrum_so4fml,FMT=chfmtf)path,l,".biso4fml"
   endif


return

end subroutine generate_name_for_write_bispectrum_so4




!<begin from pow_so4

subroutine write_pow_so4(pow_so4,pow_so4f)

  USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im,imm
  use ml_in_ndm_module, ONLY: n2_typ,j_max,path,jj_max, l => i_poscar

  implicit none

  double precision,dimension(0:int(2*j_max),imm),intent(in) :: pow_so4
  double precision,dimension(0:int(2*j_max),3,imm),intent(in) :: pow_so4f
  character*70 :: pow_so4ml, pow_so4fml
  integer :: pow_so4unit,pow_so4funit
  character (len=60) :: CHFMT, CHFMTf
   write(*,*) jj_max, imm,im, int(2*j_max)

   call generate_name_for_write_powso4(l, pow_so4ml, pow_so4fml)

   pow_so4unit=41
   pow_so4funit=42

   write(*,*) 'heeree2'

   open(pow_so4unit,file=pow_so4ml,status='unknown')
   open(pow_so4funit,file=pow_so4fml,status='unknown')
   write(CHFMT,*)'(i6, 1x,',int(2*j_max)+1,'e20.10)'
   write(CHFMTf,*)'(i6, 1x,',int(2*j_max)+1,'e20.10)'
         !do j=1,im
            !write(pow_so4unit,FMT=CHFMT) j, pow_so4(0:jj_max,j)
            !write(pow_so4funit,FMT=CHFMTf) j, pow_so4f(0:int(2*j_max),1, j)
            !write(pow_so4funit,FMT=CHFMTf) j, pow_so4f(0:int(2*j_max),2, j)
            !write(pow_so4funit,FMT=CHFMTf) j, pow_so4f(0:int(2*j_max),3, j)
         !enddo

   close(pow_so4unit)
   close(pow_so4funit)
return
end subroutine write_pow_so4


subroutine generate_name_for_write_powso4(l, pow_so4ml, pow_so4fml)
! input l
! output pow_so4ml, pow_so4fml


  USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im,imm
  use ml_in_ndm_module, ONLY: path
  integer, intent(in) :: l
  character*70, intent(out) :: pow_so4ml, pow_so4fml
 character*80 :: chfmt, chfmtf
 integer :: lpath


  lpath=len_trim(path)


   if ( l.le.9 ) then
      write(chfmt,*)'(A',lpath,',A5,I1,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I1,A10)'
      write(pow_so4ml,FMT=chfmt)path,"00000",l,".powso4ml"
      write(pow_so4fml,FMT=chfmtf)path,"00000",l,".powso4fml"
   elseif ( l.le.99 ) then
      write(chfmt,*)'(A',lpath,',A5,I2,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I2,A10)'
      write(pow_so4ml,fmt=chfmt)path,"0000",l,".powso4ml"
      write(pow_so4fml,FMT=chfmtf)path,"0000",l,".powso4fml"
   elseif ( l.le.999 ) then

      write(chfmt,*)'(A',lpath,',A5,I3,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I3,A10)'
      write(pow_so4ml,fmt=chfmt)path,"000",l,".powso4ml"
      write(pow_so4fml,FMT=chfmtf)path,"000",l,".powso4fml"
   elseif ( l.le.9999 ) then

      write(chfmt,*)'(A',lpath,',A5,I4,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I4,A10)'
      write(pow_so4ml,fmt=chfmt)path,"00",l,".powso4ml"
      write(pow_so4fml,FMT=chfmtf)path,"00",l,".powso4fml"
   elseif ( l.le.99999 ) then

      write(chfmt,*)'(A',lpath,',A5,I5,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I5,A10)'
      write(pow_so4ml,FMT=chfmt)path,"0",l,".powso4ml"
      write(pow_so4fml,FMT=chfmtf)path,"0",l,".powso4fml"
   else
      write(chfmt,*)'(A',lpath,',A5,I6,A9)'
      write(chfmtf,*)'(A',lpath,',A5,I6,A10)'
      write(pow_so4ml,FMT=chfmt)path,l,".powso4ml"
      write(pow_so4fml,FMT=chfmtf)path,l,".powso4fml"
   endif


return

end subroutine generate_name_for_write_powso4
!<end   from pow_so4



subroutine write_g3(g3,g3f)

  USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im,imm
  use ml_in_ndm_module, ONLY: n3_typ,eta => eta_desc,zeta => zeta_desc,lambda => lambda_desc,n_eta,n_zeta,n_lambda,path,l => i_poscar

  implicit none

  double precision,intent(in),dimension(n_eta,n_lambda,n_zeta,n3_typ,imm) :: g3
  double precision,intent(in),dimension(n_eta,n_lambda,n_zeta,3,n3_typ,imm) :: g3f
  character*70 :: g3ml, g3fml
  integer :: j,p1,p2,p3,ctyp,lpath,g3unit,g3funit
character(len=7) :: cnumber
integer :: fnumber, number
namelist /input_ml/ path
lpath=len_trim(path)

number=1e+6
fnumber=number+l
write(cnumber,'(i7)') fnumber
g3ml=trim(adjustl(path))//cnumber(2:7)//".g3ml"
g3fml=trim(adjustl(path))//cnumber(2:7)//".g3fml"


g3unit=33
g3funit=34

!$$ if ( l.le.9 ) then
!$$     write(g3ml,'(A<lpath>,A5,I1,A5)')path,"00000",l,".g3ml"
!$$     write(g3fml,'(A<lpath>,A5,I1,A6)')path,"00000",l,".g3fml"
!$$  elseif ( l.le.99 ) then
!$$     write(g3ml,'(A<lpath>,A4,I2,A5)')path,"0000",l,".g3ml"
!$$     write(g3fml,'(A<lpath>,A4,I2,A6)')path,"0000",l,".g3fml"
!$$  elseif ( l.le.999 ) then
!$$     write(g3ml,'(A<lpath>,A3,I3,A5)')path,"000",l,".g3ml"
!$$     write(g3fml,'(A<lpath>,A3,I3,A6)')path,"000",l,".g3fml"
!$$  elseif ( l.le.9999 ) then
!$$     write(g3ml,'(A<lpath>,A2,I4,A5)')path,"00",l,".g3ml"
!$$     write(g3fml,'(A<lpath>,A2,I4,A6)')path,"00",l,".g3fml"
!$$  elseif ( l.le.99999 ) then
!$$     write(g3ml,'(A<lpath>,A1,I5,A5)')path,"0",l,".g3ml"
!$$     write(g3fml,'(A<lpath>,A1,I5,A6)')path,"0",l,".g3fml"
!$$  else
!$$     write(g3ml,'(A<lpath>,I6,A5)')path,l,".g3ml"
!$$     write(g3fml,'(A<lpath>,I6,A6)')path,l,".g3fml"
!$$  endif
!$$
   open(g3unit,file=g3ml,status='unknown')
   open(g3funit,file=g3fml,status='unknown')

   do ctyp=1,n3_typ
      write(g3unit,'(A3,1X,I2)')"typ",ctyp
      write(g3funit,'(A3,1X,I2)')"typ",ctyp
      do p3=1,n_zeta
         do p2=1,n_lambda
            do p1=1,n_eta
               write(g3ml,'(A3,1X,F6.3,1X,A4,1X,F6.3)')"eta",eta(p1),"lambda",lambda(p2),"zeta",zeta(p3)
               write(g3fml,'(A3,1X,F6.3,1X,A4,1X,F6.3)')"eta",eta(p1),"lambda",lambda(p2),"zeta",zeta(p3)
               write(g3fml,'(A12,1X,2(A15,1X))')"g3f_x","g3f_y","g3f_z"
               do j=1,im

                  write(g3unit,'(G16.8)')g3(p1,p2,p3,ctyp,j)
                  write(g3funit,'(3(G16.8,1x))')g3f(p1,p2,p3,1,ctyp,j),g3f(p1,p2,p3,2,ctyp,j),g3f(p1,p2,p3,3,ctyp,j)

               enddo
            enddo
         enddo
      enddo
   enddo

   close(g3unit)
   close(g3funit)
return
end subroutine write_g3


subroutine write_g2(g2,g2f)

USE T_kind_param_m, ONLY:  double
use gen_com_m, ONLY: im,imm
use ml_in_ndm_module, ONLY: n2_typ,eta => eta_desc,rs => rs_desc,n_eta,n_rs,path,l => i_poscar

implicit none

double precision,dimension(n_eta,n_rs,n2_typ,imm),intent(in) :: g2
double precision,dimension(n_eta,n_rs,3,n2_typ,imm),intent(in) :: g2f
character*70 :: g2ml, g2fml
integer :: j,p1,p2,ctyp,g2unit,g2funit,lpath
integer :: number, fnumber
character(len=7) :: cnumber
namelist /input_ml/ path
lpath=len_trim(path)

number=1e+6
fnumber=number+l
write(cnumber,'(i7)') fnumber
g2ml=trim(adjustl(path))//cnumber(2:7)//".g2ml"
g2fml=trim(adjustl(path))//cnumber(2:7)//".g2fml"

g2unit=31
g2funit=32

!$$   if ( l.le.9 ) then
!$$      write(g2ml,'(A<lpath>,A5,I1,A5)')path,"00000",l,".g2ml"
!$$      write(g2fml,'(A<lpath>,A5,I1,A6)')path,"00000",l,".g2fml"
!$$   elseif ( l.le.99 ) then
!$$     write(g2ml,'(A<lpath>,A4,I2,A5)')path,"0000",l,".g2ml"
!$$     write(g2fml,'(A<lpath>,A4,I2,A6)')path,"0000",l,".g2fml"
!$$  elseif ( l.le.999 ) then
!$$     write(g2ml,'(A<lpath>,A3,I3,A5)')path,"000",l,".g2ml"
!$$     write(g2fml,'(A<lpath>,A3,I3,A6)')path,"000",l,".g2fml"
!$$  elseif ( l.le.9999 ) then
!$$     write(g2ml,'(A<lpath>,A2,I4,A5)')path,"00",l,".g2ml"
!$$     write(g2fml,'(A<lpath>,A2,I4,A6)')path,"00",l,".g2fml"
!$$  elseif ( l.le.99999 ) then
!$$     write(g2ml,'(A<lpath>,A1,I5,A5)')path,"0",l,".g2ml"
!$$     write(g2fml,'(A<lpath>,A1,I5,A6)')path,"0",l,".g2fml"
!$$  else
!$$     write(g2ml,'(A<lpath>,I6,A5)')path,l,".g2ml"
!$$     write(g2fml,'(A<lpath>,I6,A6)')path,l,".g2fml"
!$$  endif
!$$
   open(g2unit,file=g2ml,status='unknown')
   open(g2funit,file=g2fml,status='unknown')

   do ctyp=1,n2_typ
      write(g2unit,'(A3,1X,I2)')"typ",ctyp
      write(g2funit,'(A3,1X,I2)')"typ",ctyp
      do p1=1,n_rs
         do p2=1,n_eta
            write(g2unit,'(A3,1X,G16.8,1X,A4,1X,G16.8)')"eta",eta(p2),"Rs",rs(p1)
            write(g2funit,'(A3,1X,G16.8,1X,A4,1X,G16.8)')"eta",eta(p2),"Rs",rs(p1)
            write(g2funit,'(A12,1X,2(A15,1X))')"g2f_x","g2f_y","g2f_z"
            do j=1,im

               write(g2unit,'(G16.8)')g2(p2,p1,ctyp,j)
               write(g2funit,'(3(G16.8,1x))')g2f(p2,p1,1,ctyp,j),g2f(p2,p1,2,ctyp,j),g2f(p2,p1,3,ctyp,j)

            enddo
         enddo
      enddo
   enddo

   close(g2unit)
   close(g2funit)

return
end subroutine write_g2



integer function nitems(line)
character,intent(in):: line*(*)
integer i, n, toks

i = 1;
n = len_trim(line)
toks = 0
nitems = 0
do while(i <= n)
   do while(line(i:i) == ' ')
     i = i + 1
     if (n < i) return
   enddo
   toks = toks + 1
   nitems = toks
   do
     i = i + 1
     if (n < i) return
     if (line(i:i) == ' ') exit
   enddo
enddo
end function nitems
