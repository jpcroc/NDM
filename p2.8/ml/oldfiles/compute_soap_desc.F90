subroutine compute_soap(i_start_at,i_final_at,local_soap_out,local_soapf_out)

 USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im,imm,A2cm,lperiod,bg,at,indi2
  use tab_imm_m, ONLY : ityp,iwmax2,fp,xp
  use var_pot, only : ntyp
  use angular_functions, only : spherical_harm,grad_spherical_harm
  use ml_in_ndm_module, ONLY: rangml,pi,r_cut,alpha_soap,l_max,n3_typ,t_typ,tconf,write_desc,w3_rho,weighted,massat, cg_vector, linvisible

  implicit none

  integer, intent (in) :: i_start_at,i_final_at
  double precision,dimension(0:l_max,n3_typ,imm),intent(out) :: local_soap_out
  double precision,dimension(0:l_max,3,n3_typ,imm),intent(out) :: local_soapf_out

  real(double), dimension(:,:), allocatable :: xpnp
  real(double), dimension(3) :: dxp_ji,dxp_jk
  integer :: i,j,iw,iw1,iw2,iz

  integer :: k,l,m,ctyp
  double precision :: r_ji,r_jk,fcut_ji,fcut_jk,arg_bessel
  double precision,dimension(0:l_max) :: ms_bessel
  double precision :: factor

  ALLOCATE(xpnp(3,imm))
  if (lperiod) then
   xpnp(:,:)=xp(:,:)
  else
   call notperiod(xp,xpnp)
  end if
  call cryst_to_cart (imm, xpnp, bg, -1)

  local_soap_out(:,:,:)=0d0
  local_soapf_out(:,:,:,:)=0d0
  if (allocated(t_typ)) deallocate(t_typ)
  allocate(t_typ(ntyp,ntyp,ntyp))
  call tconf(t_typ(:,:,:))

  if (i_start_at==1)  iw2=0
  if (i_start_at > 1) iw2=iwmax2(i_start_at-1)


  do j=i_start_at,i_final_at
     iw1=iw2+1
     iw2=iwmax2(j)

     do iw=iw1,iw2-1
        i=indi2(iw)
        if (j==i) cycle
        do iz=iw+1,iw2
           k=indi2(iz)
           if ((j==k).or.(i==k)) cycle
           if (weighted) then
              ctyp=1
              factor=w3_rho(massat(ityp(j)),massat(ityp(i)),massat(ityp(k)))
           else
              ctyp=t_typ(ityp(i),ityp(k),ityp(j))
              factor=1.d0
           endif

           dxp_ji(1:3) = xpnp(1:3,i) - xpnp(1:3,j)
           dxp_jk(1:3) = xpnp(1:3,k) - xpnp(1:3,j)
           WHERE ( (dxp_ji.GT.0.5d0).OR.(dxp_ji.LT.-0.5d0) )
              dxp_ji(1:3) = dxp_ji(1:3) - Dble(Nint(dxp_ji(1:3)))
           END WHERE
           WHERE ( (dxp_jk.GT.0.5d0).OR.(dxp_jk.LT.-0.5d0) )
              dxp_jk(1:3) = dxp_jk(1:3) - Dble(Nint(dxp_jk(1:3)))
           END WHERE

           dxp_ji = MatMul(at,dxp_ji)/A2cm
           dxp_jk = MatMul(at,dxp_jk)/A2cm
           r_ji = dsqrt(sum( dxp_ji(1:3)**2 ))
           r_jk = dsqrt(sum( dxp_jk(1:3)**2 ))

           fcut_ji = 0.5d0*(cos(pi*dsqrt(r_ji)/r_cut)+1d0)
           fcut_jk = 0.5d0*(cos(pi*dsqrt(r_jk)/r_cut)+1d0)

           arg_bessel=alpha_soap*r_ji*r_jk
           if (arg_bessel==0d0) then
              ms_bessel(0)=1d0
              ms_bessel(1:l_max)=0d0
           else
              do l=0,l_max
                 if (l==0) then
                    ms_bessel(0)=sinh(arg_bessel)/arg_bessel
                 elseif (l==1) then
                    ms_bessel(1)=cosh(arg_bessel)/arg_bessel - sinh(arg_bessel)/(arg_bessel)**2
                 else
                    ms_bessel(l)=ms_bessel(l-2)-(2*l-1)*ms_bessel(l-1)/arg_bessel
                 endif
              enddo
           endif

           do l=0,l_max
              do m=-l,l
                 local_soap_out(l,ctyp,j) = local_soap_out(l,ctyp,j) + factor * 4d0 * pi * fcut_ji * fcut_jk * dexp(-alpha_soap*0.5*(r_ji**2+r_jk**2)) * &
                                                     ms_bessel(l) * spherical_harm(l,m,dxp_ji(1:3)) * conjg(spherical_harm(l,m,dxp_jk(1:3)))
              enddo
           enddo

        enddo !iz
     enddo !iw

  enddo !j

  do i=1,ctyp
     do l=0,l_max
        local_soap_out(l,ctyp,:)=local_soap_out(l,ctyp,:)/dsqrt(dot_product(local_soap_out(l,ctyp,:),local_soap_out(l,ctyp,:)))
     enddo
  enddo

  if (write_desc) then
     call write_soap(local_soap_out,local_soapf_out)
  endif

 deallocate(xpnp)

return
end subroutine compute_soap


subroutine write_soap(soap,soapf)

  USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im,imm
  use ml_in_ndm_module, ONLY: n3_typ,path,l_max,l => i_poscar

  implicit none

  double precision,dimension(0:l_max,n3_typ,imm),intent(in) :: soap
  double precision,dimension(0:l_max,3,n3_typ,imm),intent(in) :: soapf
  character*70 :: soapml, soapfml
  integer :: j,p,ctyp,lpath,soapunit,soapfunit

  namelist /input_ml/ path
  lpath=len_trim(path)

  soapunit=45
  soapfunit=46

   if ( l.le.9 ) then
      write(soapml,'(A<lpath>,A5,I1,A9)')path,"00000",l,".soapml"
      write(soapfml,'(A<lpath>,A5,I1,A10)')path,"00000",l,".soapfml"
   elseif ( l.le.99 ) then
      write(soapml,'(A<lpath>,A4,I2,A9)')path,"0000",l,".soapml"
      write(soapfml,'(A<lpath>,A4,I2,A10)')path,"0000",l,".soapfml"
   elseif ( l.le.999 ) then
      write(soapml,'(A<lpath>,A3,I3,A9)')path,"000",l,".soapml"
      write(soapfml,'(A<lpath>,A3,I3,A10)')path,"000",l,".soapfml"
   elseif ( l.le.9999 ) then
      write(soapml,'(A<lpath>,A2,I4,A9)')path,"00",l,".soapml"
      write(soapfml,'(A<lpath>,A2,I4,A10)')path,"00",l,".soapfml"
   elseif ( l.le.99999 ) then
      write(soapml,'(A<lpath>,A1,I5,A9)')path,"0",l,".soapml"
      write(soapfml,'(A<lpath>,A1,I5,A10)')path,"0",l,".soapfml"
   else
      write(soapml,'(A<lpath>,I6,A9)')path,l,".soapml"
      write(soapfml,'(A<lpath>,I6,A10)')path,l,".soapfml"
   endif

   open(soapunit,file=soapml,status='unknown')
   open(soapfunit,file=soapfml,status='unknown')

   do ctyp=1,n3_typ
      do p=0,l_max
         write(soapunit,'(A3,1X,I2)')"typ",ctyp
         write(soapfunit,'(A3,1X,I2)')"typ",ctyp
         write(soapfml,'(A12,1X,2(A15,1X))')"soapf_x","soapf_y","soapf_z"
         do j=1,im

            write(soapunit,'(G16.8)')soap(l,ctyp,j)
            write(soapfunit,'(3(G16.8,1x))')soapf(l,1,ctyp,j),soapf(l,2,ctyp,j),soapf(l,3,ctyp,j)

         enddo
      enddo
   enddo

   close(soapunit)
   close(soapfunit)
return
end subroutine write_soap
