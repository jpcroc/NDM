subroutine coord_soap(i)

 USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im,imm
  use tab_imm_m, ONLY : ityp,xp
  use ml_in_ndm_module, ONLY: massat,r_soap,ns_data

implicit none

  integer,intent(in) :: i
  integer :: j
  double precision,dimension(3) :: r_cm
  double precision :: sum_mass

  if (.not.allocated(r_soap)) allocate(r_soap(3,imm,ns_data))

  r_cm(:)=0d0
  sum_mass=0d0
  do j=1,im
     r_cm(:)=r_cm(:)+massat(ityp(j))*xp(:,j)
     sum_mass=sum_mass+massat(ityp(j))
  enddo
  r_cm(:)=r_cm(:)/sum_mass
  do j=1,im
     r_soap(:,j,i)=xp(:,j)-r_cm(:)
  enddo

return
end subroutine coord_soap

subroutine compute_soap(data_im,k_soap_out)

 USE T_kind_param_m, ONLY:  double
  use gen_com_m, ONLY: im,imm
  use tab_imm_m, ONLY : ityp
  use angular_functions, only : spherical_harm
  use ml_in_ndm_module, ONLY: ns_data,pi,r_cut,l_max,alpha_soap,kappa_soap,r_soap

 implicit none

  integer,dimension(ns_data),intent(in) :: data_im
  double precision,dimension(ns_data,ns_data),intent(out) :: k_soap_out

  integer :: i,k,ji,jk,l,m
  double precision :: r_ji,r_jk,fcut_ji,fcut_jk,arg_bessel
  double precision,dimension(0:l_max) :: ms_bessel
  double complex,dimension(-l_max:l_max,0:l_max) :: c
  double precision :: factor

  k_soap_out(:,:)=0d0
  do i=1,ns_data
     do k=1,i
        do ji=1,data_im(i)
           do jk=1,data_im(k)
          
              r_ji = dsqrt(sum(r_soap(1:3,ji,i)**2))
              r_jk = dsqrt(sum(r_soap(1:3,jk,k)**2))
   
              if (r_ji.gt.r_cut) cycle
              if (r_jk.gt.r_cut) cycle
              fcut_ji=0.5d0*(cos(pi*r_ji/r_cut)+1d0)
              fcut_jk=0.5d0*(cos(pi*r_ji/r_cut)+1d0)
   
   
              arg_bessel=alpha_soap*r_ji*r_jk
              if (arg_bessel==0) then
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
            
              c(:,:)=(0d0,0d0) 
              do l=0,l_max
                 do m=-l,l
                    c(m,l) = c(m,l) + 4d0 * pi * fcut_ji * fcut_jk * dexp(-alpha_soap*0.5*(r_ji**2+r_jk**2)) * &
                                      ms_bessel(l) * spherical_harm(l,m,r_soap(1:3,ji,i)) * conjg(spherical_harm(l,m,r_soap(1:3,jk,k)))
                 enddo
              enddo

           enddo !jk
        enddo !ji 
   
        do l=0,l_max
           do m=-l,l
              k_soap_out(k,i) = k_soap_out(k,i) + conjg(c(m,l)) * c(m,l)
           enddo
        enddo

     enddo !k
  enddo !i

  do i=1,ns_data
     do k=1,i
        k_soap_out(k,i)=(k_soap_out(k,i)/dsqrt(k_soap_out(k,k)*k_soap_out(i,i)))**kappa_soap
     enddo !k
  enddo !i

  
 deallocate(r_soap)

return
end subroutine compute_soap
