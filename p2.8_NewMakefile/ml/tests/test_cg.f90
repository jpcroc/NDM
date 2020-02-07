program test_cg
implicit none 
integer :: ntype 
double precision :: j_max
  
 j_max=4.0
 ntype=1
 call compute_cg_vector(j_max,ntype) 

end 

subroutine compute_cg_vector(j_max,ntype)
implicit none 
double precision, intent(in) :: j_max
integer, intent(in) ::  ntype

integer :: j, j1, j2, m, m1,m2, j0_max, j1_max, j2_max
double precision :: cg 

 if ((ntype==1) .and. dabs(int(j_max)-j_max).gt.0.1d0) then
    write(*,*) 'Probably you intend to use some SO3 or bi-SO3 -  related descriptor'
    write(*,*) 'in that case j_max can be only integer, now j_max', j_max
    stop
 end if 

 j0_max=int(j_max*ntype)
 j1_max=int(j_max*ntype)
 j2_max=int(j_max*ntype)
 
 if (allocated(cg_vector)) deallocate(cg_vector) ; allocate(cg_vector( 0:j1_max,-j1_max:j1_max, & 
                                                                       0:j2_max,-j2_max:j2_max, & 
                                                                       0:j0_max,-j0_max:j0_max )) 
 cg_vector(:,:,:,:,:,:)=0.d0
 do j =0,j0_max
 do j1=0,j1_max
 do j2=0,j2_max
 do m=-j, j,ntype
   do m1=-j1, j1,ntype
     do m2=-j2, j2,ntype
      if (m-m1-m2/=0) cycle
      if (j1+j2-j <0) cycle
      if (j-j1+j2 <0) cycle
      if (j-j2+j1 <0) cycle
      if (mod(j+j1+j2,ntype)==1) cycle
      cg_vector(j1,m1,j2,m2,j,m) =  cg(j1,m1,j2,m2,j,m, ntype)
      !debug write(*,'(7i4, e20.10)') j, j1,j2, m, m1,m2,m-m1-m2, cg(j1,m1,j2,m2,j,m, ntype)
     end do
   end do   
 end do
 end do
 end do
 end do

end 




   function cg(j1,m1,j2,m2,j,m,ntype)

      implicit none

      integer,intent(in) :: j1,j2,j,m1,m2,m,ntype
      double precision :: cg, wigner3j

      cg = (-1d0)**((j1-j2+m)/ntype) * dsqrt(2d0*j/ntype+1d0) * wigner3j(j1,m1,j2,m2,j,-m,ntype)

   end function cg


   function wigner3j(j1,m1,j2,m2,j,m,ntype)

      implicit none

      integer,intent(in) :: j1,j2,j,m1,m2,m, ntype
      integer :: t,tmin,tmax
      double precision :: factorial, wigner3j,triangle,coef,sum_coef

      coef = (-1d0)**((j1-j2-m)/ntype) * dsqrt(factorial((j1+m1)/ntype)*factorial((j1-m1)/ntype)* &
                                                        factorial((j2+m2)/ntype)*factorial((j2-m2)/ntype)* &
                                                        factorial((j+m)/ntype)*factorial((j-m)/ntype))

      triangle = dsqrt(factorial((j1+j2-j)/ntype)*factorial((j1-j2+j)/ntype)* &
                  factorial((-j1+j2+j)/ntype)/factorial((j1+j2+j+1)/ntype))

      sum_coef=0d0
      tmin = max((j2-j-m1)/ntype,(j1+m2-j)/ntype,0)
      tmax = min((j1+j2-j)/ntype,(j1-m1)/ntype,(j2+m2)/ntype)
      do t=tmin,tmax
         sum_coef = sum_coef + (-1d0)**t / (                   &
         factorial(t)*factorial((j1+j2-j)/ntype-t)*factorial((j1-m1)/ntype-t)* &
                      factorial((j2+m2)/ntype-t)*factorial((j-j2+m1)/ntype+t)*factorial((j-j1-m2)/ntype+t) )
      enddo

      wigner3j = coef * triangle * sum_coef

   end function wigner3j


   function factorial(n) result(res)

      implicit none

      integer,intent(in) :: n
      double precision :: res
      integer :: i

      if (n<0) then
         write(7,*) n
         write(7,*)  "Argument of factorial function should be positive"
         stop
       elseif (n==0) then
         res=1
      endif

      res=1d0
      do i=2,n
         res = res*i
      enddo

   end function factorial



