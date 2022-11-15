program hehe
  implicit none
  real*16 :: r, en,A,b,C(3:5),x,enm1,force
  integer::i,j,k,n,deux_n,ir
  real*16::fdeux_n
  real*16, parameter ::ev2erg=1.602d-12,rmax=12.
  integer,parameter::nmax=4000

  A=14.647429
  b=0.547926
! C(3)=1.461
! c(4)=14.11
! C(5)=183.6
 en=0
  do ir=1,nmax
     enm1=en
     r=ir*rmax/nmax
     en=A*exp(-r/b)
!     write(6,*)
!     write(6,*)'rep',r,en
!     do n=3,5
!        deux_n=2*n
!        x=b*r
!        en=en-fdeux_n(x,deux_n)*C(n)/r**deux_n
!        write(6,*)n,r**deux_n,fdeux_n(x,deux_n)
!     end do
     force=(en-enm1)/(rmax/nmax)
!     write(9,'(3D20.8)') r,en,force
     write(10,'(3D20.8)')r*1.0d-8,en*ev2erg,force*ev2erg/1d-8
     write(11,'(3D20.8)')r,en,force

  end do
!     write(6,*)'fichier 9 en a.u.'
     write(6,*)'fichier 10 en cgs pur'
     write(6,*)'fichier 11 ev Ang'

end program hehe

function fdeux_n(x,deux_n)
  implicit none
  integer :: deux_n
  real*16::x
  real*16::fdeux_n
  
  real*16::sum
  integer::fact
  integer ::k
  
  sum=0
     sum=sum+1
  do k=1,deux_n
     sum=sum+x**k/fact(k)
  end do
  fdeux_n=1-exp(-x)*sum
  return
end function fdeux_n


function fact(k)
  integer::k
  integer::fact

  integer::i
  fact=1
  if (k==0) return
  do i=1,k
     fact=fact*i
  end do
  return
end function fact
