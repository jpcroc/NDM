subroutine tersoff_zbl
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m
  use gen_com_m
  use var_pot

  integer :: i,j,k,l,m,n,iti
  real(double) ::xsp(ngrid),ysp(ngrid),bsp(ngrid),csp(ngrid),dsp(ngrid)
  real(double):: ktor,ktorho
  real(double) :: rk,rhok,rk2

  ktor=(maxval(roff2)+1.)/ngrid
!C'est ça qui va pas !!! POT=0
     pot=0.0
     write(6,*)csive,ngrid,ntyp,npair,catom,roff1,roff2
!     write(6,*)'tersoff + ziegler = probablement plante voir force_tersiff_cel commente et initialisation de pot '
!     stop

     call zieg2(pot,pot_d,csive,ngrid,ntyp,npair,catom,roff1,roff2,lu_roff_pair)
!    write(6,*)pot

  do l=1,npair
     select case (typ_pot_pair(l))
        case(13,14,15)
           do k=1,ngrid            
              rk=(k*ktor) 
              xsp(k)=rk
              ysp(1:ngrid)=pot(1,l,1:ngrid)
           end do
           call cspline (ngrid,xsp,ysp,bsp,csp,dsp)
           pot(1,l,1:ngrid)=ysp(1:ngrid)
           pot(2,l,1:ngrid)=bsp(1:ngrid)
           pot(3,l,1:ngrid)=csp(1:ngrid)
           pot(4,l,1:ngrid)=dsp(1:ngrid)
        case default
           write(6,*)'PB avec tersoff zbl, l, typ_pot_pair(l) = ',l, typ_pot_pair(l)
        end select
end do

end subroutine tersoff_zbl
