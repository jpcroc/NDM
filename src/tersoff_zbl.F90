module tersoff_zbl_mod
  USE arret_ndm_mod,only:arret_ndm
  USE zieg2_mod,only: zieg2
        USE spline_mod,only: cspline
        implicit none
        contains
subroutine tersoff_zbl
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m
  USE gen_com_m, ONLY:
  USE var_pot, ONLY:ngrid,catom,csive,ipo,lu_roff_pair,npair,ntyp,pot,pot_d,roff1,roff2,typ_pot_pair,typ_pot_pair
  USE force_tersoff_facteurs

  integer :: i,j,k,l,m,n,iti
  real(double) ::xsp(ngrid),ysp(ngrid),bsp(ngrid),csp(ngrid),dsp(ngrid)
  real(double):: ktor,ktorho
  real(double) :: rk,rhok,rk2



!interface
!
!subroutine zieg2(pot, pot_d, csive,ngrid, ntyp,npair, catom, roff1, roff2,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)
!  !-----------------------------------------------
!  !   M o d u l e s
!  !-----------------------------------------------
!  USE T_kind_param_m, ONLY:  double
!
!  implicit none
!  !-----------------------------------------------
!  !   D u m m y   A r g u m e n t s
!  !-----------------------------------------------
!  integer, dimension(:,:), allocatable  :: ipo                      ! indice des paires d'atomes
!  integer, allocatable:: typ_pot_pair(:) ! donne le type d'interaction de la paire
!  integer , intent(in) :: ngrid,ipotentiel
!  integer  :: ntyp
!  integer  :: npair
!  real(double) , intent(in) :: csive
!  real(double)  :: auxe= 23.06134575D-20 
!  real(double) , intent(inout) :: pot(4,npair,0:ngrid+1),pot_d(4,npair,0:ngrid+1)
!  real(double)  :: catom(ntyp)
!  real(double)  :: roff1(npair)
!  real(double)  :: roff2(npair)
!  logical :: lu_roff_pair(npair)
!
!end subroutine zieg2
!end interface


 write(6,*)'AJOUT ZBK a TERSOFF roff',roff1,roff2,csive
  ktor=csive
!C'est ça qui va pas !!! POT=0
     pot=0.0
     write(6,*)csive,ngrid,ntyp,npair,catom,roff1,roff2
!     write(6,*)'tersoff + ziegler = probablement plante voir force_tersiff_cel commente et initialisation de pot '
!     call arret_ndm

     call zieg2(pot,pot_d,csive,ngrid,ntyp,npair,catom,roff1,roff2,lu_roff_pair,ipotentiel,typ_pot_pair,ipo)
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

!           do k=1,ngrid            
!              rk=(k*ktor) 
!              write(648,'(I5, 4G15.7)')k,rk,pot(1,l,k)*erg2ev, fr(rk,Ater(l),lambda1(l))*erg2ev, (pot(1,l,k)+fr(rk,Ater(l),lambda1(l)))*erg2ev
!           end do


        case default
           write(6,*)'Tersoff zbl pas appliqu� � la paire , l, typ_pot_pair(l) = ',l, typ_pot_pair(l)
        end select
end do

end subroutine tersoff_zbl
end module tersoff_zbl_mod
