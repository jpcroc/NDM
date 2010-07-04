
MODULE potential_eam
! This module defines the parameters for the Stillinger-Weber potential
  implicit none
  save
 
  integer, parameter :: NSPECIES = 1
!  integer, parameter :: MAXNEI = 300
  real(8), parameter :: ONE_THIRD = (1.0d0/ 3.0d0)
  real(8), parameter :: PI = 3.14159265358979D0
  real(8)            :: rcut,rcut2,Rc
  character(len=5)   :: name_element 

END MODULE  potential_eam


!subroutine init_potential()
!  use potential_eam
!  use defs
!  implicit none

  

!  if (NSPECIES == 2) then
!  else if (NSPECIES == 1 ) then
!    
!    call READ_PARAMPOT()
!    rcut=Rc+1.0
!    rcut2 = rcut*rcut
!    write(*,*) 'rcutrcut', Rc, rcut, rcut2
!  else
!     stop
!  endif
!end subroutine 




subroutine calcforce(natoms,types,pos,box,force,pot_energy)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE T_kind_param_m, ONLY:  double
      use gen_com_m
      use art_in_ndm_module
      use tab_imm_m
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  implicit none
!input output  decalarations for the ART parameters
  integer, intent(in) :: natoms
  integer, intent(in), dimension(natoms):: types
  real(8), intent(in), dimension(3*natoms), target:: pos
  real(8), intent(in) :: box(3)
  real(8), intent(out) :: pot_energy
  real(8), intent(out), dimension(3*natoms), target:: force



  real(8), dimension(:), pointer :: x, y, z 

  integer                             ::i 
  



    
! We first set-up pointers for the x, y, z components in the position and forces

!debug points for the people
   if (NATOMS /= im) then
    write(*,*) 'calcfo_ndm: NATOMS and im', NATOMS, im
   end if
    
    x => pos(1:NATOMS)
    y => pos(NATOMS+1:2*NATOMS)
    z => pos(2*NATOMS+1:3*NATOMS)
    do i=1,NATOMS 
      xp(1,i)=x(i)/angst
      xp(2,i)=y(i)/angst
      xp(3,i)=z(i)/angst
    enddo

        
       
!debug       write(*,*) 'in calfo_ndm....it_art, itab,itetabvois,ltabvois' 	
!debug       write(*,*) it_art, itab,itetabvois,ltabvois,mod(it_art,itetabvois) 	
	
          
!	 if (it_art==1) then
	  if (itab/=0) then
            if (mod(it_art,itab)==0) then
               call caltabt
            endif
         endif
!	end if 
         if (ltabvois.and.mod(it_art,itetabvois)==0) call caltabi 
         call calfo 


    do i=1,natoms
    force(i)         = fp(1,i)*erg2ev/angst
    force(i+NATOMS)  = fp(2,i)*erg2ev/angst
    force(i+2*NATOMS)= fp(3,i)*erg2ev/angst
    end do
    
     pot_energy=potist*erg2ev
     it_art=it_art+1



end subroutine 



  SUBROUTINE READ_PARAMPOT()
! In this version the present subroutine is obsolete
      return
   end subroutine   






