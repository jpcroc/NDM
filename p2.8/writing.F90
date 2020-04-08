module writing_mod
  USE cryst_to_cart_mod
  USE gen_com_m, ONLY:imm,at,bg
  implicit none
contains
  subroutine writing (file_num, file_name, grandeur)

    ! *****************************************************
    ! THIS PROGRAM DETECT THE ATOMS THAT BELONG TO THE
    ! CRYSTAL SURFACE (case ibound = 1 or 2)
    ! *****************************************************
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double

    USE tab_imm_m

    implicit none	
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   P a r a m e t e r s
    !-----------------------------------------------
    !-----------------------------------------------
    !   L o c a l   V a r i a b l e s
    !-----------------------------------------------
    integer, intent(in) 		:: file_num
    character (LEN=*), intent(in)   :: file_name
    real (double), dimension(imm), intent(in) :: grandeur
    integer :: i
    integer :: iform
    !-----------------------------------------------
    !
    !-----------------------------------------------

    iform = 3	!(iform = 1 => CL périodiques selon X)
    !(iform = 2 => CL libres selon X; boite "rectangulaire") 

    call cryst_to_cart (imm, xp, bg, -1) !cart vers cryst

    ! writing of the output file (implemented for an unique atom type)
    Open (unit=file_num, file=file_name, status='unknown',action='write')     

    Write (file_num, '(a, i5)') 	'Number of particles = ', imm
    Write (file_num, '(a)') 	'A = 1.000 Angstrom (basic length-scale)'
    Write (file_num, '(a)') 	'# Unit cell vector #1'
    if      (iform==1) then	
       Write (file_num, '(a)') 	'H0(1,1) =    104.9104210000000     A'
       Write (file_num, '(a)') 	'H0(1,2) =    0.000000000000000     A'
       Write (file_num, '(a)') 	'H0(1,3) =    0.000000000000000     A'
       Write (file_num, '(a)') 	'# Unit cell vector #2'
       Write (file_num, '(a)') 	'H0(2,1) =    78.68281600000000     A'
       Write (file_num, '(a)') 	'H0(2,2) =    81.76958100000000     A'
       Write (file_num, '(a)') 	'H0(2,3) =    0.000000000000000     A'
    else if (iform==2) then
       Write (file_num, '(a)') 	'H0(1,1) =    133.4634209385322     A'
       Write (file_num, '(a)') 	'H0(1,2) =    0.000000000000000     A'
       Write (file_num, '(a)') 	'H0(1,3) =    0.000000000000000     A'
       Write (file_num, '(a)') 	'# Unit cell vector #2'
       Write (file_num, '(a)') 	'H0(2,1) =    0.000000000000000     A'
       Write (file_num, '(a)') 	'H0(2,2) =    83.06605379269264     A'
       Write (file_num, '(a)') 	'H0(2,3) =    0.000000000000000     A'
    else if (iform==3) then
       Write(file_num,'(a,f25.12,a)')'H0(1,1) = ', at(1,1)*1e+8, '  A'
       Write(file_num,'(a,f25.12,a)')'H0(1,2) = ', at(2,1)*1e+8, '  A'
       Write(file_num,'(a,f25.12,a)')'H0(1,3) = ', at(3,1)*1e+8, '  A'
       Write(file_num,'(a)')'# Unit cell vector #2'
       Write(file_num,'(a,f25.12,a)')'H0(2,1) = ', at(1,2)*1e+8, '  A'
       Write(file_num,'(a,f25.12,a)')'H0(2,2) = ', at(2,2)*1e+8, '  A'
       Write(file_num,'(a,f25.15,a)')'H0(2,3) = ', at(3,2)*1e+8, '  A'
       Write(file_num,'(a)')'# Unit cell vector #3'
       Write(file_num,'(a,f25.12,a)')'H0(3,1) = ', at(1,3)*1e+8, '  A'
       Write(file_num,'(a,f25.12,a)')'H0(3,2) = ', at(2,3)*1e+8, '  A'
       Write(file_num,'(a,f25.12,a)')'H0(3,3) = ', at(3,3)*1e+8, '  A'
    end if



    Write (file_num, '(a)') 	'# Unit cell vector #3'
    Write (file_num, '(a)') 	'H0(3,1) =    0.000000000000000     A'
    Write (file_num, '(a)') 	'H0(3,2) =    0.000000000000000     A'
    Write (file_num, '(a)') 	'H0(3,3) =    24.72762300000000     A'
    Write (file_num, '(a)') 	'.NO_VELOCITY.'
    Write (file_num, '(a)') 	'entry_count = 4'
    Write (file_num, '(2a)') 	'auxiliary[0] = ', file_name
    Write (file_num, '(a)') 	'0.000'
    Write (file_num, '(a)') 	'Fe '

    Do i=1,imm
       Write(file_num,'(e24.16,3x,e24.16,3x,e24.16,3x,e24.16)') xp(1,i),xp(2,i),xp(3,i),grandeur(i)
    End do

    call cryst_to_cart (imm, xp, at, 1)  !cryst vers cart

  end subroutine writing
end module writing_mod
