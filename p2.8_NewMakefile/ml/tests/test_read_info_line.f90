PROGRAM test_read_info_line_poscar
! 000  2  Fe  26  C 14  -11.48330883    5.03066499    0.00000000

implicit none

  integer               :: nb_elements
  character(len=3)      :: EFS_tag, element1, element2, element3, element4
  integer               :: mass1, mass2, mass3, mass4
  real(kind=kind(1.d0)) :: E_total, E_fit1, E_fit2
  logical               :: has_energy, has_force, has_stress
  integer               :: ios, inp
  

  open(file = 'test_info_line.poscar', unit = 50,  action = 'read',  IOSTAT = ios)
      
  inp = 50 
  read(inp, *) EFS_tag, nb_elements
  print *,  'FOUND:', EFS_tag
  
 
  if (EFS_tag(1:1) == '1') has_energy = .true.
  
  if (EFS_tag(1:2) == '1') has_force = .true. 
  
  if (EFS_tag(1:3) == '1') has_stress = .true. 
  
  
  print *, has_energy, has_force, has_stress
  print *, 'Found' , nb_elements, 'chemical elements'
  BACKSPACE(inp)
  
  if (nb_elements == 1) then
     read(inp,*) EFS_tag, nb_elements, element1, mass1, E_total, E_fit1, E_fit2
     print *,  element1, mass1
     
  elseif (nb_elements == 2) then
     read(inp,*) EFS_tag, nb_elements, element1, mass1, element2, mass2, E_total, E_fit1, E_fit2
     print *,  element1, mass1, element2, mass2
     
  elseif (nb_elements == 3) then  
     read(inp,*) EFS_tag, nb_elements, element1, mass1, element2, mass2, E_total, E_fit1, E_fit2
     print *,  element1, mass1, element2, mass2, element3, mass3
     
  elseif (nb_elements == 4) then   
     read(inp,*) EFS_tag, nb_elements, element1, mass1, element2, mass2, E_total, E_fit1, E_fit2
     print *,  element1, mass1, element2, mass2, element3, mass3, element4, mass4
     
  elseif (nb_elements > 4) then   
     print *, 'Too many chemical spieces in the file. Upgrade the reading DB subroutine'
     stop "fatal read_poscar_sasha"
     
  elseif (nb_elements < 1) then   
     print *, 'Number of chemical spieces can not be less than 1. Check you DB files'
     stop "fatal read_poscar_sasha"
 
  end if  

  
 
  close(inp)

       

  
  
  
END PROGRAM test_read_info_line_poscar
