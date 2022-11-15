module build_subdata_mod
        use read_poscar_mod
        use neighbours_mod
        !use math
        implicit none
        contains
!$-----------------------------------------------------
subroutine prepare_database
!$-----------------------------------------------------
use ml_in_ndm_module, only : rangml, db_file, debug, optimize_weights
implicit none
integer :: n_db_class
logical :: ok

if (debug) then
   if (rangml==0) write(6,'("ML: begin reading the databse .... in prepare_database")')
end if

inquire (file=db_file, exist=ok)
if (ok) then
  open (file=db_file, unit=50, action='read')
  !reading how many lines are in db_file
  call read_n_db_class(50,n_db_class)
  !allocate and fill the db_model object
  call read_db_file(50,n_db_class)
  ! in this is for the case of otimizing weigths ...
else
  if (rangml==0) write(6,*) 'ML: Fatal Error the database file is not there. The present name: ', db_file
  stop "read_database"
endif
close(50)


if (debug) then
   if (rangml==0) write(6,'("ML: end reading the databse .... in prepare_database")')
end if

call prepare_name_file_from_db

return
end subroutine prepare_database
!>-----------------------------------------------------

!$-----------------------------------------------------
subroutine read_n_db_class(inp,nlines)
! Read how many lines are in the inp file.
! The lines that have # in the firsrt character are ignored
! Input:
!          inp (the file number)
! Output:
!          nlines number of lines
!$-----------------------------------------------------
use ml_in_ndm_module, only: debug, rangml
implicit none
integer, intent(in) :: inp
integer, intent(out) :: nlines
character(len=80) :: ctmp
integer :: io

nlines = 0
DO
  read(inp,'(a)',iostat=io) ctmp
  !write(*,*) ctmp, ctmp(1:1), io
  if (io/=0) exit
  if (ctmp(1:1)=='#') cycle
  nlines = nlines + 1
END DO


if (debug) then
   if (rangml==0) write(6,'("ML: in read_n_db_class number of active lines in db file ", i5)') nlines
end if

rewind(inp)

return
end subroutine read_n_db_class
!>-----------------------------------------------------

!$-----------------------------------------------------
subroutine read_db_file(inp,n_db_class)
! Read the n_db_class lines of the db file
! Input:
!        inp (the file)
!        n_db_class  (no of lines)
! Output:
!        db_model object
!        iconf_data, iconf_data_test, iconf_data_train
!        allocate of config_real, config_desc, db_test, db_train
!$-----------------------------------------------------
use ml_in_ndm_module, only : rangml, debug, db_file, &
                             iconf_data, iconf_data_train, iconf_data_test, &
                             db_train, db_test, optimize_weights, selection_type, selection_type_first_start
use derived_types, only: db_model, config_real, config_desc
implicit none
integer, intent(in) :: inp, n_db_class
character(len=80) :: ctmp
character(len=1) :: c1,c2,c3
integer :: io,i,iconf_temp

if (allocated(db_model)) deallocate(db_model) ; allocate(db_model(n_db_class))
i=0
DO
  read(inp,'(a)',iostat=io) ctmp
  !debug write(*,*) ctmp(1:1), ctmp
  if (io/=0) exit
  if (ctmp(1:1)=='#') cycle
  backspace(inp)
  i = i + 1

  if (.not.(selection_type==selection_type_first_start)) then
    if (optimize_weights) then
      read(inp,*) db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec, c1,c2,c3, &
                  db_model(i)%w_e, db_model(i)%w_f, db_model(i)%w_s, &
                  db_model(i)%w_e_end, db_model(i)%w_f_end, db_model(i)%w_s_end
    else
      read(inp,*) db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec, c1,c2,c3, db_model(i)%w_e, db_model(i)%w_f, db_model(i)%w_s
    end if
    db_model(i)%no_start=1

    if (db_model(i)%no_total .lt. db_model(i)%no_selec) then
         if (rangml==0) then
           write(6,*) 'Some problems in read_db_file', db_file
           write(6,*) 'In the class', db_model(i)%class, 'and the KLM', db_model(i)%klm
           write(6,*)'The number of total files', db_model(i)%no_total,'is lower  than the no of selected files ', db_model(i)%no_selec
         end if
         stop "no of files no_total and no_selec in read_db_file"
    end if



  else

    if (optimize_weights) then
      read(inp,*) db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec, db_model(i)%no_start, c1,c2,c3, &
                  db_model(i)%w_e, db_model(i)%w_f, db_model(i)%w_s, &
                  db_model(i)%w_e_end, db_model(i)%w_f_end, db_model(i)%w_s_end
    else
      read(inp,*) db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec, db_model(i)%no_start, c1,c2,c3, db_model(i)%w_e, db_model(i)%w_f, db_model(i)%w_s
    end if

      if (db_model(i)%no_total .lt. (db_model(i)%no_selec + db_model(i)%no_start -1)) then
         if (rangml==0) then
           write(6,*) 'Some problems in read_db_file', db_file
           write(6,*) 'In the class', db_model(i)%class, 'and the KLM', db_model(i)%klm
           write(6,*)'The number of total files', db_model(i)%no_total,'is lower than the no of selected files  + starting no of file ', db_model(i)%no_selec, db_model(i)%no_start
         end if
         stop "no of files no_total and no_selec in read_db_file"
      end if

  end if ! selection_type


  if (c1=='T' .or. c1=='t') then
    db_model(i)%has_db_energy=.true.
    else if   (c1=='F' .or. c1=='f') then
    db_model(i)%has_db_energy=.false.
    else
      if (rangml==0) write(6,*) 'ML: Error in reading T/F energy in database inp file'
      stop 'error 1 in reading db_file.inp, subroutine read_db_file'
  end if


  if (c2=='T' .or. c2=='t') then
    db_model(i)%has_db_force=.true.
    else if   (c2=='F' .or. c2=='f') then
    db_model(i)%has_db_force=.false.
    else
      if (rangml==0) write(6,*) 'ML: Error in reading T/F force in database inp file'
      stop 'error 2 in reading db_file.inp, subroutine read_db_file'
  end if


  if (c3=='T' .or. c3=='t') then
    db_model(i)%has_db_stress=.true.
    else if   (c3=='F' .or. c3=='f') then
    db_model(i)%has_db_stress=.false.
    else
      if (rangml==0) write(6,*) 'ML: Error in reading T/F  stress in database inp file'
      stop 'error 3 in reading db_file.inp, subroutine read_db_file'
  end if

  !debug write(*,*) db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec, db_model(i)%has_db_energy, db_model(i)%has_db_force, db_model(i)%has_db_stress, db_model(i)%w_e, db_model(i)%w_f, db_model(i)%w_s
END DO
rewind(inp)

if (i .ne. n_db_class) then
  if (rangml==0) write(6,*) 'Problems in readind db_file',inp
  if (rangml==0) write(6,*) 'These values should be equal', n_db_class, i
  stop "read_db_file"
end if

if (debug) then
   if (rangml==0) write(6,*) 'ML: Reading DB file class, klm , no_total, no_selec'
   do i=1,n_db_class
      if (rangml==0) write(*,'(i3, "  ",(a)," ",(a),i6,i6)') i, db_model(i)%class, db_model(i)%klm, db_model(i)%no_total, db_model(i)%no_selec
      if (db_model(i)%no_total .lt. db_model(i)%no_selec) then
         if (rangml==0) then
           write(6,*) 'Some problems in readinf db_file', db_file
           write(6,*) 'In the class', db_model(i)%class, 'and the KLM', db_model(i)%klm
           write(6,*)'The number of total files', db_model(i)%no_total,'is lower  than the no of selected files ', db_model(i)%no_selec
         end if
         stop "no of files in read_db_file"
      end if
   end do
end if


iconf_temp=0
do i=1, size(db_model)
   iconf_temp = iconf_temp + db_model(i)%no_selec
end do
iconf_data_train = iconf_temp

iconf_temp=0
do i=1, size(db_model)
   iconf_temp = iconf_temp + db_model(i)%no_total
end do
iconf_data = iconf_temp
iconf_data_test = iconf_data - iconf_data_train

if (allocated(config_real)) deallocate(config_real) ; allocate(config_real(iconf_data))
if (allocated(config_desc)) deallocate(config_desc) ; allocate(config_desc(iconf_data))


if (allocated(db_train))    deallocate(db_train)    ; allocate(db_train(iconf_data_train))
if (allocated(db_test))     deallocate(db_test)     ; allocate(db_test(iconf_data_test))


return
end subroutine read_db_file
!>-----------------------------------------------------




!$-----------------------------------------------------
subroutine prepare_name_file_from_db
! Read some information in the database_file and set some vectors
!Input:
!       selection_type
!       db_model object
!Output:
!       db_train(1:iconf_data_train), db_test(1:iconf_data_test) - map the corresponding configuration between 1 to iconf_data
!       config_real%
!                  %class %klm %cnumber
!                  %filename %i_order %train
!$-----------------------------------------------------
use ml_in_ndm_module, only : rangml, debug, selection_type, &
                       selection_type_first, selection_type_last,  &
                       selection_type_random, selection_type_first_start, &
                      db_path,seed,iconf_data, iconf_data_train, iconf_data_test, &
                      db_train, db_test, optimize_weights
use derived_types, only: db_model, config_real
implicit none
integer :: i,i_conf, number,fnumber
character(len=180) :: ftemp
character(len=7)::cnumber
integer, dimension(:), allocatable :: anum
integer :: icount_temp, ii, i_train, i_test
!number of config to be analyzed
logical :: is_already_selected
number=1e+6
icount_temp=0

! preparing the name of the selected files for the traininf files:
icount_temp=0
do i=1,size(db_model)


   if (selection_type==selection_type_random) then
      if (allocated(anum)) deallocate(anum) ; allocate(anum(db_model(i)%no_selec))
      seed = seed + 1
      call rks2 (db_model(i)%no_total,db_model(i)%no_selec,seed,anum)
   end if
   !i_start_name_file=10
   do i_conf=db_model(i)%no_start,db_model(i)%no_selec+db_model(i)%no_start-1
     icount_temp = icount_temp + 1
     select case (selection_type)
       case (selection_type_first)
          fnumber = number + i_conf
       case (selection_type_last)
          fnumber = number + db_model(i)%no_total - db_model(i)%no_selec + i_conf
       case (selection_type_random)
          fnumber = number + anum(i_conf)
       case (selection_type_first_start)
          fnumber = number + i_conf
       case default
          if (rangml==0) write(*,*) 'ML error: Only 3 possible choices for selection_type 1, 2 or 3. Read the manual'
          stop "fatal read_db_file"
     end select

     write(cnumber,'(i7)') fnumber
     ftemp=trim(adjustl(db_path))//db_model(i)%class//"_"//db_model(i)%klm//"_"//cnumber(2:7)//'.poscar'
     config_real(icount_temp)%class=db_model(i)%class
     config_real(icount_temp)%klm=db_model(i)%klm
     config_real(icount_temp)%cnumber=cnumber(2:7)
     config_real(icount_temp)%db_line = i
     config_real(icount_temp)%no_file_in_db_line = fnumber - number
     config_real(icount_temp)%filename=ftemp
     config_real(icount_temp)%train=.true.
     config_real(icount_temp)%selected=.true.

     config_real(icount_temp)%w_e=db_model(i)%w_e
     config_real(icount_temp)%w_f=db_model(i)%w_f
     config_real(icount_temp)%w_s=db_model(i)%w_s

    if (optimize_weights) then
      config_real(icount_temp)%w_e_end=db_model(i)%w_e_end
      config_real(icount_temp)%w_f_end=db_model(i)%w_f_end
      config_real(icount_temp)%w_s_end=db_model(i)%w_s_end
    end if

     config_real(icount_temp)%has_energy=db_model(i)%has_db_energy
     config_real(icount_temp)%has_force=db_model(i)%has_db_force
     config_real(icount_temp)%has_stress=db_model(i)%has_db_stress

     if (debug) then
        if (rangml==0) write(*,'("ML: files of the database...in prepare_name_file_from_db: ",(a))') trim(ftemp)
     end if
     !call read_poscar_sasha(ftemp,icount_temp)
   end do ! iconf


   if (db_model(i)%no_selec == db_model(i)%no_total) cycle
   do i_conf=1,db_model(i)%no_total
     if (selection_type_random) then
      is_already_selected=.false.
       do ii=1,db_model(i)%no_selec
          if (i_conf == anum(ii)) is_already_selected=.true.
       end do
       if (is_already_selected) cycle
     else
       if (config_real(i_conf)%train) cycle
     end if
     icount_temp = icount_temp + 1
     fnumber = number + i_conf
     write(cnumber,'(i7)') fnumber
     ftemp=trim(adjustl(db_path))//db_model(i)%class//"_"//db_model(i)%klm//"_"//cnumber(2:7)//'.poscar'
     config_real(icount_temp)%class=db_model(i)%class
     config_real(icount_temp)%klm=db_model(i)%klm
     config_real(icount_temp)%cnumber=cnumber(2:7)

     config_real(icount_temp)%filename=ftemp
     config_real(icount_temp)%train=.false.

     config_real(icount_temp)%db_line = i
     config_real(icount_temp)%no_file_in_db_line = fnumber - number

     config_real(icount_temp)%w_e=db_model(i)%w_e
     config_real(icount_temp)%w_f=db_model(i)%w_f
     config_real(icount_temp)%w_s=db_model(i)%w_s

     if (optimize_weights) then
       config_real(icount_temp)%w_e_end=db_model(i)%w_e_end
       config_real(icount_temp)%w_f_end=db_model(i)%w_f_end
       config_real(icount_temp)%w_s_end=db_model(i)%w_s_end
     end if

     config_real(icount_temp)%has_energy=db_model(i)%has_db_energy
     config_real(icount_temp)%has_force=db_model(i)%has_db_force
     config_real(icount_temp)%has_stress=db_model(i)%has_db_stress

   end do

end do

if  (icount_temp /=iconf_data) then
   if (rangml==0) write(6,'("ML: Problems in the reading database icount_temp, iconf_data: ", 2i9)') &
                  icount_temp, iconf_data
  stop "inconsistencies in the DB reading: subroutine prepare_name_file_from_db"
end if

!Fill the db_test and db_train vectors.
i_train=0
i_test=0
do i=1,iconf_data
   if (config_real(i)%train) then
      i_train = i_train+1
      db_train(i_train) = i
    else
      i_test = i_test+1
      db_test(i_test)=i
   end if
end do

if  (i_train /=iconf_data_train) then
   if (rangml==0) write(6,'("ML: Problems in the reading database icount_temp, iconf_data_train:  ")') &
                  i_train, iconf_data_train
  stop "inconsistencies in the DB train reading: subroutine prepare_name_file_from_db"
end if


if  (i_test /=iconf_data_test) then
   if (rangml==0) write(6,'("ML: Problems in the reading database icount_temp, iconf_data_train:  ")') &
                  i_test, iconf_data_test
  stop "inconsistencies in the DB train reading: subroutine prepare_name_file_from_db"
end if


return
end subroutine prepare_name_file_from_db
!>-----------------------------------------------------




!$-------------------------------------------------------------------$!
subroutine read_poscar_sasha(name_file,ifile)                         !
! read the Milady's poscar defined by Sasha                           !
!Input:                                                               !
!          ifile, name_file                                           !
!Output:                                                              !
!         im, imm                                                     !
!         imm is resized if im > imm                                  !
!         config_real(ifile)%                                         !
!                           %pos_cart  %pos_crst                      !
!                           %cell, %nat, %itype, %ntypes              !
!                           %force %stress %spin                      !
!                           has_force, has_energy, has_stress         !
!$-------------------------------------------------------------------$!

use ml_in_ndm_module, only  : rangml,debug, im, imm, descriptor_type, &
                              descriptor_magnetic_sld, descriptor_magnetic_sld_afs, desc_forces, &
                              size_periodic_table, periodic_table_element, fix_weighted_for_element, fix_ch_elements, fix_no_of_elements, lfix_weight_auto, weighted,&
                              factor_weight_mass, linvisible, fix_ch_elements_invisible, fix_no_of_elements_invisible
use derived_types, only:config_real
use math
implicit none
character(len=*), intent(in) :: name_file
integer, intent(in) :: ifile
logical :: ok
integer :: inp
integer               :: nb_elements ! the same as nspecies???
character(len=3)      :: EFS_tag, element1, element2, element3, element4
character(len=2), dimension(:), allocatable :: element_poscar
integer, dimension(:), allocatable :: mass_poscar
integer               :: mass1, mass2, mass3, mass4
real(kind=kind(1.d0)) :: E_total, E_fit1, E_fit2
logical               :: has_energy, has_force, has_stress
real(kind=kind(1.d0)) :: alat
real(kind=kind(0.d0)), dimension(3,3) :: box, box_inv
character(len=80)   :: dummy,trimdummy
character(len=80)   :: line, header_poscar
integer :: im_local, icnt, ij, i_p, spin, i_pos_form,i,nitype ! ,nitems2, itest_there_is_a_number
integer, dimension(:), allocatable :: nspecies,ityp             ! nspecies is the same as nb_elements???
real(kind=kind(1.d0)), dimension(:,:), allocatable :: xp,xc,fp, l_spin
real(kind=kind(1.d0)) :: st(6), volume

inquire (file=trim(adjustl(name_file)), exist=ok)
if (ok) then
  open (file=trim(adjustl(name_file)), unit=50, action='read')
else
  if (rangml==0) write(6,*) 'ML: Fatal Error the poscar file is not there. The file: ', trim(adjustl(name_file))
  stop "fatal read_poscar_sasha"
endif

if (debug) then
     if (rangml==0) write(6,'("ML: Reading ..... in read_poscar_sasha ",(a))') name_file
end if
inp=50

has_energy=.false.
has_force=.false.
has_stress=.false.
read(inp, *) EFS_tag, nb_elements


if (EFS_tag(1:1) == '1') has_energy = .true.
if (EFS_tag(2:2) == '1') has_force = .true.
if (EFS_tag(3:3) == '1') has_stress = .true.



if (debug) then
  if (rangml==0) then
     if (rangml==0) print *, 'DEBUG: found EFS', EFS_tag
     if (rangml==0) print *, 'DEBUG:', has_energy, has_force, has_stress
     if (rangml==0) print *, 'DEBUG: found' , nb_elements, 'chemical elements'
  end if
end if
BACKSPACE(inp)

if (allocated(element_poscar)) deallocate(element_poscar) ; allocate(element_poscar(nb_elements))
if (allocated(mass_poscar))    deallocate(mass_poscar) ; allocate(mass_poscar(nb_elements))

read(inp,*) EFS_tag, nb_elements, ((element_poscar(ij),mass_poscar(ij)),ij=1,nb_elements), E_total, E_fit1, E_fit2

if (debug) then
  if (rangml==0) write(6,*) EFS_tag, nb_elements, ((element_poscar(ij),mass_poscar(ij)), ij=1,nb_elements), E_total, E_fit1, E_fit2
end if

if (allocated(config_real(ifile)%mass_per_type))  deallocate(config_real(ifile)%mass_per_type) ; allocate(config_real(ifile)%mass_per_type(nb_elements))
if (allocated(config_real(ifile)%weight_per_type))  deallocate(config_real(ifile)%weight_per_type) ; allocate(config_real(ifile)%weight_per_type(nb_elements))
if (allocated(config_real(ifile)%Z_per_type))  deallocate(config_real(ifile)%Z_per_type) ; allocate(config_real(ifile)%Z_per_type(nb_elements))
if (allocated(config_real(ifile)%fix_type_poscar_to_periodic)) deallocate(config_real(ifile)%fix_type_poscar_to_periodic) ; allocate(config_real(ifile)%fix_type_poscar_to_periodic(nb_elements))

!assign the appropiate mass ...
do ij=1,nb_elements
  do i_p=1,size_periodic_table
    if (element_poscar(ij)==periodic_table_element(i_p)%symbol) then
        config_real(ifile)%fix_type_poscar_to_periodic(ij)=periodic_table_element(i_p)%Z
        config_real(ifile)%mass_per_type(ij)=periodic_table_element(i_p)%mass
        config_real(ifile)%Z_per_type(ij)=periodic_table_element(i_p)%Z
    end if
  end do
  !debug  write(6,*) 'poscar', config_real(ifile)%Z_per_type(ij), config_real(ifile)%mass_per_type(ij)
end do

!assign the appropiate weigths ...
if (weighted) then
  do ij=1,nb_elements
    icnt=0
    do i_p=1, fix_no_of_elements
      if (fix_ch_elements(i_p)==element_poscar(ij)) then
        if (lfix_weight_auto) then
          config_real(ifile)%weight_per_type(ij)=config_real(ifile)%mass_per_type(ij)/factor_weight_mass
        else
          config_real(ifile)%weight_per_type(ij)=fix_weighted_for_element(i_p)
        end if
        icnt= icnt+1
      end if
    end do
    if (icnt==0) then
      if (rangml==0) write(6,*) "ML: We have detected in the database an element that is not in the input list:"
      if (rangml==0) write(6,*) "ML: The unknown element is ",  element_poscar(ij)," in the file ", name_file
      if (rangml==0) write(6,*) "ML: Please update the fix_no_of_elements and chemical_elements."
      stop 'fatal in read_poscar_sasha, unknown element'
    end if
  end do
 if (linvisible) then
  if (allocated(config_real(ifile)%invisible_per_type))  deallocate(config_real(ifile)%invisible_per_type) ; allocate(config_real(ifile)%invisible_per_type(nb_elements))
  config_real(ifile)%invisible_per_type(:)=.false.
  do ij=1,nb_elements
    icnt=0
    do i_p=1,fix_no_of_elements_invisible
      if (fix_ch_elements_invisible(i_p)==element_poscar(ij)) then
        config_real(ifile)%invisible_per_type(ij)=.true.
      end if
    end do
  end do
 end if

else
  do ij=1,nb_elements
          config_real(ifile)%weight_per_type(ij)=1.d0
  end do
end if
if ( (debug).and.(rangml==0) ) then
  do ij=1,nb_elements
    write(6,*) 'in poscar type, weights, Z_pertype, mass_per_type', ij,   config_real(ifile)%weight_per_type(ij), config_real(ifile)%Z_per_type(ij), config_real(ifile)%mass_per_type(ij)
  end do
end if
!stop 'stop in debug mode ... to be continued ... '



!C_debug if (nb_elements == 1) then
!C_debug    read(inp,*) EFS_tag, nb_elements, element1, mass1, E_total, E_fit1, E_fit2
!C_debug    if (debug) then
!C_debug      if (rangml==0) print *,  element1, mass1, E_total, E_fit1, E_fit2
!C_debug    end if
!C_debug elseif (nb_elements == 2) then
!C_debug    read(inp,*) EFS_tag, nb_elements, element1, mass1, element2, mass2, E_total, E_fit1, E_fit2
!C_debug    if (debug) then
!C_debug      if (rangml==0) print *,  'DEBUG: ', element1, mass1, element2, mass2, E_total, E_fit1, E_fit2
!C_debug    end if
!C_debug elseif (nb_elements == 3) then
!C_debug    read(inp,*) EFS_tag, nb_elements, element1, mass1, element2, mass2, E_total, E_fit1, E_fit2
!C_debug    if (debug) then
!C_debug       if (rangml==0) print *,  element1, mass1, element2, mass2, element3, mass3, E_total, E_fit1, E_fit2
!C_debug    end if
!C_debug elseif (nb_elements == 4) then
!C_debug    read(inp,*) EFS_tag, nb_elements, element1, mass1, element2, mass2, E_total, E_fit1, E_fit2
!C_debug    if (debug) then
!C_debug       if (rangml==0) print *,  element1, mass1, element2, mass2, element3, mass3, element4, mass4, E_total, E_fit1, E_fit2
!C_debug    end if

!C_debug elseif (nb_elements > 4) then
!C_debug    print *, 'Too many chemical spieces in the file. Upgrade the reading DB subroutine'
!C_debug    stop "fatal read_poscar_sasha"

if (nb_elements < 1) then
   print *, 'Number of chemical spieces can not be less than 1. Check yours DB files'
   stop "fatal read_poscar_sasha"
end if



! Lattice vector coordinates of the cell
read(inp,*) alat
read(inp,*) box(1:3,1)
read(inp,*) box(1:3,2)
read(inp,*) box(1:3,3)
box(:,:) = alat*box(:,:)

!if (rangml==0) then
!  write(*,*) box(:,:)
!end if

read(inp,'(a)') line
if (itest_there_is_a_number(line)==0) then
   !if (debug) then
   !    if (rangml==0) write(*,*) "ML: POSCAR type  VASP5 file"
   !end if
   nitype=nitems2(line)
 else
   !if (debug) then
   !     if (rangml==0) write(*,*) "ML: POSCAR type  VASP4 file"
   !end if
   nitype=nitems2(line)
   BACKSPACE(inp)
end if

if (nitype /= nb_elements) then
   if (rangml==0) then
       write(6,*) 'There are inconsistencies in the number of species in the poscar file nb_elements, nitype', nb_elements, nitype
       write(6,*) 'Check this file:  ', name_file
   end if
   stop 'species wrong defined: subroutine read_poscar_sasha'
end if

if (allocated(nspecies)) deallocate(nspecies); allocate(nspecies(nitype))

! Number of atoms in the cell
nspecies(:)=0
read(inp,*) (nspecies(i),i=1,nitype)

read(inp,'(a)') dummy
trimdummy=trim(adjustl(dummy))
if (.not.( (trimdummy(1:1).eq.'c').or.(trimdummy(1:1).eq.'C') .or. (trimdummy(1:1).eq.'d') .or. (trimdummy(1:1).eq.'D') ))  then
    if (rangml==0) then
      write(6,*) 'The POSCAR file is not in the good format. The first letter in that line is not D , d, C or c'
      write(6,*) trimdummy(1:1), " ", dummy
       write(6,*) 'Check this file:  ', name_file
    end if
    stop "vasp format read_poscar_sasha"
end if

if ((trimdummy(1:1).eq.'c').or.(trimdummy(1:1).eq.'C')) then
!this is cartesian format
i_pos_form=0
end if

if ((trimdummy(1:1).eq.'d').or.(trimdummy(1:1).eq.'D')) then
!this is crystalografic format
i_pos_form=1
end if

 ! Number of atoms in the simulation box
 im_local=SUM(nspecies(:))
if  (im_local.GT.imm) then
  imm=im_local
  if (debug) then
  if (rangml==0) then
    write(6,'("ML: WARNING the imm values was changed to a upper value of im in read_poscar_sasha ", i7)') imm
  end if
  end if
end if
im = im_local

if (allocated(xc))     deallocate(xc)     ; allocate(xc(1:3,1:im_local))
if (allocated(xp))     deallocate(xp)     ; allocate(xp(1:3,1:im_local))
if (allocated(fp))     deallocate(fp)     ; allocate(fp(1:3,1:im_local))
if (allocated(l_spin)) deallocate(l_spin) ; allocate(l_spin(1:3,1:im_local))
if (allocated(ityp))   deallocate(ityp)   ; allocate(ityp(1:im_local))
icnt=0
do ij=1,nitype
   do i=1, nspecies(ij)
     icnt=icnt+1
     ityp(icnt)=ij
     read(inp,*) xc(1:3,icnt)
   end do
end do


 ! Atom real coordinates
if (i_pos_form==1) then
    ! inthat case crystalografic coord are read
    xp(1:3,1:im_local) = MatMul( box(1:3,1:3), xc(1:3,1:im_local) )
  else if (i_pos_form==0) then
    ! in that case cartesian coordinates are read
    call matinv_gen(box,box_inv)
    xp(1:3,1:im_local) = xc(1:3,1:im_local)
    xc(1:3,1:im_local) = MatMul(box_inv(1:3,1:3), xp(1:3,1:im_local))
end if
if (desc_forces) then
  read(inp,*)
  do i=1,im_local
   read(inp,*) fp(1:3,i)
  end do

  read(inp,*)
  read(inp,*) st(1:6)
end if

!debug read(inp,*)
!debug read(inp,*) spin
spin=0

if ((descriptor_type == descriptor_magnetic_sld).or.(descriptor_type == descriptor_magnetic_sld_afs)) then
  do i=1,im_local
    read(inp,*) l_spin(1:3,i)
  end do
end if




if (allocated(config_real(ifile)%itype)) deallocate(config_real(ifile)%itype);  allocate(config_real(ifile)%itype(im_local))

if (allocated(config_real(ifile)%pos_cart)) deallocate(config_real(ifile)%pos_cart) ;   allocate(config_real(ifile)%pos_cart(3,im_local))
if (allocated(config_real(ifile)%pos_crst)) deallocate(config_real(ifile)%pos_crst) ;   allocate(config_real(ifile)%pos_crst(3,im_local))
if (allocated(config_real(ifile)%force))    deallocate(config_real(ifile)%force) ;   allocate(config_real(ifile)%force(3,im_local))
if (allocated(config_real(ifile)%atomic_spin))    deallocate(config_real(ifile)%atomic_spin) ;   allocate(config_real(ifile)%atomic_spin(3,im_local))
config_real(ifile)%ntypes=nitype
!C_debug if (nitype==1) then
!C_debug   config_real(ifile)%mass_per_type(1)=mass1
!C_debug elseif (nitype==2) then
!C_debug   config_real(ifile)%mass_per_type(1)=mass1
!C_debug   config_real(ifile)%mass_per_type(2)=mass2
!C_debug elseif (nitype==3) then
!C_debug   config_real(ifile)%mass_per_type(1)=mass1
!C_debug   config_real(ifile)%mass_per_type(2)=mass2
!C_debug   config_real(ifile)%mass_per_type(3)=mass3
!C_debug elseif (nitype==4) then
!C_debug   config_real(ifile)%mass_per_type(1)=mass1
!C_debug   config_real(ifile)%mass_per_type(2)=mass2
!C_debug   config_real(ifile)%mass_per_type(3)=mass3
!C_debug   config_real(ifile)%mass_per_type(4)=mass4
!C_debug elseif (nitype > 4) then
!C_debug    if (rangml==0) write(6,*) 'Fatal error number of types is lower than 5'
!C_debug    stop 'fatal in nitype in read_poscar_sasha'
!C_debug end if


config_real(ifile)%nat=im_local
config_real(ifile)%itype(:)=ityp
config_real(ifile)%pos_crst(1:3,1:im_local)=xc(1:3,1:im_local)
config_real(ifile)%pos_cart(1:3,1:im_local)=xp(1:3,1:im_local)
config_real(ifile)%force=fp
config_real(ifile)%energy(1)=E_total
config_real(ifile)%energy(2)=E_fit1
config_real(ifile)%energy(3)=E_fit2

config_real(ifile)%cell=box
config_real(ifile)%stress(1:6)=st(1:6)
config_real(ifile)%spin=spin
config_real(ifile)%atomic_spin=l_spin


config_real(ifile)%has_energy=has_energy.and.config_real(ifile)%has_energy
config_real(ifile)%has_force=has_force.and.config_real(ifile)%has_force
config_real(ifile)%has_stress=has_stress.and.config_real(ifile)%has_stress


call calc_volume(config_real(ifile)%cell(:,1),config_real(ifile)%cell(:,2),config_real(ifile)%cell(:,3),volume)
config_real(ifile)%volume = volume


close(inp)
return
end subroutine read_poscar_sasha
!>-----------------------------------------------------






!old subroutine from Wesley ... .....
subroutine buildsubdata(nd_data)

use ml_in_ndm_module, only : rangml,selection_type,pref,ns_data,seed,kelem

implicit none

integer,intent(in) :: nd_data

integer :: i,j
integer,dimension(ns_data,kelem) :: sample
double precision :: y_target

namelist /input_ml/ selection_type,pref,ns_data,kelem,seed

select case(selection_type)
  case(1)
     if (rangml==0) write(6,*)"ML: Selects first",ns_data,"elements of database"
     do i=1,ns_data
        call read_poscar(pref,i,y_target)
     enddo
  case(2)
     if (rangml==0) write(6,*)"ML: Selects last",ns_data,"elements of database"
     do i=nd_data-(ns_data-1),nd_data
        call read_poscar(pref,i,y_target)
     enddo
  case (3)
     if (rangml==0) write(6,*)"ML: Random selection of",ns_data,"subset of",kelem,"elements of database"
     if (ns_data .gt. kelem) then
          if (rangml==0) write(6,*) "ML error: ns_data=", ns_data,"should be lower than kelem=",kelem,"and obviously is not the case"
          stop
     end if


     call combinations(nd_data,kelem,ns_data,seed,sample)
     do i=1,ns_data
        do j=1,kelem
           call read_poscar(pref,sample(i,j),y_target)
        enddo
     enddo
   case default
      if (rangml==0) write(*,*) 'ML error: Only 3 possible choices for selection_type 1, 2 or 3. Read the manual'
      stop

end select

return
end subroutine buildsubdata


subroutine get_number_of_files(nfiles)

use ml_in_ndm_module, only : path,lpath,rangml

implicit none

integer,intent(out) :: nfiles
!oldwesley character(len=120) :: command
lpath=len_trim(path)

if (rangml==0) then
!oldwesley  write(command,*)"ls $(echo '",path,"*.poscar' | sed 's/ *//g') | wc -l > nfiles.txt"
!oldwesley call system(command)
!oldwesley open(31,file='nfiles.txt',action="read")
!oldwesley  read(31,*)nfiles
!oldwesley close(31)
!oldwesley call system('rm nfiles.txt')
nfiles=1
end if


return
end subroutine get_number_of_files


subroutine combinations(n,k,n_sample,seedin,b)

implicit none

integer(kind=4),intent(in) :: n,k,seedin,n_sample
integer,dimension(n_sample,k),intent(out) :: b

integer,dimension(k) :: a
integer(kind=4) :: i,j,c,seed
logical :: test

seed=seedin
call rks2 (n,k,seed,a)
b(1,:)=a(:)
i=2
c=0
do while (i-1.ne.n_sample)
   c=c+1
   seed=seed+c
   call rks2 (n,k,seed,a)

   test = .true.
   do j = 1,i-1
      if (all(a(:) == b(j,:))) then
         test = .false.
         exit
      endif
   end do

   if (test) then
      b(i,:)=a(:)
      i=i+1
   endif
enddo

return
end subroutine combinations


subroutine rks2 ( n, k, seed, a )

use ml_in_ndm_module, only : rangml

implicit none

integer (kind = 4) :: k,c1,c2,i,k0,n,seed
integer (kind = 4),dimension(k) :: a
real (kind = 8) :: r ! ,r8_uniform_01

if ( k < 0 .or. n < k ) then
  if (rangml==0) write ( *, '(a)' ) ''
  if (rangml==0) write ( *, '(a)' ) 'KSUB_RANDOM2 - Fatal error!'
  if (rangml==0) write ( *, '(a,i8)' ) '  N = ', n
  if (rangml==0) write ( *, '(a,i8)' ) '  K = ', k
  if (rangml==0) write ( *, '(a)' ) '  but 0 <= K <= N is required!'
  stop 1
end if

if ( k == 0 ) then
  return
end if

c1 = k
c2 = n
k0 = 0
i = 0

do i = 1,n

  r = r8_uniform_01 ( seed )

  if ( real ( c2, kind = 8 ) * r <= real ( c1, kind = 8 ) ) then

    c1 = c1 - 1
    k0 = k0 + 1
    a(k0) = i

    if ( c1 <= 0 ) then
      exit
    end if

  end if

  c2 = c2 - 1

end do
return
end subroutine rks2

function r8_uniform_01 ( seed )

use ml_in_ndm_module, only : rangml

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = 8 ) value.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^1 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Pierre LEcuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  !integer (kind = 4) :: i4_huge,k,seed
  integer (kind = 4) :: k,seed
  real (kind = 8) :: r8_uniform_01

  if ( seed == 0 ) then
    if (rangml==0) write ( *, '(a)' ) ''
    if (rangml==0) write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    if (rangml==0) write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge ( )
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function r8_uniform_01

function i4_huge ( )

!*****************************************************************************80
!
!! I4_HUGE returns a "huge" I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) I4_HUGE, a "huge" integer.
!
  implicit none

  integer ( kind = 4 ) i4_huge

  i4_huge = 2147483647

  return
end function i4_huge
end module
