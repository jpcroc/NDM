module config2data_mod
  USE gen_com_m, only:uwrt,lwrt  
   USE arret_ndm_mod,only:arret_ndm  
  USE T_kind_param_m, ONLY:  double
  USE Mat_utils_mod,only: Matinv,is_upper_triangular,convert_cell
  USE gen_com_m, ONLY : position_conversion_lammps,energy_conversion_lammps

    implicit none
    real(double), dimension(3,3)::passage,passage_inv
!  real(double), dimension(3,3)::old_passage
  real(double), dimension(3,3)::atprec
contains
!FOR ALL PROCESSES: XP MUST CONTAIN moTHE ACTUAL COORDINATES   
  subroutine config2data (imm,im,xp,ityp,at,ntyp,lwrite,filename)
    USE T_kind_param_m, ONLY:  double
    USE var_pot, ONLY:q,ipotentiel

    implicit none
    integer,intent(in)::imm,im,ntyp
    real(double),intent(in)::xp(3,imm),at(3,3)
    integer,intent(in)::ityp(imm)
    logical,intent(in)::lwrite
    character(*),optional :: filename
    character*80::file

    logical upper
    real(double)::xhi,yhi,zhi,xy,xz,yz,xlo,ylo,zlo,QTOT
    integer::ic,i
    real(double), dimension(3) :: tmp_coord_i,new_tmp_coord_i
    real(double), dimension(3,3)::at_lammps
    file='conf.lmp'
    if (present(filename))file=filename

    !  if ((at(2,1).ne.0).or.(at(3,1).ne.0).or.(at(3,2).ne.0))then

    atprec=at
    if (lwrite) then 
       open(63,file=file,status='unknown')
       write(63,*)
    end if
    call is_upper_triangular(at,upper)
    if (.not.upper) then
       call convert_cell (at,at_lammps,passage)
       call matinv(passage, passage_inv)
!!$       do ic=1,3
!!$          write(uwrt,*)passage(:,ic)
!!$       end do
       !  write(uwrt,*)
    else
       at_lammps=at
       passage(:,:)=0
       do ic=1,3
          passage(ic,ic)=1
       end do
       passage_inv(:,:)=passage(:,:)
    end if

    !    old_passage=passage
    !building  lammps header
    xlo = 0.d0
    ylo = 0.d0
    zlo = 0.d0
    xhi = at_lammps(1,1)/position_conversion_lammps
    yhi = at_lammps(2,2)/position_conversion_lammps
    zhi = at_lammps(3,3)/position_conversion_lammps
    xy = at_lammps(1,2)/position_conversion_lammps
    xz = at_lammps(1,3)/position_conversion_lammps
    yz = at_lammps(2,3)/position_conversion_lammps

    !  xhi=at(1,1)/position_conversion_lammps
    !  yhi=dsqrt(at(1,2)**2+at(2,2)**2)/position_conversion_lammps
    !  zhi=dsqrt(at(1,3)**2+at(2,3)**2+at(3,3)**2)/position_conversion_lammps
    !  xz=at(1,1)*at(1,3)/(xhi*zhi*position_conversion_lammps*position_conversion_lammps)
    !  xy=at(1,1)*at(1,2)/(xhi*yhi*position_conversion_lammps*position_conversion_lammps)
    !  yz=(at(1,2)*at(1,3)+at(2,2)*at(2,3))/(yhi*zhi*position_conversion_lammps*position_conversion_lammps)
    if (lwrite) then 
       if(IM.lt.10)then
          write(63,"(I1,A)") IM,' atoms'
       elseif((IM.lt.100).and.(IM.gt.10))then
          write(63,"(I2,A)") IM,' atoms'
       elseif((IM.lt.1000).and.(IM.gt.100))then
          write(63,"(I3,A)") IM,' atoms' 
       elseif((IM.lt.10000).and.(IM.gt.1000))then
          write(63,"(I4,A)") IM,' atoms'       
       elseif((IM.lt.100000).and.(IM.gt.10000))then
          write(63,"(I5,A)") IM,' atoms'         
       elseif((IM.lt.1000000).and.(IM.gt.100000))then
          write(63,"(I6,A)") IM,' atoms'         
       elseif((IM.lt.10000000).and.(IM.gt.1000000))then
          write(63,"(I7,A)") IM,' atoms'        
       elseif((IM.lt.100000000).and.(IM.gt.10000000))then
          write(63,"(I8,A)") IM,' atoms'         
       else
          write(uwrt,*)'add format'
          call arret_ndm (.true.)
       endif
       write(uwrt,*)
       write(uwrt,*) " a = (xhi-xlo,0,0); b = (xy,yhi-ylo,0); c = (xz,yz,zhi-zlo). "
       write(uwrt,'(f22.16,a,f22.16,a)')xlo,' ',xhi,' xlo xhi'
       write(uwrt,'(f22.16,a,f22.16,a)')ylo,' ',yhi,' ylo yhi'
       write(uwrt,'(f22.16,a,f22.16,a)')zlo,' ',zhi,' zlo zhi'
       write(uwrt,'(f22.16,a,f22.16,a,f22.16,a)')xy,' ',xz,' ',yz,' xy xz yz'
       write(uwrt,*)


       write(63,"(I1,A)") ntyp, ' atom types' 
       write(63,*) 
       write(63,'(f22.16,a,f22.16,a)')xlo,' ',xhi,' xlo xhi'
       write(63,'(f22.16,a,f22.16,a)')ylo,' ',yhi,' ylo yhi'
       write(63,'(f22.16,a,f22.16,a)')zlo,' ',zhi,' zlo zhi'
       write(63,'(f22.16,a,f22.16,a,f22.16,a)')xy,' ',xz,' ',yz,' xy xz yz'
       write(63,*)
       write(63,"(A)")'Atoms'
       write(63,*)
    end if

    select case (ipotentiel)

    case(-10)
       do i=1,im
          tmp_coord_i = xp(:,i)
          new_tmp_coord_i = matmul(passage,tmp_coord_i)/position_conversion_lammps
          if (lwrite)  write(63,'(I8,a,I2,a,f22.15,a,f22.15,a,f22.15)')i,' ',ityp(i),' ',new_tmp_coord_i(1),' '&
               &,new_tmp_coord_i(2),' ',new_tmp_coord_i(3)

          !        write(63,"(I8,I6,F21.12,F20.12,F20.12)") i,ityp(i),&
          !             & xp(1,i)/position_conversion_lammps,xp(2,i)/position_conversion_lammps,xp(3,i)/position_conversion_lammps
       end do
    case(-11)
       QTOT=0
       do i=1,im
          QTOT=qtot+Q(ITYP(I))
       end do
       do i=1,im
          tmp_coord_i = xp(:,i)
          new_tmp_coord_i = matmul(passage,tmp_coord_i)/position_conversion_lammps
          if (lwrite)    write(63,'(I8,a,I2,a,F20.12,a,f22.15,a,f22.15,a,f22.15)')i,' ',ityp(i),' ',q(ityp(i))-QTOT/im,&
               &' ',new_tmp_coord_i(1),' ',new_tmp_coord_i(2),' ',new_tmp_coord_i(3)
          !        write(63,"(I8,I6,F20.12,F20.12,F20.12,F20.12)") i,ityp(i),&
          !             & q(ityp(i)), xp(1,i)/position_conversion_lammps,xp(2,i)/position_conversion_lammps,xp(3,i)/position_conversion_lammps
       end do
    end select
    if (lwrite)then
       call flush(63)
       close (63)
    end if
  end subroutine config2data


  subroutine flmp2fndm (fp,force_lammps,im,imm)

    real(double),intent(in)::force_lammps(3*im)
    real(double),intent(out)::fp(3,imm)
    integer,intent(in)::im,imm
    real(double), dimension(3) :: tmp_coord_i,new_tmp_coord_i
    integer::ic,ip,i
    ip=0
    do i=1,im
       do iC=1,3
          ip=ip+1
          tmp_coord_i(ic) = force_lammps(ip)
       end do
       new_tmp_coord_i = matmul(passage_inv,tmp_coord_i)*energy_conversion_lammps/position_conversion_lammps
       do ic=1,3
          fp(ic,i)= new_tmp_coord_i(ic)
       end do
    end do

  end subroutine flmp2fndm
end module config2data_mod
