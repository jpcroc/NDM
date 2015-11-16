subroutine init_neb_reaction
!   NDM part ...
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY : imm
! MAB part  ...
    USE mab_in_ndm_module, ONLY : nimage_neb,nhisto
    USE reaction_neb_module, ONLY: neb_images,atoms_on_spline,lambda_images, &
                                   llambda,coord_neb,  &
                                   npoints,lambda_min,lambda_max, &
                                   delta_lambda, delta_neb, &
                                   force_defect, free_energy_force_defect
    implicit none
!local definitions 

    integer :: ia,ix,ip

!object allocation: the full image of the real NEB path
    allocate (neb_images(nimage_neb))
    do ip=1,nimage_neb
          allocate(neb_images(ip)%xp_neb(3,imm))
          allocate(neb_images(ip)%fp_neb(3,imm))
          allocate(neb_images(ip)%ityp(imm))
          allocate(neb_images(ip)%ielat(imm))
          allocate(neb_images(ip)%iwmax(imm))
    end do

!object allocation: the coefficient of the spline for each atom (1 to im) and direction (1 to 3)

  allocate(atoms_on_spline(3,imm))
  do ia=1,imm 
    do ix=1,3
      allocate(atoms_on_spline(ix,ia)%b(nimage_neb))
      allocate(atoms_on_spline(ix,ia)%c(nimage_neb))
      allocate(atoms_on_spline(ix,ia)%d(nimage_neb))
   end do
 end do


!object allocation: the  interpolated path on grid of dimension npoints 
! this should be changed ... when we will patch with 3 images back and forth this should be changed ....
  npoints=nhisto

 allocate(lambda_images(0:npoints))
 do ip=0,npoints  
    allocate(lambda_images(ip)%xp_neb(3,imm))
    allocate(lambda_images(ip)%dxp_neb(3,imm))
 end do 




! this should be changed ... when we will patch with 3 images back and forth this should be changed ....
  lambda_min=0.d0
  lambda_max=1.d0
  delta_lambda=(lambda_max-lambda_min)/dble(npoints)
  allocate(llambda(0:npoints))
  forall (ia=0:npoints) llambda(ia)=lambda_min+delta_lambda*dble(ia)
 
!  when we patch with 3 images back and forth this should be changed ....
!  for the moment is assumed the the first image neb (the image number 1)  has the reaction coordinate 0.d0 
!  and the last image (the image number =  nimage_neb) has the reaction coordinate 1.d0
   
!  the images are defined from 1 to nimage_neb. 1 being the starting_min and nimage_neb being the final_min
 delta_neb=(1.d0-0.d0)/(dble(nimage_neb-1))
 allocate (coord_neb(nimage_neb))
 forall (ia=1:nimage_neb) coord_neb(ia)=0.d0+(dble(ia-1))*delta_neb


 allocate (force_defect(0:npoints),free_energy_force_defect(0:npoints))


return

end subroutine init_neb_reaction

subroutine read_neb_from_ndm()
!   NDM part ...
    USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY : im,imm,im_glob,fnam,lenfnam,at,bg,normat, zl,zls2,nzl,zero
    USE var_pot, ONLY: ntyp
! MAB part  ...
    USE mab_in_ndm_module, ONLY : nimage_neb
    USE reaction_neb_module, ONLY: neb_images
    implicit none

! few local variables used to read files .... 
    integer  :: ielat(imm)
    integer  :: iwmax(imm)
    integer  :: ityp(imm)
    real(double)  :: xp(3,imm)
    real(double)  :: xpp(3,imm)
    real(double)  :: vp(3,imm)
    real(double)  :: ax(3,imm)
    real(double)  :: fp(3,imm)


    integer :: ic, ip,lucin,icintype,typmax,i,typmin
    character :: extension*9
    character :: fnamneb*80
    integer :: npath
    integer, dimension(:), allocatable :: natoms_type

    npath=nimage_neb 

    write(6,*) 'reading the follwing images ....'
       do ip=1, npath
          write(extension,'(i9.9)') ip
          fnamneb=fnam(1:lenfnam)//'.coutposition.'//extension
          write(6,'(2a)')'image = ',fnamneb
          lucin = 93
          open(unit=lucin, file=fnamneb, form='unformatted', status='old')
          read (lucin) icintype
          if (icintype>=2) then
             read (lucin) at
             call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
             do ic = 1, 3
                normat(ic) = 0
                normat(ic) = normat(ic)+sum(at(:,ic)**2)
                normat(ic) = sqrt(normat(ic))
                zl(ic) = normat(ic)
             end do
             zls2 = zl/2.0

          else
             read (lucin) zl                      !size of the box
             write (6, *) 'zl ', zl
             at(1,1)=zl(1)
             at(2,2)=zl(2)
             at(3,3)=zl(3)
             at(1,2)=zero
             at(1,3)=zero
             at(2,1)=zero
             at(2,3)=zero
             at(3,1)=zero
             at(3,2)=zero
             call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3))
             zls2 = zl/2.0

          endif

          read (lucin) im                         !number of atoms in the box
          im_glob=im

          if ( (im>imm).OR.(im.LE.0) ) then
              WRITE(6,'(2a)') 'File: ', Trim(fnamneb)
              WRITE(6, '(2(a,i0,1x))') 'im = ', im, ' - imm = ', imm
             stop
          endif


          read (lucin) ityp                       !types
          if (allocated(natoms_type)) deallocate(natoms_type)
          allocate(natoms_type(ntyp))
          natoms_type(:ntyp) = 0
          typmax = -2
          typmin = 100

          typmax = max(maxval(ityp(:im)),typmax)

          i = 1
          if (im>0) then
             typmin = min(minval(ityp(:im)),typmin)
             i = im+1
          endif

          if (typmax>ntyp.or.typmin<1) then
             write (6, *) 'wrong ityp(', i, ')= ', ityp(i)
             stop
          endif

          do i=1,im
             natoms_type(ityp(i))=natoms_type(ityp(i))+1
          enddo
          read (lucin) xp
          allocate(neb_images(ip)%natoms_type(ntyp))
          neb_images(ip)%xp_neb(1:3,1:imm)=xp(1:3,1:imm)
          neb_images(ip)%fp_neb(1:3,1:imm)=fp(1:3,1:imm)
          neb_images(ip)%iwmax(1:imm)=iwmax(1:imm)
          neb_images(ip)%ielat(1:imm)=ielat(1:imm)
          neb_images(ip)%ityp(1:imm)=ityp(1:imm)
          neb_images(ip)%ntyp=ntyp
          neb_images(ip)%natoms_type(1:ntyp)=natoms_type(1:ntyp)
!          ielat_n   (:,ip) = ielat   (:) 
!          iwmax_n   (:,ip) = iwmax   (:)
!          ityp_n    (:,ip) = ityp    (:)
!          xp_n (:,:,ip)    = xp (:,:)
!          xpp_n(:,:,ip)    = xpp(:,:)
!          vp_n (:,:,ip)    = vp (:,:)
!          ax_n (:,:,ip)    = ax (:,:)
!          fp_n (:,:,ip)    = fp (:,:)
          close(lucin)
       end do

return


end subroutine read_neb_from_ndm




subroutine interpolate_the_neb_images()
! NDM part ...
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY : im,imm,im_glob,fnam,lenfnam,at,bg,normat, zl,zls2,nzl,zero,angst
  USE var_pot, ONLY: ntyp
! MAB part ...
  USE mab_in_ndm_module, ONLY : nimage_neb,nimage_lambda,nhisto
  USE reaction_neb_module, ONLY: neb_images,lambda_images,atoms_on_spline, & 
                                llambda, coord_neb,& 
                                lambda_min,lambda_max,delta_lambda,delta_neb,npoints

  implicit none
  integer :: i,ia,ix,ival
  real(double) :: rbuffer(nimage_neb)

 do ia=1,im
  do ix=1,3
     forall (i=1:nimage_neb) rbuffer(i)= neb_images(i)%xp_neb(ix,ia)
     call cspline(nimage_neb, coord_neb(1:nimage_neb), rbuffer(1:nimage_neb), & 
                   atoms_on_spline(ix,ia)%b(1:nimage_neb), &
                   atoms_on_spline(ix,ia)%c(1:nimage_neb), &
                   atoms_on_spline(ix,ia)%d(1:nimage_neb))
  end do
 end do



 do ia=1,im
    do ix=1,3
      do i=0,npoints

        ival=int(llambda(i)/delta_neb)+1  
        ! pay attention to fact that the grid of neb images is from  1 to nimage_neb
        !                        and the grid of lambda if from 0 to npoints ...... 
        lambda_images(i)%xp_neb(ix,ia)=neb_images(ival)%xp_neb(ix,ia)   + &
                   atoms_on_spline(ix,ia)%b(ival)*(llambda(i)-coord_neb(ival))    +    & 
                   atoms_on_spline(ix,ia)%c(ival)*(llambda(i)-coord_neb(ival))**2 + & 
                   atoms_on_spline(ix,ia)%d(ival)*(llambda(i)-coord_neb(ival))**3  
        lambda_images(i)%dxp_neb(ix,ia)= atoms_on_spline(ix,ia)%b(ival)   +    & 
              2.d0*atoms_on_spline(ix,ia)%c(ival)*(llambda(i)-coord_neb(ival)) + & 
              3.d0*atoms_on_spline(ix,ia)%d(ival)*(llambda(i)-coord_neb(ival))**2  
      end do

    end do
  end do

!some debug ....
!  do i = 0, npoints
!    ival=int(lambda(i)/delta_neb)+1
!
!    func=neb_images(ival)%xp_neb(1,1)   + &
!         atoms_on_spline(1,1)%b(ival)*(lambda(i)-coord_neb(ival))    +    & 
!         atoms_on_spline(1,1)%c(ival)*(lambda(i)-coord_neb(ival))**2 + & 
!         atoms_on_spline(1,1)%d(ival)*(lambda(i)-coord_neb(ival))**3  
!   
!   write(775,*) lambda(i), func*angst
!
!  end do
!
! do i=1, nimage_neb
!   write(766,*) coord_neb(i),neb_images(i)%xp_neb(1,1)*angst
! end do

return



end subroutine interpolate_the_neb_images



subroutine test_energy_along_reaction()
! NDM part 
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use tab_imm_m
  use var_pot
 ! MAB part ...
  USE mab_in_ndm_module, ONLY : ene0
  USE reaction_neb_module, ONLY: neb_images,lambda_images,atoms_on_spline, &
                                 llambda, &
                                lambda_min,lambda_max,delta_lambda,npoints, &
                                force_defect, free_energy_force_defect

! use mab_in_ndm_module
  implicit none
  !-----------------------------------------------
  real(double), dimension(3,imm) :: fp0
  integer :: i, ia,ix,i_loop
  real(double) :: ene_ini
 
force_defect(0:npoints)=0.d0

do i=0,npoints
  imd = im
  nad(:ntyp) = na(:ntyp)
  vp=0.0
  fp=0.
  do ia=1,im
   do ix=1,3
      xp(ix,ia)=lambda_images(i)%xp_neb(ix,ia)
   end do
  end do

  call caltabt
  call caltabi
  call calfo
  if (i==0) ene_ini=potist
  ene0=potist
  do ia=1,im
   do ix=1,3
      force_defect(i)=force_defect(i)+lambda_images(i)%dxp_neb(ix,ia)*fp(ix,ia)
   end do 
  end do
!debug_defect   
write (779,*) llambda(i), (ene0-ene_ini)*erg2ev
!debug_defect   
write (780,*) llambda(i), force_defect(i)
 fp0=fp
end do

free_energy_force_defect(0)=0.d0
   do i_loop=1,npoints
    !
    free_energy_force_defect(i_loop)=free_energy_force_defect(i_loop-1)+0.5d0*delta_lambda*(force_defect(i_loop-1)+force_defect(i_loop))
!debug_defect     
write(781,*) llambda(i_loop), -free_energy_force_defect(i_loop)*erg2ev 
   enddo
!maybe a simple ideea for precond with the temperature:
!   forall(i_loop=0:nhisto) Free_temp(i_loop)=exp(-Free_energy(i_loop)/temperature)
!   renorm_f=temperature*log(sum(Free_temp)*delta_z)
!   renormalise par rapport a l'aire
!   forall(i_loop=-nhisto1:nhisto+nhisto1) Free_energy(i_loop)=Free_energy(i_loop)+renorm_f




end subroutine test_energy_along_reaction



