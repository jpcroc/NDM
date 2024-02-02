module plottpcel_mod
  use boxconfig,only:box_config
  USE cellconfig,only:cell_config, cell_config_arps
  USE T_kind_param_m, ONLY:  double
  USE arret_ndm_mod,only: arret_ndm
  use gen_com_m,only:iteration,unitP
  use newunit_mod,only:newunit
  use Tpara,only:para_space_config    
#ifdef PARA
  USE Tpara,only:myidsp,nprocs,mpi_comm_world,comm_space
#else
  USE Tpara,only:myidsp,nprocs,nprocspace
#endif

   
     
!  USE atomconfig,only:atom_config,atom_config_d,atom_config_e


  implicit none
  

  
  integer::iplotcel

  logical::  lpsph,lppl
  integer:: slxyz(3)
  real(double)::Rplt, posplt(3)
  
  
  type,extends(cell_config)::slice_config
     integer,allocatable::indc(:,:)
     integer,allocatable::ncs
   contains
     procedure, pass::build=>build_slice
     procedure, pass::merge=>merge_slice
  end type slice_config
  
  contains

  subroutine build_slice(slice,celcf,box,slxyz)
    class(slice_config)::slice
    class(cell_config)::celcf
    type(box_config)::box
    integer,intent(in)::slxyz(3)

    integer::nx(3),natperc,ixc(3),ixs(3),ratiox(3),iko,kos
    integer,allocatable::ncs(:)
    
    if (slxyz(1).ne.0) then
       if (mod(celcf%nox,slxyz(1)).ne.0) then
          write(6,*)'choose a slxyz(1) which divides nox',celcf%nox
          call arret_ndm
       end if
       nx(1)=celcf%nox/slxyz(1)
    else
       nx(1)=1
    end if

    if (slxyz(2).ne.0) then
       if (mod(celcf%noy,slxyz(2)).ne.0) then
          write(6,*)'choose a slxyz(2) which divides noy',celcf%noy
          call arret_ndm
       end if
       nx(2)=celcf%noy/slxyz(2)
    else
       nx(2)=1
    end if
    
    if (slxyz(3).ne.0) then
       if (mod(celcf%noz,slxyz(3)).ne.0) then
          write(6,*)'choose a slxyz(3) which divides noz',celcf%noz
          call arret_ndm
       end if
       nx(3)=celcf%noz/slxyz(3)
    else
       nx(3)=1
    end if
    ratiox(1)=celcf%nox/nx(1);    ratiox(2)=celcf%nox/nx(2);    ratiox(3)=celcf%nox/nx(3)
    natperc=celcf%noxyz*celcf%natperc/nx(1)*nx(2)*nx(3)
!    write(6,*)'POINT1'
    call slice%cell_config%init(box,nx(1),nx(2),nx(3),natperc,ltpc=.true.,latomalloc=.false.)
!    write(6,*)'POINT2'
    slice%ncs=celcf%noxyz/slice%noxyz
    allocate(slice%indc(slice%ncs,slice%noxyz))
    allocate(ncs(slice%noxyz))
    ncs(:)=0
    slice%nato(:)=0
    do iko=1,celcf%noxyz
       ixc=celcf%koxyz(iko)
       ixs(:)=int(ixc(:)/float(ratiox(:)))
       kos=1+ixs(1)+slice%nox*(ixs(2)+slice%noy*ixs(3))
       ncs(kos)=ncs(kos)+1
       slice%indc(ncs(kos),kos)=iko
          slice%nato(kos)=slice%nato(kos)+celcf%nato(iko)
       
       end do
    write(6,*)'RANG',slice%nato
    
    do iko=1,slice%noxyz
       if (slice%ncs.ne.ncs(iko)) then
          write(6,*)'erreur constr slice'
       end if
    end do


    
  end subroutine build_slice
  
  subroutine merge_slice(slice,celcf)
    class(slice_config)::slice
    class(cell_config)::celcf

    integer::is,ik,isc
    slice%tempc(:)=0
    slice%sigc(:,:,:)=0
    do is=1,slice%noxyz
       do isc=1,slice%ncs
          ik=slice%indc(isc,is)
          slice%tempc(is)=slice%tempc(is)+ (celcf%tempc(ik)*celcf%nato(ik))/slice%nato(is)
          slice%sigc(:,:,is)=slice%sigc(:,:,is)+ (celcf%sigc(:,:,ik)*float(slice%noxyz)/celcf%noxyz)
       end do
    end do
  end subroutine merge_slice



  subroutine plottpcel(celcf,boxcf,itapp,psc)
    type(para_space_config)::psc    
    integer,optional::itapp
    class(cell_config)::celcf
    class(box_config)::boxcf
    integer::itp,unitlt,i,koxyz(3),unitlp,iultp
    character :: extension*9
    character*80::namef
    real(double)::Pcell
    integer,save::icall=0
    namelist /ltpc/lpsph,lppl,posplt,Rplt,slxyz
    lpsph=.false.
    lppl=.false.
    posplt(:)=0.5
    Rplt=6.0
    slxyz(:)=0 ! slx(1)=3 => average along 3 cells along x; slx(1)=0 (defalut = average along all X 
    icall=icall+1
    if (.not.celcf%ltpcel) then
       write(6,*)'coding error in plottpcel call, %ltpcel.ne.true'
       call arret_ndm
    end if
    
    if (present(itapp)) then
       itp=itapp
    else
       itp=iteration
    end if

    if (mod(iplotcel,2)==0) then
       if (myidsp==0) then
          call actualplot(celcf,'TEMPC','PRESSC',itp)
       end if
    end if
    
    if (iplotcel.ge.1) then

       block
         type(cell_config)::cellcomp
         call celcf%constrcomp(cellcomp,boxcf,latomcp=.false.,psc=psc)
         
         
         if (myidsp==0) then
            call newunit(iultp)
            open(unit=iultp,name='ltpcel.in')
            read(iultp,nml=ltpc)
            close(iultp)
            if (lppl) then
               if (all(slxyz==0)) then
                  write(6,*)'inconsistent slxyz=0 and lppl'
                  call arret_ndm
               else
                  call plotslice(cellcomp,boxcf,slxyz,itp)
                  
               end if
            end if
            if (lpsph)  then
               
            end if
         end if
       end block
    end if
  end subroutine plottpcel

  
  subroutine actualplot(celcf,ctemp,cpress,itp)
    class(cell_config)::celcf
    character(len=*)::ctemp,cpress
    integer::i,unitlt,itp,koxyz(3),unitlp
    character :: extension*9,namef*80
    real(double)::pcell
    
    call newunit(unitlt)
    write(extension,'(i9.9)')itp
    namef=trim(ctemp)//trim(extension)
    open(unitlt,file=namef,form='formatted')
    
    do i=1,celcf%noxyz
       koxyz=celcF%koxyz(i)
       select type (celcf)
       type is (cell_config)
          write(unitlt,'(I12,3I5,G15.5,I5)'),i,koxyz(1:3),celcf%tempc(i),celcf%nato(i)
       type is (slice_config)
          write(unitlt,'(I12,3I5,G15.5,I5)'),i,koxyz(1:3),celcf%tempc(i),celcf%nato(i)
       type is (cell_config_arps)
          write(unitlt,'(I12,3I5,G15.5,4I5)')i,koxyz(1:3),celcf%tempc(i),celcf%nato(i),celcf%nmov(0,i),celcf%nmov(1,i),celcf%nmov(2,i)
       end select
    end do
    close(unitlt)
    
    call newunit(unitlp)
    write(extension,'(i9.9)')itp
    namef=trim(cpress)//trim(extension)
    open(unitlp,file=namef,form='formatted')
    
    do i=1,celcf%noxyz
       koxyz=celcF%koxyz(i)
       pcell=(celcf%sigc(1,1,i)+celcf%sigc(2,2,i)+celcf%sigc(3,3,i))*unitP
       select type (celcf)
       type is (cell_config)
          write(unitlp,'(I12,3I5,G15.5,I5,9E15.5)'),i,koxyz(1:3),pcell,celcf%nato(i),celcf%sigc(:,:,i)*unitP
       type is (slice_config)
          write(unitlp,'(I12,3I5,G15.5,I5,9E15.5)'),i,koxyz(1:3),pcell,celcf%nato(i),celcf%sigc(:,:,i)*unitP
       type is (cell_config_arps)
          write(unitlp,'(I12,3I5,G15.5,4I5,9E15.5)')i,koxyz(1:3),pcell,celcf%nato(i),celcf%nmov(0,i),&
               &celcf%nmov(1,i),celcf%nmov(2,i),celcf%sigc(:,:,i)*unitP
       end select
    end do
    close(unitlt)
  end subroutine actualplot

    
  subroutine plotslice (celcf,box,slxyz,itp)
    type(cell_config)::celcf
    type(box_config)::box
    integer::slxyz(3),itp
    type(slice_config)::slice
   
    call slice%build(celcf,box,slxyz)
    call slice%merge(celcf)
!    call slice%print(unit=100,mess='SLICE')
    call actualplot(slice,'TEMPS','PRESSS',itp)
    
  end subroutine plotslice
  

         
       
end module plottpcel_mod
