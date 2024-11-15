module plottpcel_mod
  use boxconfig,only:box_config
  USE cellconfig,only:cell_config, cell_config_arps
  USE T_kind_param_m, ONLY:  double
  USE arret_ndm_mod,only: arret_ndm
  use gen_com_m,only:iteration,unitP,timel,rang
  use newunit_mod,only:newunit
  use Tpara,only:para_space_config    
#ifdef PARA
  USE Tpara,only:myidsp,nprocs,mpi_comm_world,comm_space
#else
  USE Tpara,only:myidsp,nprocs,nprocspace
#endif

  USE cryst_to_cart_mod,only: cryst_to_cart

 implicit none
 integer::iplotcel ! triggers the detailled analysis of ltpcel (0or 2  = std; 1 or 2 =specific)

  logical::  lpsph,lppl
  integer:: slxyz(3)
  integer::iplotformat=0
  real(double)::Rplt, posplt(3)
  integer::iteplotcomp=0
  
  type,extends(cell_config)::slice_config
     integer,allocatable::indc(:,:)
     integer::ncs
   contains
     procedure, pass::build=>build_slice
     procedure, pass::merge=>merge_slice
  end type slice_config
  
  type,abstract::domain_config
     integer,allocatable::indc(:,:),nato(:)
     integer,allocatable::ncs(:)
     integer::ndom,maxncs
     real(double),allocatable::tempc(:),Pr(:),vol(:)
   contains
     procedure(build_i), deferred,pass::build
     procedure, pass::merge=>merge_domain
     procedure,pass:: plotd=>plotdomain
  end type domain_config

  abstract interface
     subroutine build_i(sphere,celcf,box)
       import domain_config
       import cell_config
       import box_config
       class(domain_config)::sphere
       class(cell_config)::celcf
       class(box_config)::box
     end subroutine build_i
  end interface

  type,extends(domain_config)::sphere_config
     real(double)::posplt(3),rplt
   contains
     procedure, pass::build=>build_sphere
  end type sphere_config
  

!!$  type::sphere_config
!!$
!!$     integer,allocatable::indc(:,:),nato(:)
!!$     integer,allocatable::ncs(:)
!!$     integer::ndom,maxncs
!!$     real(double),allocatable::tempc(:),Pr(:),vol(:)
!!$     real(double)::posplt(3),rplt
!!$   contains
!!$     procedure, pass::build=>build_sphere
!!$     procedure, pass::merge=>merge_domain
!!$     procedure,pass:: plotd=>plotdomain
!!$
!!$  end type sphere_config
  
  contains


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


    namelist /ltpc/lppl,posplt,Rplt,slxyz ,lpsph,iteplotcomp,iplotformat,iplotcel !est enlve de la namelist pour déacriver cette partien qui cree un bug (écrase cm pour une raison inconnue)
    lpsph=.false.
    lppl=.false.
    posplt(:)=0.5
    Rplt=6.0
    slxyz(:)=0 ! slx(1)=3 => average along 3 cells along x; slx(1)=0 (defalut = average along all X
    iplotcel=0 
    icall=icall+1
    if (.not.celcf%ltpcel) then
       write(6,*)'coding error in plottpcel call, %ltpcel.ne.true'
       call arret_ndm
    end if
!    if (icall==1) then 
       call newunit(iultp)
       open(unit=iultp,file='ltpcel.in')
       read(iultp,nml=ltpc)
       close(iultp)
!    end if
    
    if (present(itapp)) then
       itp=itapp
    else
       itp=iteration
    end if
    if (rang==0)write(6,*)'PLOTCEL',itp,iplotcel,iteplotcomp
    if (mod(iplotcel,2)==0) then
       if (iteplotcomp.gt.0) then
          if (mod(itp,iteplotcomp)==0)then
             if (myidsp==0) then
                call actualplot(celcf,boxcf,'TEMPC','PRESSC',itp)
             end if
          end if
       end if
    end if
    
    if (iplotcel.ge.1) then

       block
         type(cell_config)::cellcomp
         call celcf%constrcomp(cellcomp,boxcf,latomcp=.false.,psc=psc)
         
         
         if (myidsp==0) then
            if (lppl) then
               if (all(slxyz==0)) then
                  write(6,*)'inconsistent slxyz=0 and lppl'
                  call arret_ndm
               else
                  call plotslice(cellcomp,boxcf,slxyz,itp)
                  
               end if
            end if
            if (lpsph)  then
               call plotsph(cellcomp,boxcf,posplt,rplt,itp)
            end if
         end if
       end block
    end if
  end subroutine plottpcel


    subroutine build_sphere(sphere,celcf,box)
    class(sphere_config)::sphere
    class(cell_config)::celcf
    class(box_config)::box

!!  end subroutine build_sphere
  
    integer::iko,kos,i1,i2,i3,ic,ndom,idom
    real(double)::dist,posr(3),cv(1,3),distm,edge(3)
    if ((any(sphere%posplt(:).gt.1)).or.(any(sphere%posplt(:).lt.0))) then
       write(6,*)' check 0<=POSPLT <=1'
       call arret_ndm
    end if

    cv(1,:)=sphere%posplt
    call cryst_to_cart(1,cv,box%at,1)
    posr=cv(1,:)
    distm=0.    
    do i1=0,1
       do i2=0,1
          do i3=0,1
             do ic=1,3
                edge(ic)=i1*box%at(ic,1)+i2*box%at(ic,2)+i3*box%at(ic,3)
             end do
             cv(1,:)=posr(:)-edge(:)
             call cryst_to_cart (1, cv, box%bg, -1) !cart vers cryst cryst vers cart sur cv
             do ic=1,3
                if (box%ipbc(ic)==1) then
                   if ( (cv(1,ic).GT.0.5d0).OR.(cv(1,ic).LT.-0.5d0) )then
                      cv(1,ic) = cv(1,ic) - Dble(Nint(cv(1,ic)))
                   end if
                end if
             end do
             call cryst_to_cart (1, cv, box%at, 1) !cryst vers cart sur c
             dist=norm2(cv(1,:)) ! dist en cm
             distm=max(distm,dist)

          end do
       end do
    end do
    ndom=1+int(distm/sphere%rplt)
    allocate(sphere%ncs(ndom))
    allocate(sphere%nato(ndom))
    sphere%nato=0
    allocate(sphere%Pr(ndom))
    allocate(sphere%vol(ndom))
    allocate(sphere%tempc(ndom))
    sphere%ncs=0
    sphere%ndom=ndom
    sphere%vol=0.
    
    do iko=1,celcf%noxyz
       edge(:)=celcf%edge(iko,box)
       cv(1,:)=edge(:)-posr(:)
       call cryst_to_cart (1, cv, box%bg, -1) !cart vers cryst cryst vers cart sur cv
       do ic=1,3
          if (box%ipbc(ic)==1) then
             if ( (cv(1,ic).GT.0.5d0).OR.(cv(1,ic).LT.-0.5d0) )then
                cv(1,ic) = cv(1,ic) - Dble(Nint(cv(1,ic)))
             end if
          end if
       end do
       call cryst_to_cart (1, cv, box%at, 1) !cryst vers cart sur c
       dist=norm2(cv(1,:)) ! dist en cm
       idom=1+int(dist/sphere%rplt)
       sphere%ncs(idom)=sphere%ncs(idom)+1
    end do
    sphere%maxncs=maxval(sphere%ncs(:))
    allocate(sphere%indc(sphere%maxncs,sphere%ndom))
    sphere%ncs=0    
    do iko=1,celcf%noxyz
       edge(:)=celcf%edge(iko,box)
       cv(1,:)=edge(:)-posr(:)
       call cryst_to_cart (1, cv, box%bg, -1) !cart vers cryst cryst vers cart sur cv
       do ic=1,3
          if (box%ipbc(ic)==1) then
             if ( (cv(1,ic).GT.0.5d0).OR.(cv(1,ic).LT.-0.5d0) )then
                cv(1,ic) = cv(1,ic) - Dble(Nint(cv(1,ic)))
             end if
          end if
       end do
       call cryst_to_cart (1, cv, box%at, 1) !cryst vers cart sur c
       dist=norm2(cv(1,:)) ! dist en cm
       idom=1+int(dist/sphere%rplt)
       sphere%ncs(idom)=sphere%ncs(idom)+1
!       write(6,*)'idom',idom,ndom,sphere%ncs(idom),sphere%maxncs
       sphere%indc(sphere%ncs(idom),idom)=iko

       sphere%nato(idom)=sphere%nato(idom)+celcf%nato(iko)
       sphere%vol(idom)=sphere%vol(idom)+box%volu/celcf%noxyz
    end do
    
  end subroutine build_sphere
  
  subroutine merge_domain(domain,celcf,box)
    class(domain_config)::domain
    class(cell_config)::celcf
    class(box_config)::box

    integer::is,ik,isc
    real(double)::pc
    domain%tempc(:)=0
    domain%Pr(:)=0
    do is=1,domain%ndom
       do isc=1,domain%ncs(is)
          ik=domain%indc(isc,is)
          domain%tempc(is)=domain%tempc(is)+ (celcf%tempc(ik)*celcf%nato(ik))/domain%nato(is)
          Pc=((box%volu/celcf%noxyz)/domain%vol(is))*&
               &(celcf%sigc(1,1,ik)+celcf%sigc(2,2,ik)+celcf%sigc(3,3,ik))/3.
          domain%Pr(is)=domain%Pr(is)+ Pc
       end do
    end do
  end subroutine merge_domain

  subroutine plotdomain(domain,ctemp,cpress,itp)
    class(domain_config)::domain
    character(len=*)::ctemp,cpress
!    real(double)::,optional::r0
    integer::itp

    integer::i,unitlt,unitlp
    character :: extension*9,namef*80

    call newunit(unitlt)
    write(extension,'(i9.9)')itp
    namef=trim(ctemp)//trim(extension)
    open(unit=unitlt,file=namef,form='formatted')
    
    do i=1,domain%ndom
       write(unitlt,'(I12,G17.5,I8)')i,domain%tempc(i),domain%nato(i)
    end do
    close(unitlt)
    
    call newunit(unitlp)
    write(extension,'(i9.9)')itp
    namef=trim(cpress)//trim(extension)
    open(unitlp,file=namef,form='formatted')
    
    do i=1,domain%ndom
       write(unitlp,'(I12,G18.5,I8)'),i,domain%pr(i)*unitP,domain%nato(i)
    end do
    close(unitlt)
  end subroutine plotdomain


    
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

    ratiox(1)=celcf%nox/nx(1);    ratiox(2)=celcf%noy/nx(2);    ratiox(3)=celcf%noz/nx(3)
    natperc=celcf%noxyz*celcf%natperc/nx(1)*nx(2)*nx(3)

    call slice%cell_config%init(box,nx(1),nx(2),nx(3),natperc,ltpc=.true.,latomalloc=.false.)

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
    !    write(6,*)'RANG',slice%nato
    
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




  subroutine plotsph(celcf,box,posplt,rplt,itp)
    type(cell_config)::celcf
    type(box_config)::box
    integer::itp
    real(double)::posplt(3),rplt

    type(sphere_config)::sphere
    
    rplt=rplt*1d-8
    posplt=posplt
    sphere%posplt=posplt
    sphere%rplt=rplt
    call sphere%build(celcf,box)
    call sphere%merge(celcf,box)
    call sphere%plotd('Tsph','Psph',itp)

  end subroutine plotsph
  
  subroutine actualplot(celcf,boxcf,ctemp,cpress,itp)
    class(cell_config)::celcf
    class (box_config),intent(in)::boxcf
    character(len=*)::ctemp,cpress
    integer::itp
    
    integer::i,unitlt,koxyz(3),unitlp,kx,ky,kz,koo,il
    character :: extension*9,namef*80
    real(double)::pcell,opedg(3)

    
    select case (iplotformat)
    case(1)
       call newunit(unitlt)
       write(extension,'(i9.9)')itp
       namef=trim(ctemp)//trim(extension)//'.cube'
       open(unitlt,file=namef,form='formatted')
       write(unitlt,*)iteration,timel
       write(unitlt,*)'TEMP PER CELL'
       write(unitlt,*) '2 0.0 0.0 0.0'
       write(unitlt,'(I6,3G17.9)')celcf%nox,1d8*boxcf%at(1,1)/celcf%nox,1d8*boxcf%at(2,1)/celcf%nox,1d8*boxcf%at(3,1)/celcf%nox
       write(unitlt,'(I6,3G17.9)')celcf%noy,1d8*boxcf%at(1,2)/celcf%noy,1d8*boxcf%at(2,2)/celcf%noy,1d8*boxcf%at(3,2)/celcf%noy
       write(unitlt,'(I6,3G17.9)')celcf%noz,1d8*boxcf%at(1,3)/celcf%noz,1d8*boxcf%at(2,3)/celcf%noz,1d8*boxcf%at(3,3)/celcf%noz
       write(unitlt,'(A)') '1 0.0 0.0 0.0'
       opedg(:)=boxcf%at(:,1)+boxcf%at(:,2)+boxcf%at(:,3)
       write(unitlt,'(A,3G17.9)') '1 ', 1d8*opedg(1:3)

       il =0
       do kx=0,celcf%nox-1
          do ky=0,celcf%noy-1
             do kz=0,celcf%noz-1
                koo = 1+kx+celcf%nox*(ky+celcf%noy*kz)
                il=il+1
                write(unitlt,'(G15.7)',advance='no')celcf%tempc(koo)
                if (mod(il,6)==0 )write(unitlt,*)
             end do
          end do
       end do
       
       close(unitlt)


       call newunit(unitlp)
       write(extension,'(i9.9)')itp
       namef=trim(cpress)//trim(extension)//'.cube'
       open(unitlp,file=namef,form='formatted')
       write(unitlt,*)iteration,timel
       write(unitlt,*)'PRESS PER CELL'
       write(unitlt,*) '2 0.0 0.0 0.0'
       write(unitlt,'(I6,3G17.9)')celcf%nox,1d8*boxcf%at(1,1)/celcf%nox,1d8*boxcf%at(2,1)/celcf%nox,1d8*boxcf%at(3,1)/celcf%nox
       write(unitlt,'(I6,3G17.9)')celcf%noy,1d8*boxcf%at(1,2)/celcf%noy,1d8*boxcf%at(2,2)/celcf%noy,1d8*boxcf%at(3,2)/celcf%noy
       write(unitlt,'(I6,3G17.9)')celcf%noz,1d8*boxcf%at(1,3)/celcf%noz,1d8*boxcf%at(2,3)/celcf%noz,1d8*boxcf%at(3,3)/celcf%noz
       write(unitlt,*) '1 0.0 0.0 0.0'
       opedg(:)=boxcf%at(:,1)+boxcf%at(:,2)+boxcf%at(:,3)
       write(unitlt,'(A,3G17.9)') '1 ', opedg(1:3)

       il =0
       do kx=0,celcf%nox-1
          do ky=0,celcf%noy-1
             do kz=0,celcf%noz-1
                koo = 1+kx+celcf%nox*(ky+celcf%noy*kz)
                il=il+1
                pcell=0.33333333333333333*(celcf%sigc(1,1,koo)+celcf%sigc(2,2,koo)+celcf%sigc(3,3,koo))*unitP
                write(unitlp,'(G15.7)',advance='no')pcell
                if (mod(il,6)==0 )write(unitlp,*)
             end do
          end do
       end do

       
       close(unitlp)

    case(0)
       
       call newunit(unitlt)
       write(extension,'(i9.9)')itp
       namef=trim(ctemp)//trim(extension)
       open(unitlt,file=namef,form='formatted')
     
       do i=1,celcf%noxyz
          koxyz=celcF%koxyz(i)
#ifdef PARA
          select type (celcf)
          type is (cell_config)
             write(unitlt,'(I12,3I5,G21.5)'),i,koxyz(1:3),celcf%tempc(i)
          type is (slice_config)
             write(unitlt,'(I12,3I5,G21.5)'),i,koxyz(1:3),celcf%tempc(i)
          type is (cell_config_arps)
             write(unitlt,'(I12,3I5,G21.5)')i,koxyz(1:3),celcf%tempc(i)
          end select
#else
          select type (celcf)
          type is (cell_config)
             write(unitlt,'(I12,3I5,G21.5,I8)'),i,koxyz(1:3),celcf%tempc(i),celcf%nato(i)
          type is (slice_config)
             write(unitlt,'(I12,3I5,G21.5,I8)'),i,koxyz(1:3),celcf%tempc(i),celcf%nato(i)
          type is (cell_config_arps)
             write(unitlt,'(I12,3I5,G21.5,4I8)')i,koxyz(1:3),celcf%tempc(i),celcf%nato(i)&
                  &,celcf%nmov(0,i),celcf%nmov(1,i),celcf%nmov(2,i)
          end select

#endif
          
       end do
       close(unitlt)
       
       call newunit(unitlp)
       write(extension,'(i9.9)')itp
       namef=trim(cpress)//trim(extension)
       open(unitlp,file=namef,form='formatted')
       
       do i=1,celcf%noxyz
          koxyz=celcF%koxyz(i)
          pcell=0.33333333333333333*(celcf%sigc(1,1,i)+celcf%sigc(2,2,i)+celcf%sigc(3,3,i))*unitP
#ifdef PARA
          select type (celcf)
          type is (cell_config)
             write(unitlp,'(I12,3I5,G21.5,I8,9E18.5)'),i,koxyz(1:3),pcell,celcf%nato(i),celcf%sigc(:,:,i)*unitP
          type is (slice_config)
             write(unitlp,'(I12,3I5,G21.5,I8,9E18.5)'),i,koxyz(1:3),pcell,celcf%nato(i),celcf%sigc(:,:,i)*unitP
          type is (cell_config_arps)
             write(unitlp,'(I12,3I5,G21.5,4I8,9E18.5)')i,koxyz(1:3),pcell,celcf%nato(i),celcf%nmov(0,i),&
                  &celcf%nmov(1,i),celcf%nmov(2,i),celcf%sigc(:,:,i)*unitP
          end select
#else
          select type (celcf)
          type is (cell_config)
             write(unitlp,'(I12,3I5,G21.5,9E18.5)'),i,koxyz(1:3),pcell,celcf%sigc(:,:,i)*unitP
          type is (slice_config)
             write(unitlp,'(I12,3I5,G21.5,9E18.5)'),i,koxyz(1:3),pcell,celcf%sigc(:,:,i)*unitP
          type is (cell_config_arps)
             write(unitlp,'(I12,3I5,G21.5,9E18.5)')i,koxyz(1:3),pcell,celcf%sigc(:,:,i)*unitP
          end select
#endif
          
       end do
       close(unitlt)
    end select
       
  end subroutine actualplot

    
  subroutine plotslice (celcf,box,slxyz,itp)
    type(cell_config)::celcf
    type(box_config)::box
    integer::slxyz(3),itp
    type(slice_config)::slice
    call slice%build(celcf,box,slxyz)
    call slice%merge(celcf)
!    call slice%print(unit=100,mess='SLICE')
    call actualplot(slice,box,'TSlice','PSlice',itp)
    
  end subroutine plotslice
  

         
       
end module plottpcel_mod
