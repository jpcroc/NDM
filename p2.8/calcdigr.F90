module calcdigr_mod
  USE arret_ndm_mod,only:arret_ndm
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY:lperiod,rang,rcrdf,it,pi,timel,lspacendm
  use atomconfig,only: atom_config
  use boxconfig,only:box_config
  USE cellconfig,only:cell_config, caltabtC
  USE var_pot, ONLY:ntyp,ty !nkmax,ntyp,digr,gdertot
  use vect_dist_mod,only:vect_dist

#ifdef PARA
    USE mpi
    USE Tpara,only:COMM_space,nprocspace

#endif
  implicit none
  
  type rdf_typ
     real(double)::rcrdf
     real(double),allocatable::gdertot(:),digr(:,:,:),coord(:,:,:)
     integer::nkmax
     logical::linstantrdf
     integer::nrdf
     character::nom(2)
     real(double)::volu
     integer::im
     integer,allocatable::nad(:)
     
  end type rdf_typ
    type(rdf_typ)::rdf0
contains

  subroutine initrdf(rdfc,nkmax,rcrdf,linstantrdf,nom)
    type(rdf_typ),intent(out)::rdfc
    integer,intent(in)::nkmax
    real(double),intent(in)::rcrdf
    logical,intent(in)::linstantrdf
    character::nom(2)
    
    allocate(rdfc%digr(ntyp,ntyp,nkmax))
    allocate(rdfc%coord(ntyp,ntyp,nkmax))
    allocate(rdfc%nad(ntyp))
    rdfc%digr=0
    rdfc%coord=0
    allocate(rdfc%gdertot(nkmax))
    rdfc%gdertot=0
    rdfc%rcrdf=rcrdf
    rdfc%nkmax=nkmax
    rdfc%nrdf=0
    rdfc%linstantrdf=linstantrdf
    rdfc%volu=0
    rdfc%im=0
    rdfc%nom=nom
  end subroutine initrdf

  subroutine calcdigr(atrdf,celrdf,boxrdf,rdfc)
    !-----------------------------------------------
    !   M o d u l e s
    !-----------------------------------------------
    USE T_kind_param_m, ONLY:  double




    implicit none
    class(atom_config),intent(in)::atrdf
    type(box_config),intent(in)::boxrdf
    type(cell_config),intent(in):: celrdf
    type(rdf_typ)::rdfc

    
!!$    integer,intent(in)::im
!!$    real(double),intent(in),allocatable::xp(:,:)
!!$    integer,allocatable,intent(in)::ityp(:),ielat(:)

    integer :: i, iti, itj, i1, i2, icell, kx, ky, kz, koo, ko1, j, &
         ic, k, m,m1,n,iti1,iti2
    real(double) :: rij,  x1, x2, x3, rmax2,rmax, incre
    real(double) :: rspace2,invincre
    real(double) :: aaa, bbb, ccc,cv(1,3)
    real(double),allocatable::digrtemp(:,:,:)
    logical::linter
    real(double)::cx(3)
    
    allocate(digrtemp(ntyp,ntyp,rdfc%nkmax))
    digrtemp=0
    rmax=rdfc%rcrdf

    if (rdfc%im==0) then
       rdfc%im=atrdf%im
    else
       if (rdfc%im.ne.atrdf%im) then
          write(6,*)'IM atrdf incosistent with rdfc STOP'
          call arret_ndm
       end if
    end if
    if (rdfc%volu==0) then
       rdfc%volu=boxrdf%volu
    else
       if (rdfc%volu.ne.boxrdf%volu) then
          write(6,*)'VOLU atrdf inconistent with rdfc STOP'
          call arret_ndm
       end if
    end if
    
    if (rmax.gt.minval(celrdf%celsize)) then
       write(6,*)'diminuer nox, noy, noz'
       call arret_ndm
    end if
    !      write(6,*)'rmax ',rmax
    rdfc%nrdf=rdfc%nrdf+1
    rmax2 = rmax**2
    incre = rmax/rdfc%nkmax
    invincre = 1/incre
    if (rang==0) write(6,*) 'nkmax incre',rdfc%nkmax,incre
    do i = 1, atrdf%im
       koo = atrdf%ielat(i)
       do i1 = 0, 26
          ko1 = celrdf%ncel(koo,i1)
          do i2 = 1, celrdf%nato(ko1)
             j = celrdf%atincel(i2,ko1)
             if(j==i) cycle 

             call vect_dist(atrdf,celrdf,boxrdf,i,j,rum=rmax,dist=rij,linter=linter,lperiod=lperiod)
             if (.not.linter) cycle

             k= int(rij*invincre)

             m=k+1
             digrtemp(atrdf%ityp(i),atrdf%ityp(j),m)= digrtemp(atrdf%ityp(i),atrdf%ityp(j),m)+1
             !                 write(6,*)' digr ', digr(ityp(i),ityp(j),m)
          end do
       end do
   end do
    !------------------------------------------------
    ! Calcul de la RDF totale
    !------------------------------------------------
   do iti1=1,ntyp
       rdfc%nad(iti1)=count(atrdf%ityp==iti1)
    end do
#ifdef PARA
    if ((nprocspace.gt.1).and.(lspacendm.eqv..true.)) then
       call comm_space%sum(digrtemp)
       call comm_space%sum(rdfc%nad)
    end if
#endif

   rdfc%digr(:,:,:) = rdfc%digr(:,:,:)+digrtemp(:,:,:)
    do iti1=1,ntyp
       if(rdfc%nad(iti1)==0) cycle
       do iti2=1,ntyp
          if(rdfc%nad(iti2)==0) cycle
          do m=1,rdfc%nkmax
             rdfc%gdertot(m)=rdfc%gdertot(m)+rdfc%digr(iti1,iti2,m)/(rdfc%nad(iti1)*rdfc%nad(iti2))
          enddo

       enddo
    enddo

    !  if (rang==0) write(6,*) 'PARA-T sortie calcdigr'

    return
  end subroutine calcdigr


  
  subroutine rdfT(rdfc)
    implicit none
    type(rdf_typ)::rdfc

    integer :: i, iti, itj, i1, i2, icell, kx, ky, kz, koo, ko1, j, &
         ic, k, m,m1,lucoord,n
    real(double) :: rij, c1, c2, c3, x1, x2, x3, rmax, incre
    real(double) :: rspace2,invincre
    real(double) :: aaa,bbb,ccc,ddd
    real(double),allocatable :: digrt(:),digrt2(:)


    character :: paire1*20,paire2*20,paire3*20
    character :: fpaire1*80, fpaire2*80,fpaire3*80
    character :: fpairecoord*80
    integer :: lenfpaire1,lenfpaire2,lenfpaire3
    integer :: lusauvrdf
    character :: charsauvrfdc*8  

    !-----------------------------------------------
    !
    ! local variables
    !
    allocate(digrt(rdfc%nkmax),digrt2(rdfc%nkmax))
    rdfc%coord(:,:,:)=0.D0
    lucoord = 11
    !  rmax=minval(celsize)
    rmax=rcrdf*1.0d-8
    !       write(6,*)'rmax ',rmax
    incre = rmax/rdfc%nkmax
    invincre = 1/incre
    if(rang==0) then
       write(6,*)
       write(6,*)'------------------------------------------'
       write(6,*)'--------Calcul des Fonctions de correlation---------'
       write(6,*)'nrdf ',rdfc%nrdf

       write (6, '(A,I5,A,D10.3)') '*  ITERATION  = ', it, '  time = ', timel
    end if

    if(.not.rdfc%linstantrdf) then

       do i1=1,ntyp
          digrt(:)=0. ; digrt2(:)=0.
          if(rdfc%nad(i1)==0) cycle
          do i2=1,ntyp
             if(rdfc%nad(i2)==0) cycle

             paire1=ty(i1)
             paire2=ty(i2)
             fpaire1=paire1
             fpaire2=paire2
             lenfpaire1=index(fpaire1,' ')-1
             lenfpaire2=index(fpaire2,' ')-1
             fpairecoord = fpaire1(1:lenfpaire1)// &
                  '_'//fpaire2(1:lenfpaire2)//'.coord'
             if(rang==0)               open(lucoord, file = fpairecoord, status = 'unknown')

             do m1=1,rdfc%nkmax
                !                  write(6,*)i1,i2,n,digr(i1,i2,n)
                m = m1+1

                do k = 1, m1
                   rdfc%coord(i1,i2,m1)=rdfc%coord(i1,i2,m1)+rdfc%digr(i1,i2,k)
                   digrt2(m1)=digrt2(m1)+rdfc%digr(i1,i2,k)/(rdfc%nrdf*rdfc%nad(i1)*rdfc%im)
                end do

                rspace2 = (m1*incre)**2

                digrt(m1)=digrt(m1)+rdfc%digr(i1,i2,m1)*rdfc%volu/(4*pi*incre*rspace2* &
                     rdfc%nad(i1)*rdfc%im*rdfc%nrdf)


                bbb = rdfc%digr(i1,i2,m)*rdfc%volu/(4*pi*incre*rspace2* &
                     rdfc%nad(i1)*rdfc%nad(i2)*rdfc%nrdf)
                aaa=m1*incre 
                ccc=rdfc%coord(i1,i2,m1)/(rdfc%nrdf*rdfc%nad(i1))
                if(rang==0) write (lucoord,*) aaa,bbb,ccc
             end do
             if(rang==0)         close(lucoord)

          end do

          fpaire1=ty(i1); lenfpaire1=index(fpaire1,' ')-1
          fpairecoord = fpaire1(1:lenfpaire1)// &
               '.coord'
          if(rang==0) then
             open(lucoord, file = fpairecoord, status = 'unknown')
             do m1=1,rdfc%nkmax
                aaa=m1*incre
                write(lucoord,*)aaa,digrt(m1),digrt2(m1)
             end do
             close(lucoord)
          end if
       end do


       !------------------------------------------------
       ! RDF totale moyenne
       !------------------------------------------------
       !         if(lrdftot) then
       if(rang==0)         open(35,file='rdftot.moy',status='unknown')

       do m1=1,rdfc%nkmax
          m = m1+1
          rspace2 = (m1*incre)**2
          ddd=rdfc%gdertot(m)*rdfc%volu/(4*pi*incre*rspace2)
          if (rang==0) write(35,*) m1*incre,ddd/rdfc%nrdf
       enddo
       if(rang==0)         close(35)
       !         endif

    else                     !linstantrdf=.TRUE.
       if(rang==0)then
          open(unit=33,file='tamprfdc',form='formatted',status='unknown')

          if (it<=9) write (33, '(I1)') it
          if (it<=99.and.it>9) write (33, 200) it
          if (it<=999.and.it>99) write (33, 300) it
          if (it<=9999.and.it>999) write (33, 400) it
          if (it<=99999.and.it>9999) write (33, 500) it
          if (it<=999999.and.it>99999) write(33, 600) it
          if (it<=9999999.and.it>999999) write(33, 700) it
          if (it<=99999999.and.it>9999999) write(33, 800) it
          if (it<=999999999.and.it>99999999) write(33, 900) it
          if  (it>999999999) then
             write (6, *) 'probleme de format dans calccoordo.f90'
             call arret_ndm
          endif
          rewind 33

          read (33, 1000) charsauvrfdc

          lusauvrdf=34

          write(6,*) ' sauvegarde RDF partielle it=',it
          write(6,*)
       end if
       do i1=1,ntyp
          if(rdfc%nad(i1)==0) cycle
          digrt(:)=0. ; digrt2(:)=0.               
          do i2=1,ntyp
             if(rdfc%nad(i2)==0) cycle

             paire1=ty(i1)
             paire2=ty(i2)
             fpaire1=paire1
             fpaire2=paire2
             lenfpaire1=index(fpaire1,' ')-1
             lenfpaire2=index(fpaire2,' ')-1
             fpairecoord = fpaire1(1:lenfpaire1)// &
                  '_'//fpaire2(1:lenfpaire2)//'.'//charsauvrfdc
             if(rang==0)   open(lusauvrdf, file = fpairecoord, status = 'unknown')

             do m1=1,rdfc%nkmax-1
                m = m1+1
                do k = 1, m1
                   rdfc%coord(i1,i2,m1)=rdfc%coord(i1,i2,m1)+rdfc%digr(i1,i2,k)
                   digrt2(m1)=digrt2(m1)+rdfc%digr(i1,i2,k)/(rdfc%nrdf*rdfc%nad(i1)*rdfc%im)
                end do
                rspace2 = (m1*incre)**2
                digrt(m1)=digrt(m1)+rdfc%digr(i1,i2,m1)*rdfc%volu/(4*pi*incre*rspace2* &
                     rdfc%nad(i1)*rdfc%im*rdfc%nrdf)
                bbb = rdfc%digr(i1,i2,m)*rdfc%volu/(4*pi*incre*rspace2* &
                     rdfc%nad(i1)*rdfc%nad(i2))
                aaa=m1*incre 
                ccc=rdfc%coord(i1,i2,m1)/(rdfc%nad(i1))
                if(rang==0)            write (lusauvrdf,*) aaa,bbb,ccc
             end do
             if(rang==0) close(lusauvrdf)
          end do



          fpaire1=ty(i1); lenfpaire1=index(fpaire1,' ')-1
          fpairecoord = fpaire1(1:lenfpaire1)// &
               '.coord.'//charsauvrfdc
          if(rang==0) then
             open(lucoord, file = fpairecoord, status = 'unknown')
             do m1=1,rdfc%nkmax
                aaa=m1*incre
                write(lucoord,*)aaa,digrt(m1),digrt2(m1)
             end do
             close(lucoord)
          end if

       end do
       !------------------------------------------------
       ! RDF totale instantanee
       !------------------------------------------------
       if(rang==0)            write(6,*) ' sauvegarde RDF totale it=',it
       if(rang==0)            write(6,*)
       if(rang==0)            open(35,file='rdftot.'//charsauvrfdc,status='unknown')
       do m1=1,rdfc%nkmax
          m = m1
          rspace2 = (m1*incre)**2
          ddd=rdfc%gdertot(m)*rdfc%volu/(4*pi*incre*rspace2)
          if(rang==0)             write(35,*) m1*incre,ddd
       enddo
       if(rang==0)           close(35)
    endif
200 format(i2)
300 format(i3)
400 format(i4)
500 format(i5)
600 format(i6)
700 format(i7)
800 format(i8)
900 format(i9)
1000 format(a6)
    if(rang==0)        write(6,*)'--------------------------------------'
    if(rang==0)        write(6,*) '--------------------------------------'
!!!  if(rang==0)        close(lusauvrdf)
    if(rang==0)        close(33,status='DELETE')




    !      coorpart(:,:)=0.D0
    rdfc%digr(:,:,:)=0.D0
    rdfc%gdertot(:)=0.D0
    return
  end subroutine rdfT

end module calcdigr_mod
