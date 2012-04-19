subroutine rdf
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m

  !      USE coordo_m
  !******************************************************************
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------

  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, iti, itj, i1, i2, icell, kx, ky, kz, koo, ko1, j, &
       ic, k, m,m1,lucoord,n
  real(double) :: rij, c1, c2, c3, x1, x2, x3, rmax, incre
  real(double) :: rspace2,invincre
  real(double) :: aaa,bbb,ccc,ddd
  real(double) :: digrt(kmax),digrt2(kmax)

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
  coord(:,:,:)=0.D0
  lucoord = 11
!  rmax=minval(celsize)
  rmax=rcrdf*1.0d-8
  !       write(6,*)'rmax ',rmax
  incre = rmax/kmax
  invincre = 1/incre
  if(rang==0) then
     write(6,*)
     write(6,*)'------------------------------------------'
     write(6,*)'--------Calcul des Fonctions de correlation---------'
     write(6,*)'nrdf ',nrdf

     write (6, '(A,I5,A,D10.3)') '*  ITERATION  = ', it, '  time = ', timel
  end if

  if(.not.linstantrdf) then

     do i1=1,ntyp
        digrt(:)=0. ; digrt2(:)=0.
        if(nad(i1)==0) cycle
        do i2=1,ntyp
           if(nad(i2)==0) cycle

           paire1=ty(i1)
           paire2=ty(i2)
           fpaire1=paire1
           fpaire2=paire2
           lenfpaire1=index(fpaire1,' ')-1
           lenfpaire2=index(fpaire2,' ')-1
           fpairecoord = fpaire1(1:lenfpaire1)// &
                '_'//fpaire2(1:lenfpaire2)//'.coord'
           if(rang==0)               open(lucoord, file = fpairecoord, status = 'unknown')

           do m1=1,kmax
              !                  write(6,*)i1,i2,n,digr(i1,i2,n)
              m = m1+1

              do k = 1, m1
                 coord(i1,i2,m1)=coord(i1,i2,m1)+digr(i1,i2,k)
                 digrt2(m1)=digrt2(m1)+digr(i1,i2,k)/(nrdf*nad(i1)*im)
              end do

              rspace2 = (m1*incre)**2

              digrt(m1)=digrt(m1)+digr(i1,i2,m1)*volu/(4*pi*incre*rspace2* &
                   nad(i1)*im*nrdf)


              bbb = digr(i1,i2,m)*volu/(4*pi*incre*rspace2* &
                   nad(i1)*nad(i2)*nrdf)
              aaa=m1*incre 
              ccc=coord(i1,i2,m1)/(nrdf*nad(i1))
              if(rang==0) write (lucoord,*) aaa,bbb,ccc
           end do
           if(rang==0)         close(lucoord)

        end do

        fpaire1=ty(i1); lenfpaire1=index(fpaire1,' ')-1
        fpairecoord = fpaire1(1:lenfpaire1)// &
             '.coord'
        if(rang==0) then
           open(lucoord, file = fpairecoord, status = 'unknown')
           do m1=1,kmax
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

     do m1=1,kmax
        m = m1+1
        rspace2 = (m1*incre)**2
        ddd=gdertot(m)*volu/(4*pi*incre*rspace2)
        if (rang==0) write(35,*) m1*incre,ddd/nrdf
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
           stop
        endif
        rewind 33

        read (33, 1000) charsauvrfdc

        lusauvrdf=34

        write(6,*) ' sauvegarde RDF partielle it=',it
        write(6,*)
     end if
     do i1=1,ntyp
        if(nad(i1)==0) cycle
        digrt(:)=0. ; digrt2(:)=0.               
        do i2=1,ntyp
           if(nad(i2)==0) cycle

           paire1=ty(i1)
           paire2=ty(i2)
           fpaire1=paire1
           fpaire2=paire2
           lenfpaire1=index(fpaire1,' ')-1
           lenfpaire2=index(fpaire2,' ')-1
           fpairecoord = fpaire1(1:lenfpaire1)// &
                '_'//fpaire2(1:lenfpaire2)//'.'//charsauvrfdc
           if(rang==0)   open(lusauvrdf, file = fpairecoord, status = 'unknown')

           do m1=1,kmax
              m = m1+1
              do k = 1, m1
                 coord(i1,i2,m1)=coord(i1,i2,m1)+digr(i1,i2,k)
                 digrt2(m1)=digrt2(m1)+digr(i1,i2,k)/(nrdf*nad(i1)*im)
              end do
              rspace2 = (m1*incre)**2
              digrt(m1)=digrt(m1)+digr(i1,i2,m1)*volu/(4*pi*incre*rspace2* &
                   nad(i1)*im*nrdf)
              bbb = digr(i1,i2,m)*volu/(4*pi*incre*rspace2* &
                   nad(i1)*nad(i2))
              aaa=m1*incre 
                   ccc=coord(i1,i2,m1)/(nad(i1))
              if(rang==0)            write (lusauvrdf,*) aaa,bbb,ccc
           end do
           if(rang==0) close(lusauvrdf)
        end do



        fpaire1=ty(i1); lenfpaire1=index(fpaire1,' ')-1
        fpairecoord = fpaire1(1:lenfpaire1)// &
             '.coord.'//charsauvrfdc
        if(rang==0) then
           open(lucoord, file = fpairecoord, status = 'unknown')
           do m1=1,kmax
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
     do m1=1,kmax
        m = m1+1
        rspace2 = (m1*incre)**2
        ddd=gdertot(m)*volu/(4*pi*incre*rspace2)
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
  digr(:,:,:)=0.D0
  gdertot(:)=0.D0
  return
end subroutine rdf
