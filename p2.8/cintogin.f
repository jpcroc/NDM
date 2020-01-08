      implicit none
      integer imM,ntyp
      parameter(ntyp=8)
      integer im,icintype,iti,icg,iv,itii,natav
      real*8 zl(3),nzl(3),xii(3),qtot,oldtstep,xt1,xt2,dilat
      real*8, pointer:: xpp(:,:),vp(:,:),ax(:,:),xp(:,:)
      integer, pointer ::ityp(:),num_at_glob(:)
      integer na(ntyp),lenfnam,igen,iout
      integer icintypemod
      integer lucin,lugin
      integer i,ic
      integer la,lb,lc
      integer nmod
      character *80 a1,fnam,fnamcin,fnamcout,fnamgin,
     &     fnamgout
      character*1 rep
      real*8 at(3,3),bg(3,3),normat(3)
      logical ltriclin


                                !      write(6,*)'CE PROGRAMME DOIT ETRE COMPILE AVEC '
                                !      write(6,*)'et avec  recips.o cryst_to_car. !!!!'
      write(6,*)'CINTOGIN =1 X to Z=2'
      read(5,*)nmod
      write(6,*)'nom ?'
      read(5,*)a1
      write(6,*)'imm ?'
      read(5,*) imm
      allocate(xp(3,imm))
      allocate(ax(3,imm))
      allocate(xpp(3,imm))
      allocate(vp(3,imm))
      allocate (ityp(imm))
      allocate (num_at_glob(imm))
                                !      write(6,*)'type de fichiers d_entree 0=>.gin 1=>.cin'
c     read(5,*)igen
      igen=1
      fnam=a1
      lenfnam=index(fnam,' ')-1
      lucin=8            
      fnamcin =fnam(1:lenfnam)//'.cout'
      fnamgin =fnam(1:lenfnam)//'.gin'
      fnamcout =fnam(1:lenfnam)//'.newcin'
      fnamgout =fnam(1:lenfnam)//'.newgin'

c     fichiers .cin
      open(unit=9,file=fnamcin,form='unformatted')
c     open(unit=9,file=fnamcout,form='unformatted')
      lucin=9
      read(lucin)icintype
      write(6,*)'icintype',icintype
      if ((icintype.gt.3).or.(icintype.lt.0)) then
         write(6,*)'wrong icintype'
         stop
      endif
      icintypemod=mod(icintype,2)

      if (icintype.ge.2) then
         ltriclin=.true.
         read(lucin) at
      write(6,*)'at',at
         call recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
         do ic=1,3
            normat(ic)=0
            do i=1,3
               normat(ic)=normat(ic)+at(i,ic)**2
            enddo
            normat(ic)=dsqrt(normat(ic))
            zl(ic)=normat(ic)
         enddo   
      else
         ltriclin=.false.
         read(lucin)zl          !size of the box
         write(6,*)'zl ',zl


c     do i=1,3
c     zls2(i)=zl(i)/2.0
c     enddo

      endif                     !icintype=2

      read(lucin)im             !number of atoms in the box
      write(6,*)'im =',im
      if (im.gt.imM) then
         write (6,*) 'stop'
         stop
      endif

      read(lucin)ityp		!types 
      write(6,*)'ityp'
      do i=1,ntyp
         na(i)=0
      enddo
      do i=1,im
         if ((ityp(i).lt.1).or.(ityp(i).gt.ntyp)) then
            write(6,*)'wrong ityp(',i,')= ',ityp(i)
            stop
         endif
         na(ityp(i))=na(ityp(i))+1
      enddo
      
      read(lucin)xp
      la=1
      lb=1
      lc=1
      if(icintypemod.eq.1) then
       read (lucin) num_at_glob
         read(lucin)xpp		!former positions
         read(lucin)vp          !velocities
         read(lucin)ax          !original positions
         read(lucin)oldtstep    !time step
      endif
      



      
      
      write (6,*)'-------- boite de simulation ------'
      write(6,*)'nombre d atomes =', im
      if (ltriclin) then
         write(6,*)' TRICLINIQUE '
      else
         write(6,*)' TETRAGONALE '
      endif
      write(6,'(A,3F11.4)')'taille de la boite ZL ',
     +     1d+08*zl(1),1d+08*zl(2),1d+08*zl(3)
      write(6,*)'nombre d atomes de chaque type'
      do iti=1,ntyp
         if (na(iti).ne.0) then
            write(6,*)na(iti),' atomes de type',iti
         endif
      enddo
      
      select case (nmod)
      case(1)
         write(6,*)'dilat ?'
         read(5,*)dilat
         icintype=0
         open(unit=18,file=fnamgout,form='formatted')
         write(18,*)la,lb,lc    !number of cells in 3 directions
         
         if (ltriclin) then
            at=at*1.d8*dilat
            write(6,*)'triclin'
            write(18,*)at(1,1),at(2,1),at(3,1) !a
            write(18,*)at(1,2),at(2,2),at(3,2) !b
            write(18,*)at(1,3),at(2,3),at(3,3) !c
            write(18,*)im       !number of atoms in UC
            call cryst_to_cart(imm,xp,bg,-1) !cart vers cryst
            do i=1,im
               do ic=1,3 
                  xp(ic,i)=xp(ic,i)
               enddo
            enddo
            do i=1,im
               write(18,'(3F18.11,I6)')xp(1,i),xp(2,i),xp(3,i),ityp(i) !coordonnes reduites des
            enddo
         else
            write(6,*)'PAS triclin'
            do ic=1,3
               zl(ic)=zl(ic)*1.0d8
            enddo
            write(18,*) zl(1)*dilat,'0.0 ','0.0 '
            write(18,*) '0.0 ',zl(2)*dilat,'0.0 '
            write(18,*) '0.0 ','0.0 ',zl(3)*dilat 
            write(18,*)im       !number of atoms in UC
            do i=1,im
               do ic=1,3 
                  xp(ic,i)=xp(ic,i)*1.0d8/zl(ic) +0.5
               enddo
               write(18,'(3F17.11,I2)')xp(1,i),xp(2,i),xp(3,i),ityp(i) !positions of atoms in UC (orthonormal coordinates assumed)
            enddo
         endif
      case(2)
         open(unit=18,file=fnamcout,form='unformatted')
         do i=1,im

            xt1=xp(1,i) ; xt2=xp(3,i)

            xp(1,i)=-1*xt2 ; xp(3,i)=xt1
            if(icintypemod.eq.1) then

               xt1=xpp(1,i) ; xt2=xpp(3,i)
               xpp(1,i)=-1*xt2 ; xpp(3,i)=xt1

               xt1=vp(1,i) ; xt2=vp(3,i)
               vp(1,i)=-1*xt2 ; vp(3,i)=xt1

               xt1=ax(1,i) ; xt2=ax(3,i)
               ax(1,i)=-1*xt2 ; ax(3,i)=xt1

            endif
         enddo
         normat(:)=at(:,1)
         at(:,1)=-1*at(:,3)
         at(:,3)=normat(:)
         
         write (18) icintype
         if (icintype<=1) then
            write (18) zl
         else
            write (18) at
         endif

         write (18) im
         
         write (18) ityp
         write (18) xp
         if (icintypemod==1) then
            write(18,*)num_at_glob
            write (18) xpp
            write (18) vp
            write (18) ax
            write(18)oldtstep
         endif
      endselect

      stop
      end
