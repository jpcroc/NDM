module bondval_mod
  USE notperiod_mod,only: notperiod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE gen_com_m, ONLY:at,bg,it,atincel,ncel,rang,fnam,deltadist,lperiod,noxyz,nato,lenfnam
  USE var_pot, ONLY:ntyp,ty
  implicit none

contains

  subroutine bondval(im,xp,ityp,ielat,num_at_glob)

    USE T_kind_param_m, ONLY:  double

    integer,intent(in)::im
    real(double),intent(in),allocatable::xp(:,:)
    integer,allocatable::ityp(:),ielat(:),num_at_glob(:)

    ! variables locales

    integer :: i,k,i1,j,iti,itj,i2
    integer :: ncelvois,koo,ko1,lenfn2

    real(double)::c1p,c2p,c3p,cv(1,3),ra(3),c1,c2,c3

    real(double),allocatable::bdv(:)
    real(double) :: R,xx,dcut2,dis

    real(double),allocatable :: xpnp(:,:)

    character :: extension*9
    !APARA
    allocate(xpnp(3,im))

    if(rang==0) then

       lenfn2 = 9
       write(extension,'(i9.9)') it

       ! -------------------------------------------------------------
       !     creation du  fichier positions pour Rasmol
       ! -------------------------------------------------------------
       open(65, file=fnam(1:lenfnam)//'.'//extension(1:lenfn2)//'.bdv', form='formatted', &
            status='unknown')


    end if

    dcut2=(5.0d-8)**2

    allocate(bdv(im))


    bdv(:)=0.

    if (lperiod) then
       xpnp(:,:)=xp(:,:)
    else 
       call notperiod(im,xp,xpnp)
    end if



    do i = 1, im
       koo = ielat(i)                          ! Numero de la cellule
       iti=ityp(i)
       ncelvois = min(noxyz,27)-1

       do i1 = 0, ncelvois
          ko1 = ncel(koo,i1)

          c1p = xpnp(1,i)+sum(at(1,:)*deltadist(:,i1,koo))
          c2p = xpnp(2,i)+sum(at(2,:)*deltadist(:,i1,koo))
          c3p = xpnp(3,i)+sum(at(3,:)*deltadist(:,i1,koo))
          ! pour chaque atome ds la cel. voisine
          do i2 = 1, nato(ko1)
             j = atincel(i2,ko1)

             itj=ityp(j) 
             if (iti.eq.itj) cycle 

             if((trim(ty(iti)).ne.'O').and.((trim(ty(itj)).ne.'O'))) cycle
             c1 = c1p-xpnp(1,j)
             c2 = c2p-xpnp(2,j)
             c3 = c3p-xpnp(3,j)
             if (noxyz.eq.1) then
                cv(1,1) = c1
                cv(1,2) = c2
                cv(1,3) = c3
                call cryst_to_cart (1, cv, bg, -1) !cryst vers cart sur cv
                WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                   cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
                END WHERE
                call cryst_to_cart (1, cv, at, 1) !cryst vers cart sur cv
                c1=cv(1,1)
                c2=cv(1,2)
                c3=cv(1,3)

             end if

             !           c1 = c1p-xpnp(1,j)
             !           c2 = c2p-xpnp(2,j)
             !           c3 = c3p-xpnp(3,j)

             dis = c1*c1+c2*c2+c3*c3


             if (dis>dcut2) cycle
             dis=sqrt(dis)



             R = 0.0


             ! selection de la constante pour chaque atome central


             IF ( trim(ty(ityp(j)))  .EQ.  'H' .OR. trim(ty(ityp(i))) .EQ. 'H') &
                  &        R = 0.95
             IF ( trim(ty(ityp(j)))  .EQ.  'Ba' .OR. trim(ty(ityp(i))) .EQ. 'Ba') &
                  &        R = 2.285
             IF ( trim(ty(ityp(j)))  .EQ.  'Na' .OR. trim(ty(ityp(i))) .EQ. 'Na') &
                  &        R = 1.803
             IF ( trim(ty(ityp(j)))  .EQ.  'Rb' .OR. trim(ty(ityp(i))) .EQ. 'Rb') &
                  &        R = 2.26
             IF ( trim(ty(ityp(j)))  .EQ.  'Cs' .OR. trim(ty(ityp(i))) .EQ. 'Cs') &
                  &        R = 2.42
             IF ( trim(ty(ityp(j)))  .EQ.  'Sr' .OR. trim(ty(ityp(i))) .EQ. 'Sr') &
                  &        R = 2.118
             IF ( trim(ty(ityp(j)))  .EQ.  'B' .OR. trim(ty(ityp(i))) .EQ. 'B') &
                  &        R = 1.371
             IF ( trim(ty(ityp(j)))  .EQ.  'P'  .OR. trim(ty(ityp(i))) .EQ. 'P') &
                  &         R = 1.617
             IF ( trim(ty(ityp(j)))  .EQ.  'Si' .OR. trim(ty(ityp(i))) .EQ. 'Si') &
                  &        R = 1.624
             IF ( trim(ty(ityp(j)))  .EQ.  'Al' .OR. trim(ty(ityp(i))) .EQ. 'Al') &
                  &        R = 1.651
             IF ( trim(ty(ityp(j)))  .EQ.  'Ca' .OR. trim(ty(ityp(i))) .EQ. 'Ca') &
                  &        R = 1.967
             IF ( trim(ty(ityp(j)))  .EQ.  'K'  .OR. trim(ty(ityp(i))) .EQ. 'K') &
                  &         R = 2.132
             IF ( trim(ty(ityp(j)))  .EQ.  'Li' .OR. trim(ty(ityp(i))) .EQ. 'Li') &
                  &        R = 1.466
             IF ( trim(ty(ityp(j)))  .EQ.  'Mg' .OR. trim(ty(ityp(i))) .EQ. 'Mg') &
                  &        R = 1.693
             IF ( trim(ty(ityp(j)))  .EQ.  'Y'  .OR. trim(ty(ityp(i))) .EQ. 'Y') &
                  &         R = 1.965
             IF ( trim(ty(ityp(j)))  .EQ.  'Zr' .OR. trim(ty(ityp(i))) .EQ. 'Zr') &
                  &        R = 1.937
             IF ( trim(ty(ityp(j)))  .EQ.  'La' .OR. trim(ty(ityp(i))) .EQ. 'La') &
                  &        R = 2.172
             IF ( trim(ty(ityp(j)))  .EQ.  'Th' .OR. trim(ty(ityp(i))) .EQ. 'Th') &
                  &        R = 2.18
             IF ( trim(ty(ityp(j)))  .EQ.  'Ni' .OR. trim(ty(ityp(i))) .EQ. 'Ni') &
                  &        R = 1.654
             IF ( trim(ty(ityp(j)))  .EQ.  'Co2' .OR. trim(ty(ityp(i))) .EQ. 'Co2') &
                  &        R = 1.692
             IF ( trim(ty(ityp(j)))  .EQ.  'Co3' .OR. trim(ty(ityp(i))) .EQ. 'Co3') &
                  &        R = 1.70
             IF ( trim(ty(ityp(j)))  .EQ.  'Cl' .OR. trim(ty(ityp(i))) .EQ. 'Cl') &
                  &        R = 1.632
             IF ( trim(ty(ityp(j)))  .EQ.  'Be' .OR. trim(ty(ityp(i))) .EQ. 'Be') &
                  &        R = 1.381
             IF ( trim(ty(ityp(j)))  .EQ.  'Zn' .OR. trim(ty(ityp(i))) .EQ. 'Zn') &
                  &        R = 1.704
             IF ( trim(ty(ityp(j)))  .EQ.  'Ge' .OR. trim(ty(ityp(i))) .EQ. 'Ge') &
                  &        R = 1.748
             IF ( trim(ty(ityp(j)))  .EQ.  'Ga' .OR. trim(ty(ityp(i))) .EQ. 'Ga') &
                  &        R = 1.73
             IF ( trim(ty(ityp(j)))  .EQ.  'Cr3' .OR. trim(ty(ityp(i))) .EQ. 'Cr3') &
                  &        R = 1.724
             IF ( trim(ty(ityp(j)))  .EQ.  'Cr6' .OR. trim(ty(ityp(i))) .EQ. 'Cr6') &
                  &        R = 1.794

             IF ( trim(ty(ityp(j)))  .EQ.  'Nd' .OR. trim(ty(ityp(i))) .EQ. 'Nd') &
                  &        R = 2.117
             IF ( trim(ty(ityp(j)))  .EQ.  'Eu2' .OR. trim(ty(ityp(i))) .EQ. 'Eu2') &
                  &        R = 2.147
             IF ( trim(ty(ityp(j)))  .EQ.  'Eu3' .OR. trim(ty(ityp(i))) .EQ. 'Eu3') &
                  &        R = 2.076
             IF ( trim(ty(ityp(j)))  .EQ.  'Ce3' .OR. trim(ty(ityp(i))) .EQ. 'Ce3') &
                  &        R = 2.151
             IF ( trim(ty(ityp(j)))  .EQ.  'Ce4' .OR. trim(ty(ityp(i))) .EQ. 'Ce4') &
                  &        R = 2.028


             ! cas titanyl (Farges et al. 1996 - Part I) ; cutoff distance entre Ti=O et Ti-O = 1.84 Â¡
             IF ( trim(ty(ityp(j))) .EQ. 'Ti' .and. dis .lt. 1.84)   R = 1.89
             IF ( trim(ty(ityp(j))) .EQ. 'Ti' .and. dis .ge. 1.84)   R = 1.815
             IF ( trim(ty(ityp(i)))  .EQ. 'Ti' .and. dis .lt. 1.84)   R = 1.89
             IF ( trim(ty(ityp(i)))  .EQ. 'Ti' .and. dis .ge. 1.84)   R = 1.815


             ! Fe2+ et ensuite Fe3+
             IF ( trim(ty(ityp(j)))  .EQ.  'Fe2' .OR. trim(ty(ityp(i))) .EQ. 'Fe2') &
                  &        R = 1.734 
             IF ( trim(ty(ityp(j)))  .EQ.  'Fe3' .OR. trim(ty(ityp(i))) .EQ. 'Fe3') &
                  &        R = 1.759


             ! Mn2+, puis Mn3+ et Mn4+ et Mn7+
             IF ( trim(ty(ityp(j)))  .EQ.  'Mn2' .OR. trim(ty(ityp(i))) .EQ. 'Mn2') &
                  &        R = 1.79          
             IF ( trim(ty(ityp(j)))  .EQ.  'Mn3' .OR. trim(ty(ityp(i))) .EQ. 'Mn3') &
                  &        R = 1.76  
             IF ( trim(ty(ityp(j)))  .EQ.  'Mn4' .OR. trim(ty(ityp(i))) .EQ. 'Mn4') &
                  &        R = 1.753  
             IF ( trim(ty(ityp(j)))  .EQ.  'Mn7' .OR. trim(ty(ityp(i))) .EQ. 'Mn7') &
                  &        R = 1.79 


             ! As ici est As3+, puis 5+
             IF ( trim(ty(ityp(j)))  .EQ.  'As3' .OR. trim(ty(ityp(i))) .EQ. 'As3') &
                  &        R = 1.789 
             IF ( trim(ty(ityp(j)))  .EQ.  'As5' .OR. trim(ty(ityp(i))) .EQ. 'As5') &
                  &        R = 1.767


             ! Mo6+	 
             IF ( trim(ty(ityp(j)))  .EQ.  'Mo' .OR. trim(ty(ityp(i))) .EQ. 'Mo') &
                  &        R = 1.907 
             ! Ta5+
             IF ( trim(ty(ityp(j)))  .EQ.  'Ta' .OR. trim(ty(ityp(i))) .EQ. 'Ta') &
                  &        R = 1.92
             ! U6+, puis 4+ (valeur necessitant rafinnage + approfondi)
             IF ( trim(ty(ityp(j)))  .EQ.  'U6' .OR. trim(ty(ityp(i))) .EQ. 'U6') &
                  &        R = 2.075  
             IF ( trim(ty(ityp(j)))  .EQ.  'U4' .OR. trim(ty(ityp(i))) .EQ. 'U4') &
                  &        R = 2.112 	

             IF ( trim(ty(ityp(j)))  .EQ.  'W6' .OR. trim(ty(ityp(i))) .EQ. 'W6') &
                  &        R = 1.921 
             IF ( trim(ty(ityp(j)))  .EQ.  'Mo6' .OR. trim(ty(ityp(i))) .EQ. 'Mo6') &
                  &        R = 1.907

             IF ( trim(ty(ityp(j)))  .EQ.  'Pb2' .OR. trim(ty(ityp(i))) .EQ. 'Pb2') &
                  &        R = 2.112
             IF ( trim(ty(ityp(j)))  .EQ.  'Pb4' .OR. trim(ty(ityp(i))) .EQ. 'Pb4') &
                  &        R = 2.042

             IF ( trim(ty(ityp(j)))  .EQ.  'Au3' .OR. trim(ty(ityp(i))) .EQ. 'Au3') &
                  &        R = 1.833



             ! cas ou le symbole n a pas ete reconnu (R reste a zero)
             if ( R == 0 ) then
                write(6,*) rang,'Atome non reconnu : i j iti itj ', num_at_glob(i) ,num_at_glob(j),iti,itj
             endif

             dis=dis*1.0d8

             ! calcul de la force de liaison individuelle (xx) 
             xx = exp( ( R - dis ) / 0.37 )
             !               write(6,*)R,dis,xx



             ! ainsi que la force de liaison cumulee autour du cluster i
             bdv(i) = bdv(i) + xx
          end do
       end do
    end do
    j=0

    if(rang==0) then

       do iti=1,ntyp
          do i=1,im
             if(ityp(i).ne.iti)cycle
             j=j+1
             write(365,'(2I6,I3,F12.2)')j,i,iti,bdv(i)
          end do
       end do
       close(65)
    end if

    deallocate(bdv)
    deallocate (xpnp)

200 format(i2)
300 format(i3)
400 format(i4)
500 format(i5)
600 format(i6)
700 format(i7)
800 format(i8)

101 format(a1)
201 format(a2)
301 format(a3)
401 format(a4)
501 format(a5)
601 format(a6)
701 format(a7)
801 format(a8)
901 format(a9)
900 format(i9)




  end subroutine bondval
end module bondval_mod
