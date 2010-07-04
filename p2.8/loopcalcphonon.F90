! ************************************************
!           Sous-programme dmloop.f
!          Version MPI du 21 fevrier 2001
! ************************************************

subroutine loopcalcphonon(xp, xpp, vp, ax, fp, ielat, iwmax, ityp)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  implicit none
  !-----------------------------------------------
  !   G l o b a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer  :: ielat(imm)
  integer  :: iwmax(imm)
  integer  :: ityp(imm)
  real(double)  :: xp(3,imm)
  real(double)  :: xpp(3,imm)
  real(double)  :: vp(3,imm)
  real(double)  :: ax(3,imm)
  real(double)  :: fp(3,imm)
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, iti,ic,is,j,u,v,icj,i3m
  real(double) :: unitE,unitP
  character*4 :: cunitE, cunitP
  real(double) :: epot0,deltaE,fps(3,imm)

  real(double),pointer :: matfor(:,:),fp0(:,:),d(:),d2(:),matforsym(:,:)
  !-----------------------------------------------
  !
  !
  ! calcul des phonons
  ! construit la matrice dynamique par calcul des forces pour de petits deplacements des atomes autour de positions
  ! d'equilibre. 
  ! diagonalise cette matrice. La diagonalisation sort les carre des pulsations
  ! les frequences sont les racines carrées divisées par 2pi des valeurs propres
  !
  !
  !
  i3m=3*im
  allocate (matfor(3*im,3*im))
  allocate (matforsym(3*im,3*im))
  allocate(fp0(3,imm))
  allocate(d(i3m))
  allocate(d2(i3m))
  unitE=1.0
  cunitE=' erg'

  unitP=1.0
  cunitP='d/cm2'


  ! MPI
  if (rang==0) write (6, *) '***** calcul des fréquences de phonons  ****'

  imd = im
  nad(:ntyp) = na(:ntyp)










  ! appel de la routine generale des forces
  !position de départ
  vp=0.0
  !      do i=1,im
  !         write(6,*)i,xp(1,i),xp(2,i),xp(3,i)
  !      end do
  fp=0.
  call calfo
  write (6, '(A,D21.12)') '*Epot = ', potist
  epot0=potist
  write(6,*)'forces ; doivent etre nulles'
  fp0=fp
  do i=1,im
     write(6,'(I2,3D21.12)')i,fp0(1,i),fp0(2,i),fp0(3,i)
  end do
  !      fps=fp


  do i=1,im
     do ic=1,3
        u=3*(i-1)+ic
        is=1
        !               do is=1,-1,-2

        xp(ic,i)=xp(ic,i)+is*deltax
        !            write(6,*)
        !            write(6,*) 'i,X is', i, ic,is

        call calfo

        do j=1,im
           do icj=1,3
              v=3*(j-1)+icj
              matfor(u,v)=-1.*(fp(icj,j)-fp0(icj,j))/(dsqrt(cm(ityp(i))*cm(ityp(j)))*deltax)
              !                  write(6,*)i,j,ic,icj
              !                  write(6,*)matfor(u,v),fp(icj,j)-fp0(icj,j)
           end do
        end do

        xp(ic,i)=xp(ic,i)-is*deltax
        !               end do
     end do
  end do


  !close(1)

  do i=1,i3m
     do j=1,i3m
        matforsym(i,j)=0.5*(matfor(i,j)+matfor(j,i))
     end do
  end do

  !      do i=1,i3m
  !         do j=1,i3m
  !            write(6,*)i,j,matfor(i,j),matforsym(i,j)
  !         end do
  !      end do

  call tred2(matforsym,i3m,i3m,D,D2) 
  call tqli(D,D2,i3m,i3m,matforsym)

  call EIGSRT(D,matforsym,i3m,i3m)
  write(6,*)'carré des pulsations'
  do i=1,i3m
     write(6,*)i,d(i)
  end do
  write(6,*)
  !         pi=3.141592654d0
  write(6,*)'fréquences'
  do i=1,i3m
     write(6,*)sqrt(d(i))/(2.*pi),' 1'
  end do

  stop
  return


end subroutine loopcalcphonon





SUBROUTINE Diag(A,N,NP,D,V) 
  implicit real*8 (a-h,o-z) 
  real*8 A(NP,NP),D(NP),V(NP,NP) 
  real*8 D2(NP)
  integer :: i,j 
  call tred2(A,NP,NP,D,D2) 
  call tqli(D,D2,NP,NP,A)
  ! call DCSMAA(A,NP,NP,D,D2,IERR)
  do i=1,NP 
     do j=1,NP 
        V(i,j)=a(i,j) 
     enddo
  enddo
  call EIGSRT(D,V,NP,NP)
END SUBROUTINE Diag

SUBROUTINE tqli(d,e,n,np,z) 
  INTEGER n,np 
  DOUBLE PRECISION d(np),e(np),z(np,np) 
  !U USES pythag 
  INTEGER i,iter,k,l,m 
  DOUBLE PRECISION b,c,dd,f,g,p,r,s,pythag,xprec

  do 11 i=2,n 
     e(i-1)=e(i) 
11   continue 
     e(n)=0.d0
     do 15 l=1,n 
        iter=0 
1       do 12 m=l,n-1 
           dd=dabs(d(m))+dabs(d(m+1)) 
           if (dabs(e(m))+dd.eq.dd) goto 2 
12         continue 
           m=n 
2          if(m.ne.l)then 
              if(iter.gt.300)then 
                 write(*,*) 'too many iterations in tqli' 
                 return 
              endif
              iter=iter+1 
              g=(d(l+1)-d(l))/(2.d0*e(l)) 
              r=pythag(g,1.d0) 
              g=d(m)-d(l)+e(l)/(g+dsign(r,g)) 
              s=1.d0 
              c=1.d0 
              p=0.d0 
              do 14 i=m-1,l,-1 
                 f=s*e(i) 
                 b=c*e(i) 
                 r=pythag(f,g) 
                 e(i+1)=r 
                 if(r.eq.0.d0)then 
                    d(i+1)=d(i+1)-p 
                    e(m)=0.d0 
                    goto 1 
                 endif
                 s=f/r 
                 c=g/r 
                 g=d(i+1)-p 
                 r=(d(i)-g)*s+2.d0*c*b 
                 p=s*r 
                 d(i+1)=g+p 
                 g=c*r-b 
                 ! Omit lines from here ... 
                 do 13 k=1,n 
                    f=z(k,i+1) 
                    z(k,i+1)=s*z(k,i)+c*f 
                    z(k,i)=c*z(k,i)-s*f 
13                  continue 
                    ! ... to here when finding only eigenvalues. 
14                  continue 
                    d(l)=d(l)-p 
                    e(l)=g 
                    e(m)=0.d0 
                    goto 1 
                 endif
15               continue
                 return 
              END do

              SUBROUTINE tred2(a,n,np,d,e) 
                INTEGER n,np 
                DOUBLE PRECISION a(np,np),d(np),e(np) 
                INTEGER i,j,k,l 
                DOUBLE PRECISION f,g,h,hh,scale
                do 18 i=n,2,-1 
                   l=i-1 
                   h=0.d0 
                   scale=0.d0 
                   if(l.gt.1)then 
                      do 11 k=1,l 
                         scale=scale+dabs(a(i,k)) 
11                       continue 
                         if(scale.eq.0.d0)then 
                            e(i)=a(i,l) 
                         else 
                            do 12 k=1,l 
                               a(i,k)=a(i,k)/scale 
                               h=h+a(i,k)**2 
12                             continue 
                               f=a(i,l) 
                               g=-dsign(dsqrt(h),f) 
                               e(i)=scale*g 
                               h=h-f*g 
                               a(i,l)=f-g 
                               f=0.d0 
                               do 15 j=1,l 
                                  ! Omit following line if finding only eigenvalues 
                                  a(j,i)=a(i,j)/h 
                                  g=0.d0 
                                  do 13 k=1,j 
                                     g=g+a(j,k)*a(i,k) 
13                                   continue 
                                     do 14 k=j+1,l 
                                        g=g+a(k,j)*a(i,k) 
14                                      continue 
                                        e(j)=g/h 
                                        f=f+e(j)*a(i,j) 
15                                      continue 
                                        hh=f/(h+h) 
                                        do 17 j=1,l 
                                           f=a(i,j) 
                                           g=e(j)-hh*f 
                                           e(j)=g 
                                           do 16 k=1,j 
                                              a(j,k)=a(j,k)-f*e(k)-g*a(i,k) 
16                                            continue 
17                                            continue 
                                           endif 
                                        else 
                                           e(i)=a(i,l) 
                                        endif 
                                        d(i)=h 
18                                      continue 
                                        ! Omit following line if finding only eigenvalues. 
                                        d(1)=0.d0 
                                        e(1)=0.d0 
                                        do 24 i=1,n 
                                           ! Delete lines from here ... 
                                           l=i-1 
                                           if(d(i).ne.0.d0)then 
                                              do 22 j=1,l 
                                                 g=0.d0 
                                                 do 19 k=1,l 
                                                    g=g+a(i,k)*a(k,j) 
19                                                  continue 
                                                    do 21 k=1,l 
                                                       a(k,j)=a(k,j)-g*a(k,i) 
21                                                     continue 
22                                                     continue 
                                                    endif 
                                                    ! ... to here when finding only eigenvalues. 
                                                    d(i)=a(i,i) 
                                                    ! Also delete lines from here ... 
                                                    a(i,i)=1.d0 
                                                    do 23 j=1,l 
                                                       a(i,j)=0.d0 
                                                       a(j,i)=0.d0 
23                                                     continue 
                                                       ! ... to here when finding only eigenvalues. 
24                                                     continue 
                                                       return 
                                                    END do

                                                    FUNCTION pythag(a,b) 
                                                      DOUBLE PRECISION a,b,pythag 
                                                      DOUBLE PRECISION absa,absb 
                                                      absa=dabs(a) 
                                                      absb=dabs(b) 
                                                      if(absa.gt.absb)then 
                                                         pythag=absa*dsqrt(1.d0+(absb/absa)**2) 
                                                      else 
                                                         if(absb.eq.0.d0)then 
                                                            pythag=0.d0 
                                                         else 
                                                            pythag=absb*dsqrt(1.d0+(absa/absb)**2) 
                                                         endif
                                                      endif
                                                      return 
                                                    END FUNCTION pythag
                                                    ! (C) Copr. 1986-92 Numerical Recipes Software 'k'1k30m,t+W.

                                                    SUBROUTINE EIGSRT(D,V,N,NP) 
                                                      integer N,NP 
                                                      Real*8 D(NP),V(NP,NP),P 
                                                      DO 13 I=1,N-1 
                                                         K=I 
                                                         P=D(I) 
                                                         DO 11 J=I+1,N 
                                                            IF(D(J).LE.P)THEN 
                                                               K=J 
                                                               P=D(J) 
                                                            ENDIF
11                                                          CONTINUE 
                                                            IF(K.NE.I)THEN 
                                                               D(K)=D(I) 
                                                               D(I)=P 
                                                               DO 12 J=1,N 
                                                                  P=V(J,I) 
                                                                  V(J,I)=V(J,K) 
                                                                  V(J,K)=P 
12                                                                CONTINUE 
                                                               ENDIF 
13                                                             CONTINUE 
                                                               RETURN 
                                                            END IF
