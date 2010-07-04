subroutine correlvp(xp, xpp, vp, ax, fp, ityp)
  use gen_com_m
  integer , intent(in)  :: ityp(imm)
  real(double) , intent(in)  :: xp(3,imm)
  real(double) , intent(in)  :: xpp(3,imm)
  real(double) , intent(in)  :: vp(3,imm)
  real(double) , intent(in)  :: ax(3,imm)
  real(double) , intent(in)  :: fp(3,imm)

  real(double),save :: mvp0vp0
  real(double) :: mvptvp0
  real(double) :: ratiocor
  real(double),pointer ::vptvp0(:),ratiocortyp(:)
  real(double),pointer,save ::  vp0vp0(:) 
  integer :: i,j,k,l,ic,m,n,iti,unitch
  real(double) :: cmat
  character :: chit*2,fch*80

  if (it==0) then
     if(rang==0) then
        open(unit=64,file='autocor')
        !        do iti=1,ntyp
        !           if(na(iti).ne.0) then
        !              unitch=12+iti
        !              write(6,*)'corvp',ty(iti)
        !              fch='autocor.'//ty(iti)
        !              open(unit=unitch,file=fch)
        !           end if
        !        end do
     end if

     !     vp0(:,:)=vp(:,:)
     allocate(vptvp0(ntyp)) ; allocate(vp0vp0(ntyp)) ; allocate(ratiocortyp(ntyp))
     mvp0vp0=0. ; vp0vp0(:)=0.
     do i=1,im
        cmat=cm(ityp(i))
        iti=ityp(i)
        do ic=1,3
           mvp0vp0=mvp0vp0+cmat*ax(ic,i)**2
           vp0vp0(iti)=vp0vp0(iti)+cmat*ax(ic,i)**2
           !           vp0vp0(ntyp+1)=vp0vp0(ntyp+1)+ax(ic,i)**2
           !  write(6,*)i,cmat,ax(ic,i)
        end do
     end do

     ratiocor=1.0
     if(rang==0)     write(64,'(I8,F11.7)')it,ratiocor
     return
  else
     mvptvp0=0. ; vptvp0(:)=0.0
     do i=1,im
        cmat=cm(ityp(i))
        iti=ityp(i)
        do ic=1,3
           mvptvp0=mvptvp0+cmat*vp(ic,i)*ax(ic,i)
           vptvp0(iti)=vptvp0(iti)+cmat*vp(ic,i)*ax(ic,i)
           !          vptvp0(3)=vptvp0(3)+vp(ic,i)*ax(ic,i)
           !        write(6,*)i,cmat,vp(ic,i),ax(ic,i)
        end do
     end do
     ratiocor=mvptvp0/mvp0vp0

     do iti=1,ntyp
        if(na(iti).ne.0) then
           unitch=12+iti
           ratiocortyp(iti)=vptvp0(iti)/vp0vp0(iti)
           !          if(rang==0) write(unitch,'(I7,F11.7)')it,ratiocortyp(iti)
        end if
     end do

     !        write(6,*)mvptvp0,mvp0vp0
     if(rang==0)     write(64,'(I7,D12.5,F11.7)')it,timel,ratiocor
     return
  end if



end subroutine correlvp
