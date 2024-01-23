! ************************************************
!         Sous-programme modecalc
! ************************************************

subroutine modecalc(im,xp,vp,ax)
  USE T_kind_param_m, ONLY:  double
  implicit none

    integer,intent(in)::im
    real(double),allocatable,intent(in),dimension(:,:)::xp,vp,ax
    
    integer :: i
    real(double)::sq1,asd,xnu,sca,scaa
    real(double),allocatable,save::dxp(:,:),eigval(:),xp_t(:,:)
    integer::iu,ju,ir,il,iut,j,k1,k2,k3,l1,l2,l3
    integer,save::icall=0
    
    icall=icall+1
    if (icall==1) then
       allocate (dxp(3*im,3*im))
       allocate(eigval(3*im))
       allocate(xp_t(3,3*im))
       
       open(50,file=&
               & 'vecteurs.dat',status='unknown') ! Choix du mode a exciter

     write(*,*) 'Loading vecteurs.dat'

     do ir=1,3*im

        sq1=0
        do iu=1,im
           do ju=1,3
              read(50,'(e14.6)') asd
              dxp(iu+(ju-1)*im,ir)=asd
              sq1=sq1+asd**2
           enddo
        enddo

        sq1=sqrt(sq1)
        do iu=1,3*im
           dxp(iu,ir)=dxp(iu,ir)/sq1
        enddo
     enddo
     write(6,*)'post vecteurs'
     close(50)
     open(51,file='eigenValuesREAL.dat',&
     &       status='unknown')
     open(52,file='eigenValuesIMAG.dat',&
     &       status='unknown')
     il=0
     do ir=1,3*im
        eigVal(ir)=0.d0
     enddo

     do ir=1,3*im
        read(51,*,END=36) iut,xnu
        if(xnu.gt.1) eigVal(iut)=xnu*1.d12
        il=il+1
     enddo
36   continue
     close(51)

     write(6,*)'post eigval'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! je choisis 3 modes parmis le total , ici les modes 10, im+10 et 3im-10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     DO I=1,IM
        do J=1,3

           K1=10
           L1=1
           XP_t(J,I)=dxp((J-1)*IM+I,(L1-1)*IM+K1)

           K2=10
           L2=2
           XP_t(J,I+im)=dxp((J-1)*IM+I,(L2-1)*IM+K2)

           K3=im-10
           L3=3
           XP_t(J,I+im*2)=dxp((J-1)*IM+I,(L3-1)*IM+K3)

        ENDDO
     ENDDO

     open(82,file='Freq_selec.dat',status='unknown')
     write(82,*) eigVal((L1-1)*im+K1)
     write(82,*) eigVal((L2-1)*im+K2)
     write(82,*) eigVal((L3-1)*im+K3)
     close(82)
  endif


       sca=0.
       scaa=0.
       do j=1,3
          do i=1,im
             sca=sca+xp_t(j,i)*VP(j,i)*1e8
             scaa=scaa+xp_t(j,i)*(xp(j,i)-Ax(j,i))*1e8
          enddo
       enddo
!       write(6,*)'P1'


       open(42,file='sca.dat',status='unknown',position='append')
       write(42,*) sca,scaa
       close(42)



       sca=0.
       scaa=0.
       do j=1,3
          do i=1,im
             sca=sca+xp_t(j,i+im)*VP(j,i)*1e8
             scaa=scaa+xp_t(j,i+im)*(xp(j,i)-Ax(j,i))*1e8
          enddo
       enddo

       open(42,file='scb.dat',status='unknown',position='append')
       write(42,*) sca,scaa
       close(42)




       sca=0.
       scaa=0.
       do j=1,3
          do i=1,im
             sca=sca+xp_t(j,i+im*2)*VP(j,i)*1e8
             scaa=scaa+xp_t(j,i+im*2)*(xp(j,i)-Ax(j,i))*1e8
          enddo
       enddo

       open(42,file='scc.dat',status='unknown',position='append')
       write(42,*) sca,scaa
       close(42)


  
  return


end subroutine modecalc
