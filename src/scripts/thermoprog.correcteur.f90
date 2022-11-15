program thermoprog
  implicit none
  integer, parameter :: double = 8
  integer::i,im3,idum
  real(double), allocatable::eigval(:)
  integer::ntemp,nT
  real(double)::temperature,dtemp,hz_to_ev,temperature_to_ev,kbt,kbt2,xp2,eig,freqlim,two_pi
  real(double),allocatable,dimension(:)::Fmin,Smin,Fcla,Scla
  integer::jmat,nskip
  logical::lskip
  character*80::fnamthout
  freqlim =5d11
  read(5,*)im3
  allocate (eigval(im3))
!lit ds fréquences sans le facteuer /2PI et corrige
  nskip=0
  ntemp=300

  allocate(Fmin(ntemp))
  allocate(Fcla(ntemp))
  allocate(Smin(ntemp))
  allocate(Scla(ntemp))
  Fmin=0;Smin=0;Fcla=0;Scla=0
  hz_to_ev=4.135665538536d-15
  temperature_to_ev=8.6173303E-05


  dtemp=10.

  two_pi=2.0d0*4.d0*datan(1.d0)
  nskip=0
  do jmat=1,im3
     read(5,*)idum,eig
     
     lskip=.false.
     if ( (eig < 0).or.(dabs(eig)<freqlim)) then
        nskip=nskip+1
        cycle
     else
        eig=eig/two_pi

        do nT=1,ntemp
           temperature=nt*dtemp
           kBT=temperature*temperature_to_ev
           kBT2=2.0d0*kBT
           xp2=eig*hz_to_ev/kBT2
           Fmin(nT)=Fmin(nT)+kBT*DLOG(2.0d0*DSINH(xp2) )
           Smin(nT)=Smin(nT)+(xp2/DTANH(xp2)-DLOG(2.0d0*DSINH(xp2)))
           Fcla(nT)=Fcla(nT)+kBT*DLOG(2.0d0*xp2 )
           Scla(nT)=Scla(nT)-(DLOG(2.d0*xp2)+1.0d0)
           !             write(6,'(I5,f12.5, 5e25.15)') jmat,temperature, Fmin(nT), Fcla(nT), Smin(nT), Scla(nT),xp2

        end do
     end if  !jmat
  end do
  write(6,*)'NSKIP',nskip,im3-nskip
  fnamthout = 'PROG.thermo.dat'
  open(unit=123,file= fnamthout, form='formatted', status='unknown')
  do nT=1,ntemp
     temperature=nT*dtemp
     write(123,'(f12.5, 4e25.15)') temperature, Fmin(nT), Fcla(nT), Smin(nT), Scla(nT)
  end do
  return
  !    close(123)
end program thermoprog
