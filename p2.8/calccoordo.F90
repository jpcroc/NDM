module calccoordo_mod
  USE T_kind_param_m, ONLY:  double
  USE notperiod_mod,only: notperiod
  USE cryst_to_cart_mod,only: cryst_to_cart
  USE gen_com_m, ONLY: rang,it,timel,lperiod
  use atomconfig,only: atom_config
  use boxconfig,only:box_config
  USE cellconfig,only:cell_config, caltabtC
  USE var_pot, ONLY:ntyp!,ty !nkmax,ntyp,digr,gdertot

  implicit none
  real(double)::rclu(20)
contains
  subroutine calccoordo(atcf,celcf,boxcf)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m, ONLY:  double
#ifdef PARA
    USE mpi
    USE Tpara,only:MPI_COMM_space,status,ierr,nprocs,myidsp,NDM_MPI_REAl_DOUBLE

#endif

  implicit none
    class(atom_config),intent(in)::atcf
    type(box_config),intent(in)::boxcf
    type(cell_config),intent(in):: celcf
    integer,allocatable::na(:)
  integer :: i, iti, itj, i1, i2, koo, ko1, j, ic,nci,ip,ll
  real(double), dimension(ntyp,ntyp) :: dnco
  real(double) :: r2, c1, c2, c3,cv(1,3),x1,x2,x3
  real(double):: a1, a2, a3

#ifdef PARA
  real(double), dimension(ntyp,ntyp) :: dnco_glob
  integer :: nci_glob
#endif

  real(double),allocatable :: xpnp(:,:),rccoordo(:)

  integer,save::icall=0
  icall=icall+1
  if (icall==1) then
     rccoordo(1:ntyp)=rclu(1:ntyp)*1d-8
  end if
     allocate(xpnp(3,atcf%imm))
     allocate (na(ntyp))

    do i=1,atcf%im
       na(atcf%ityp(i))= na(atcf%ityp(i))+1
    end do
 !

  !      write(6,*)'entree calcoordo'
  write (6, *)
  write (6, *) '--------- Coordinations ----------------'
  dnco(:ntyp,:ntyp) = 0



  if (lperiod) then
     xpnp(:,:)=atcf%xp(:,:)
  else 
     call notperiod(atcf%imm,atcf%xp,xpnp,boxcf%at,boxcf%bg)
  end if

  do i = 1,atcf%im
     koo = atcf%ielat(i)
     nci=0
     iti=atcf%ityp(i)
     do i1 = 0, 26
        ko1 = celcf%ncel(koo,i1)

        do i2 = 1, celcf%nato(ko1)
           j = celcf%atincel(i2,ko1)
           if (i==j) cycle
           c1 = xpnp(1,i)-xpnp(1,j)
           c2 = xpnp(2,i)-xpnp(2,j)
           c3 = xpnp(3,i)-xpnp(3,j)
              c1 = c1+sum(boxcf%at(1,:)*celcf%deltadist(:,i1,koo))
              c2 = c2+sum(boxcf%at(2,:)*celcf%deltadist(:,i1,koo))
              c3 = c3+sum(boxcf%at(3,:)*celcf%deltadist(:,i1,koo))
           if (celcf%noxyz.eq.1) then
              cv(1,1) = c1
              cv(1,2) = c2
              cv(1,3) = c3
              call cryst_to_cart (1, cv, boxcf%bg, -1) !cryst vers cart sur cv
              WHERE ( (cv.GT.0.5d0).OR.(cv.LT.-0.5d0) )
                 cv(:,1:3) = cv(:,1:3) - Dble(Nint(cv(:,1:3)))
              END WHERE
              call cryst_to_cart (1, cv, boxcf%at, 1) !cryst vers cart sur cv
              c1=cv(1,1)
              c2=cv(1,2)
              c3=cv(1,3)
              
           end if
              
           r2 = sqrt(c1*c1+c2*c2+c3*c3)

           !               if (r2.le.3.0d-15)               write(6,*)i,j,r2

           if (r2<rccoordo(atcf%ityp(i)))then
              dnco(atcf%ityp(i),atcf%ityp(j)) = dnco(atcf%ityp(i),atcf%ityp(j))+1
              nci=nci+1
           end if
           !               if (r2>=rc(ityp(j))) cycle
           !               dnco(ityp(j),ityp(i)) = dnco(ityp(i),ityp(j))+1

        end do
     end do
           !   write(6,*)i,nci

  end do

#ifdef PARA
  call MPI_ALLREDUCE(dnco(:,:),dnco_glob(:,:),ntyp*ntyp,NDM_MPI_REAL_DOUBLE,MPI_SUM,MPI_COMM_space,ierr)
  dnco = dnco_glob
  ! nci n'est pas utilise dans la suite, je laisse en commentaire la reduction
  !  call MPI_ALLREDUCE(nci,nci_glob,1,MPI_INTEGER,MPI_SUM,MPI_COMM_space,ierr)
  !  nci = nci_glob
#endif

  do i1 = 1, ntyp
     if (na(i1)==0) cycle
     where (na(:ntyp)/=0) dnco(i1,:ntyp) = dnco(i1,:ntyp)/na(i1)
  enddo
  if (rang==0) then
     write (6, '(A,I5,A,D10.3)') '*  ITERATION  = ', it, '  time = ', timel
     write (6, *) 'nombre de coor autour de :'
     do i1 = 1, ntyp
        if (na(i1)==0) cycle
        write (6, *) i1, (dnco(i1,j),j=1,ntyp)
        !   28 FORMAT (31H Nombre de coor autour de ITYP=,I2,2X,7(F5.2,1X))
     end do
     do i1 = 1, ntyp
        if (na(i1)==0) cycle
        write (6, *) 'rayon autour des type ', i1, ' = ', rccoordo(i1)*1D+8
     end do
  endif

  deallocate (xpnp)
  return

end subroutine calccoordo
end module
