module calccoordo_mod
  USE T_kind_param_m, ONLY:  double
  USE gen_com_m, ONLY: rang,iteration,timel,lperiod
  use atomconfig,only: atom_config
  use boxconfig,only:box_config
  USE cellconfig,only:cell_config
  USE var_pot, ONLY:ntyp!,ty !nkmax,ntyp,digr,gdertot
  use vect_dist_mod,only:vect_dist

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
    USE Tpara,only:MPI_COMM_space,status,ierr,NDM_MPI_REAl_DOUBLE

#endif

  implicit none
    class(atom_config),intent(in)::atcf
    class(box_config),intent(in)::boxcf
    type(cell_config),intent(in):: celcf
    integer,allocatable::na(:)
  integer :: i, iti,  i1, i2, koo, ko1, j,nci,ll
  real(double), dimension(ntyp,ntyp) :: dnco

#ifdef PARA
  real(double), dimension(ntyp,ntyp) :: dnco_glob

#endif

  real(double),allocatable :: rccoordo(:)

  integer,save::icall=0

  logical::linter

  allocate(rccoordo(ntyp))
  icall=icall+1
  if (icall==1) then
     rccoordo(1:ntyp)=rclu(1:ntyp)*1d-8
  end if
     allocate (na(ntyp))

    do i=1,atcf%im
       na(atcf%ityp(i))= na(atcf%ityp(i))+1
    end do
 !

  !      write(6,*)'entree calcoordo'
  write (6, *)
  write (6, *) '--------- Coordinations ----------------'
  dnco(:ntyp,:ntyp) = 0


  do i = 1,atcf%im
     koo = atcf%ielat(i)
     nci=0
     iti=atcf%ityp(i)
     do i1 = 0, 26
        ko1 = celcf%ncel(koo,i1)

        do i2 = 1, celcf%nato(ko1)
           j = celcf%atincel(i2,ko1)
           if (i==j) cycle

           call vect_dist(atcf,celcf,boxcf,i,j,rum=rccoordo(atcf%ityp(i)),linter=linter,lperiod=lperiod)
           if (linter)then
              dnco(atcf%ityp(i),atcf%ityp(j)) = dnco(atcf%ityp(i),atcf%ityp(j))+1
              nci=nci+1
           end if
        end do
     end do
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
     write (6, '(A,I5,A,D10.3)') '*  ITERATION  = ', iteration, '  time = ', timel
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


  return

end subroutine calccoordo
end module calccoordo_mod
