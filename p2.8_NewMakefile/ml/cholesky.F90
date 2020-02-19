module cholesky_interface 
        interface 

                subroutine init_my_matrix_ml (context,iam,nprocs,n,a,d_a,nprow,npcol,myrow,mycol,desca,descd_a)
                implicit none
                integer :: n
                integer :: context,iam,nprocs
                integer :: nprow,npcol,myrow,mycol
                integer :: desca(:), descd_a(:)
                real(kind=8)    :: a(:,:), d_a(:,:)
                end subroutine init_my_matrix_ml

        end interface

end module cholesky_interface 


module cholesky_mod
        implicit none
        contains
subroutine  solve_kernel_matrix(rangml)
use ml_in_ndm_module, only: isave_ml
USE T_kind_param_m, ONLY:  double
implicit none
integer, intent(in) :: rangml
integer :: count1,count2,count_rate,count_max
real(double) :: time

#ifdef PARAML


  call system_clock (count1,count_rate,count_max)

  if (isave_ml==2) then
     call cholesky_threading()
  else
     write(*,*) 'scal'
     call cholesky_scalapack()
  end if
#else
    call cholesky_threading()
#endif

  call system_clock (count2,count_rate,count_max)

  time=real((count2-count1))/real(count_rate)
  if (rangml==0)  write(*,"(' ML:< solve_kernel_matrix > inversion was  ................:  ',f16.8,' s')") time


return

end subroutine  solve_kernel_matrix





subroutine cholesky_threading()
USE T_kind_param_m, ONLY:  double
use ml_in_ndm_module
use temporary_data_cov, ONLY: yfunc_train, dim_train, &
                              mlocal, u_local,v_local, i_local_cov, &
                              kbfunc,matfor
#ifdef PARAML
!  include "mkl_service.h"
  use mpi
  use mod_mpi_ml

#endif
implicit none
integer :: i,j,u,v
!debug real(double), dimension(:,:), allocatable :: matfor
real(double)::time
integer :: count1,count2,count_rate,count_max    ! time counting.
!integer :: valeur, somme

!cholesky driver
integer :: info,nrhs
!end cholesky driver

!diago driver
real(double), dimension(:), allocatable :: work
real(double), dimension(:), allocatable :: w
character*1  :: uplo
integer :: lwork
!diago driver



if (allocated(kbfunc)) deallocate(kbfunc)
allocate(kbfunc(dim_train))

                       ! 0 through threading diagonalization.

! In threading part ... whtever happens matfor is allocated only in the
! node 0.
if (rangml == 0 ) then
   if (allocated(matfor)) deallocate(matfor)
   allocate(matfor(dim_train,dim_train))
end if


#ifdef PARAML

if (isave_ml==2) then  ! this correspond to the case when DYN matrix is
                       ! computed/readed by MPI and then is packed only to the procs number
                       ! 0 through threading diagonalization.
  if (rangml == 0 ) then
   matfor(:,:)=0.d0
   do i=1,i_local_cov
     u=u_local(i)
     v=v_local(i)
     matfor(u,v)=mlocal(i)
   end do
  end if
  if (rangml /= 0 ) then
       call MPI_SEND(mlocal,i_local_cov,MPI_DOUBLE_PRECISION,0,10000+rangml,MPI_COMM_WORLD,codeml)
       call MPI_SEND(u_local,i_local_cov,MPI_INTEGER,0,20000+rangml,MPI_COMM_WORLD,codeml)
       call MPI_SEND(v_local,i_local_cov,MPI_INTEGER,0,30000+rangml,MPI_COMM_WORLD,codeml)
  end if

  if (rangml == 0)  then
      iproc=0
      do iproc=1,nb_procsml-1
       call MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD,statut,codeml)
       call MPI_GET_COUNT(statut,MPI_DOUBLE_PRECISION,nb_elements,codeml)
       !debug write(*,*) 'on proc  ', iproc, 'imax is', nb_elements
       deallocate(mlocal,u_local,v_local)
       allocate(u_local(nb_elements),v_local(nb_elements), mlocal(nb_elements))
       itempproc=statut(MPI_SOURCE)
       call MPI_RECV(mlocal,nb_elements,MPI_DOUBLE_PRECISION,statut(MPI_SOURCE), 10000+statut(MPI_SOURCE),MPI_COMM_WORLD,statut,codeml)
       call MPI_RECV(u_local,nb_elements,MPI_INTEGER,itempproc,20000+itempproc,MPI_COMM_WORLD,statut,codeml)
       call MPI_RECV(v_local,nb_elements,MPI_INTEGER,itempproc,30000+itempproc,MPI_COMM_WORLD,statut,codeml)
       do i=1,nb_elements
         u=u_local(i)
         v=v_local(i)
         matfor(u,v)=mlocal(i)
       end do
      end do
   end if !rangml==0
   !
end if ! isave_ml=2
  call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#else


 if (iread_ml==0) then
    if (allocated(matfor)) deallocate(matfor)
    allocate(matfor(dim_train,dim_train))
    matfor(:,:)=0.d0
    do i=1,i_local_cov
      u=u_local(i)
      v=v_local(i)
      matfor(u,v)=mlocal(i)
    end do
 end if

#endif

  call system_clock (count2,count_rate,count_max)

#ifdef PARAML


 if (isave_ml==2) then

 if (rangml==0) then
#endif

! this symmetrization is done only in the case DYN_MAT - MPI + diago in threading
  do i=1,dim_train
    do j=i,dim_train

    matfor(i,j)=0.5d0*(matfor(i,j)+matfor(j,i))
    matfor(j,i)=matfor(i,j)
    end do
  end do
if (debug) then
   if (rangml==0) write(*,*) 'ML: MATFOR (1,2) (2,1)', matfor (1,2) ,matfor(2,1)
   if (rangml==0) write(*,*) 'ML: MATFOR (1,5) (5,1)', matfor (1,5) ,matfor(5,1)
end if

#ifdef PARAML


end if  !rang_ml==0

  call MPI_BARRIER(MPI_COMM_WORLD,codeml)
  call system_clock (count1,count_rate,count_max)
  time=real((count1-count2))/real(count_rate)
if (debug) then
  if (rangml==0) write(*,"(' ML: MATFOR  was symmetrized in...:  ',f16.8,' s')") time
end if

end if  !isave_ml=2
#else

  call system_clock (count1,count_rate,count_max)
  time=real((count1-count2))/real(count_rate)
if (debug) then
  if (rangml==0) write(*,"(' ML: MATFOR  was symmetrized in...:  ',f16.8,' s')") time
end if
#endif

call system_clock (count1,count_rate,count_max)

!do i=1,dim_train
!  write(61,'(10g12.4)') matfor(i,:)
!end do


if (rangml==0) then

 uplo='u'
 info=0
 lwork=max(1,3*dim_train-1)+dim_train

 if (allocated(work)) deallocate(work)
 allocate(work(lwork))
 if (allocated(w)) deallocate(w)
 allocate(w(dim_train))

 do i=1,dim_train
  matfor(i,i)=matfor(i,i)+lambda_krr
 end do

!call DSYEV( 'N', 'L', dim_train, matfor, dim_train, w, work, lwork, info )

!do i=1,dim_train
!   if (ik==1) then
!     write(75,*) i, w(i)
!   end if
!end do


 kbfunc(:)=yfunc_train(:)
!kbfunc(:)=yfunc(:)
 nrhs=1

 call  DPOSV( 'L', dim_train, nrhs, matfor, dim_train, kbfunc, dim_train, info)

 if (debug) then
  write(*,'("ML: kbfunc(1:4)", 4f20.7)') kbfunc(1:4)
 end if

! call  DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i of A is not
!                positive definite, so the factorization could not be
!                completed, and the solution has not been computed.

 if (info < 0 ) then
   if (rangml==0) write(6,*) ' if INFO = -i, the i-th argument had an illegal value', info
   if (rangml==0) write(6,*) 'stop in cholesky_factorization'
   stop
 end if
 if (info > 0 ) then
   if (rangml==0) write(6,*) 'if INFO = i, the leading minor of order i of A is not positive definite',info
   if (rangml==0) write(6,*) 'stop in cholesky_factorization'
   stop
 end if


end if  ! rangml=0



#ifdef PARAML

call MPI_BARRIER(MPI_COMM_WORLD,codeml)

!valeur=rangml
!call MPI_ALLREDUCE(valeur,somme,1,MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD,codeml)
!write(*,*) 'somme', somme, rangml

!debug if (rangml==0) write(*,*) 'before  broadcast', rangml,kbfunc(1:4)
call MPI_BCAST(kbfunc,dim_train,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,codeml)
!debug write(*,*) 'after broadcast', rangml,kbfunc(1:4)
call MPI_BARRIER(MPI_COMM_WORLD,codeml)
#endif


  call system_clock (count2,count_rate,count_max)
  time=real((count2-count1))/real(count_rate)
  if (debug) then
  if (rangml==0) write(*,"(' ML: MATFOR  was factorized/solved in...:  ',f16.8,' s')") time
  end if

end subroutine cholesky_threading




#ifdef PARAML

subroutine cholesky_scalapack
 use mpi
 use ml_in_ndm_module, ONLY: debug, marginal_likelihood, lambda_krr,rangml
 use temporary_data_cov, ONLY: dim_train,yfunc_train,kbfunc
 implicit none
 integer   :: istat,info,i
 real(kind=8),dimension(:,:),allocatable  :: a
 real(kind=8),dimension(:,:),allocatable  :: a_ini
 real(kind=8),dimension(:,:),allocatable  :: d_a
 real(kind=8),dimension(:,:),allocatable  :: a_temp
 real(kind=8),dimension(:,:),allocatable  :: b_temp
 real(kind=8),dimension(:,:),allocatable  :: c
 real(kind=8),dimension(:,:),allocatable  :: c_temp
 real(kind=8),dimension(:,:),allocatable  :: c_scalar
 real(kind=8),dimension(:,:),allocatable  :: z


 real(kind=8),dimension(:),allocatable  :: trace_vect
 real(kind=8)   :: val_trace,PDLATRA


 integer   ::       context, iam
 integer   ::       mycol, myrow, nb
 integer   ::       npcol, nprocs, nprow
 integer   ::       l_nrowsa,l_ncolsa
 integer   ::       l_nrowsa_ini,l_ncolsa_ini
 integer   ::       l_nrowsa_temp,l_ncolsa_temp
 integer   ::       l_nrowsb_temp,l_ncolsb_temp
 integer   ::       l_nrowsc,l_ncolsc
 integer   ::       l_nrowsc_temp,l_ncolsc_temp
 integer   ::       l_nrowsc_scalar,l_ncolsc_scalar

 integer,parameter :: descriptor_len=9
 integer   ::       desca     ( descriptor_len )
 integer   ::       desca_ini ( descriptor_len )
 integer   ::       descd_a   ( descriptor_len )
 integer   ::       descz     ( descriptor_len )
 integer   ::       descc     ( descriptor_len )
 integer   ::       descc_temp( descriptor_len )
 integer   ::       descc_scalar( descriptor_len )
 integer   ::       desca_temp( descriptor_len )
 integer   ::       descb_temp( descriptor_len )

 integer   ::       n, numroc, ITEM
 integer :: ilocal, jlocal, iglobal, jglobal
 integer :: INDXL2G, INDXG2L
 real(kind=8) :: alpha

      ! interface
      !  !
      !subroutine init_my_matrix_ml (context,iam,nprocs,n,a,d_a,nprow,npcol,myrow,mycol,desca,descd_a)
      !implicit none
      !integer :: n
      !integer :: context,iam,nprocs
      !integer :: nprow,npcol,myrow,mycol
      !integer :: desca(:), descd_a(:)
      !real(kind=8)    :: a(:,:), d_a(:,:)
      !end subroutine init_my_matrix_ml

      !end interface



 n=dim_train
 call blacs_pinfo( iam, nprocs )
 ITEM = 1000 + iam

 if (debug) then
  if  (iam == 0)    write (6,'(" nprocs = ",i4)') nprocs
 end if

 call gridsetup_ml(nprocs,nprow,npcol)
 call blacs_get( -1, 0, context )
 call BLACS_GRIDINIT( context, 'r', nprow, npcol )
 call BLACS_GRIDINFO( context, nprow, npcol, myrow, mycol )
 call blockset_ml( nb, 64, n, nprow, npcol)

 if (debug) then
   if (iam==0) write(*,*) 'ML: The block was set to ....',nb
 end if


!descriptors for 2D matrix ...
 l_nrowsa = numroc(n,nb,myrow,0,nprow)
 l_ncolsa = numroc(n,nb,mycol,0,npcol)

 l_nrowsa_temp = numroc(n,nb,myrow,0,nprow)
 l_ncolsa_temp = numroc(n,nb,mycol,0,npcol)

 l_nrowsa_ini = numroc(n,nb,myrow,0,nprow)
 l_ncolsa_ini = numroc(n,nb,mycol,0,npcol)



 l_nrowsb_temp = numroc(n,nb,myrow,0,nprow)
 l_ncolsb_temp = numroc(n,nb,mycol,0,npcol)



 call descinit( desca,      n, n, nb, nb, 0, 0, context, l_nrowsa,      info )
 call descinit( desca_ini,  n, n, nb, nb, 0, 0, context, l_nrowsa_ini,  info )
 call descinit( descd_a,    n, n, nb, nb, 0, 0, context, l_nrowsa,      info )
 call descinit( desca_temp, n, n, nb, nb, 0, 0, context, l_nrowsa_temp, info )
 call descinit( descb_temp, n, n, nb, nb, 0, 0, context, l_nrowsb_temp, info )
 call descinit( descz,      n, n, nb, nb, 0, 0, context, l_nrowsa,      info )
! descriptors for 1D matrix ...
 l_nrowsc = numroc(n,nb,myrow,0,nprow)
!debug l_ncolsc = numroc(1,nb,mycol,0,npcol)
 l_ncolsc = numroc(1,1,mycol,0,npcol)

 l_nrowsc_temp = numroc(n,nb,myrow,0,nprow)
!debug l_ncolsc_temp = numroc(1,nb,mycol,0,npcol)
 l_ncolsc_temp = numroc(1,1,mycol,0,npcol)

 call descinit( descc,      n, 1, nb, 1, 0, 0, context, l_nrowsc, info )
 call descinit( descc_temp, n, 1, nb, 1, 0, 0, context, l_nrowsc, info )

 l_nrowsc_scalar = numroc(1,nb,myrow,0,1)
 l_ncolsc_scalar = numroc(1,nb,mycol,0,1)

 call descinit( descc_scalar, 1, 1, nb, 1, 0, 0, context, l_nrowsc_scalar, info )



 if (allocated(a)) deallocate(a); allocate (a(l_nrowsa, l_ncolsa ), stat=istat)
 if (istat/=0) stop "ERR:ALLOCATE FAILS for A"

 if (allocated(a_ini)) deallocate(a_ini); allocate (a_ini(l_nrowsa_ini, l_ncolsa_ini ), stat=istat)
 if (istat/=0) stop "ERR:ALLOCATE FAILS for A_INI"

 if (allocated(d_a)) deallocate(d_a); allocate (d_a(l_nrowsa, l_ncolsa ), stat=istat)
 if (istat/=0) stop "ERR:ALLOCATE FAILS for D_A"

 if (allocated(a_temp)) deallocate(a_temp); allocate (a_temp(l_nrowsa_temp, l_ncolsa_temp ), stat=istat)
 if (istat/=0) stop "ERR:ALLOCATE FAILS for A_TEMP"

 if (allocated(b_temp)) deallocate(b_temp); allocate (b_temp(l_nrowsb_temp, l_ncolsb_temp ), stat=istat)
 if (istat/=0) stop "ERR:ALLOCATE FAILS for B_TEMP"



 if (allocated(z)) deallocate(z); allocate (z(l_nrowsa, l_ncolsa ), stat=istat)
 if (istat/=0) stop "ERR:ALLOCATE FAILS for Z"

 if (allocated(c)) deallocate(c); allocate (c(l_nrowsc,l_ncolsc), stat=istat)
 if (istat/=0) stop "ERR:ALLOCATE FAILS for C"
 if (allocated(c_temp)) deallocate(c_temp); allocate (c_temp(l_nrowsc_temp,l_ncolsc_temp), stat=istat)
 if (istat/=0) stop "ERR:ALLOCATE FAILS for C_TEMP"
 if (allocated(c_scalar)) deallocate(c_scalar); allocate (c_scalar(l_nrowsc_scalar,l_ncolsc_scalar), stat=istat)
 if (istat/=0) stop "ERR:ALLOCATE FAILS for C_SCALAR"




 call init_my_matrix_ml (context,iam,nprocs,dim_train,A,D_A,nprow,npcol,myrow,mycol,desca,descd_a)
 A_INI(:,:) =A(:,:)

 if (debug) then
   if (iam==0) write(*,*) 'ML: Matrix was initiated ....'
 end if
 Z(:,:) =0.d0
!symmetrization of the hessian matrix A = 1/2(A + A^T)
 ! Z(:,:)=A(:,:)
 !call pdtran(N,N,0.5d0,A,1,1,desca,0.5d0,Z,1,1,descz)

!construct the matrix A + \lambda_krr I
  ! indentity matixwith lambda_krr on diagonal
  do i=1,dim_train
    call PDELSET(Z,i,i,DESCZ,lambda_krr)
  end do
  ! do the sum in A <- K + \lambda_krr I
  call PDTRAN(N,N,1.d0,A,1,1,desca,1.d0,Z,1,1,descz)
  A(:,:) = Z(:,:)

 do i =1,dim_train
  call PDELSET(C,i,1,DESCC, yfunc_train(i))
 end do
 !call  PDPOSV( 'L', N, NRHS, A, 1, 1, DESCA, C, 1, 1, DESCC, INFO )
 call  PDPOSV( 'L', N, 1, A, 1, 1, DESCA, C, 1, 1, DESCC, INFO )
 ! From here in A is stored   ( K + \lambda_krr I )^{-1}

 if (marginal_likelihood) then
    ! compute the inverse matrix (K+lambda_krr I)^{-1} from the L factor of Cholesky factorization
    call PDPOTRI('L',N,A,1,1,DESCA,INFO)

    !building a_temp(n,n)=alpha(n)*alpha(n)^T + (K+lambda_krr I)^{-1} (here K+lambaI-1 is A)
    a_temp(:,:)=a(:,:)
    call PDSYR('L',N,1.d0,C,1,1,DESCC,1,A_TEMP,1,1,DESCA_TEMP)

    ! put the upper part of A_temp to zero ...
    do ilocal=1,l_nrowsa_temp
      do jlocal=1,l_ncolsa_temp
         iglobal=INDXL2G (ilocal,nb,myrow,0,nprow)
         jglobal=INDXL2G (jlocal,nb,mycol,0,npcol)
          if (iglobal>jglobal) then
           a_temp(ilocal,jlocal)=0.d0
          end if
      end do
    end do
    !multiply the diagonal with 0.5d0

    write(*,*) 'ini0'
    do i=1,n
      call PDELGET('A',' ',alpha,a_temp,i,i,desca_temp)
      call PDELSET(a_temp,i,i,desca_temp,0.5d0*alpha)
    end do
    write(*,*) 'ini1'
    z(:,:)=a_temp(:,:)
    call PDTRAN(N,N,0.5d0,Z,1,1,desca,0.5d0,a_temp,1,1,desca_temp)
!   grad_lambda_krr =

    val_trace=PDLATRA(N,a_temp,1,1,desca_temp)
    write(*,*) 'val_trace01', val_trace



    !perform the matrix multiplication ...
    ! grad for lambda_krr
    Z(:,:)=0.d0
    do i=1,dim_train
     call PDELSET(Z,i,i,DESCZ,1.d0)
    end do
    !  Z <- d (K + lamba_krr I)/d lambda_krr  = K + I = A_INI + Z
    call PDTRAN(N,N,1.d0,A_INI,1 ,1, DESCA_TEMP,1.d0,Z,1,1,descz)
    val_trace=PDLATRA(N,Z,1,1,descz)
    write(*,*) 'val_trace02', val_trace

    a_temp=z
    if (allocated(trace_vect)) deallocate(trace_vect); allocate(trace_vect(n))
    c_temp(:,:)=0.d0
    do i=1,n
     call PDGEMM('N','N',1,1,N,1.d0,a_temp,i,1,desca_temp,Z,1,i,descz,0.d0,c_temp,i,1,descc_temp)

    ilocal=INDXG2L (i,nb,myrow,0,nprow)
     write(*,*) 'i',rangml, i,c_temp(ilocal,1)
    ! call PDELGET('A',' ',trace_vect(i),c_temp,i,1,descc_temp)
    end  do

    write(*,*) 'here'

    if (rangml==0) then
       do i=1,n
         write(41,*) trace_vect(i)
       end do
    end if

    write(*,*) 'val_trace', SUM(trace_vect(1:n))



   !  do the product Z<-atemp(:,:) x (K+I) + 0*Z
    call PDSYMM('L','L',n,n,1.d0,a_temp,1,1,desca_temp,z,1,1,descz,0.d0,b_temp,1,1,descb_temp)
    val_trace=PDLATRA(N,b_temp,1,1,descb_temp)
!   grad_lambda_krr =

    write(*,*) 'val_trace2', val_trace
    stop



    stop
 end if


 if (allocated(kbfunc)) deallocate(kbfunc); allocate(kbfunc(dim_train))

  do i= 1, n
    call pdelget('A',' ',kbfunc(i),c,i,1,descc)
  enddo

  if (debug) then
    write(*,'("ML: procs, kbfunc(1:4) ",i4,4f20.7)') iam, kbfunc(1:4)
  end if


return
end subroutine cholesky_scalapack


subroutine gridsetup_ml(nproc,nprow,npcol)
!
! This subroutine factorizes the number of processors (nproc)
! into nprow and npcol,  that are the sizes of the 2d processors mesh.
!
! Written by Carlo Cavazzoni
!
      implicit none
      integer nproc,nprow,npcol
      integer sqrtnp,i

      sqrtnp = int( sqrt( dble(nproc) ) + 1 )
      do i=1,sqrtnp
        if(mod(nproc,i).eq.0) nprow = i
      end do
      npcol = nproc/nprow

      return
end subroutine gridsetup_ml
!-----------------------------------------------------------------------
subroutine blockset_ml( nb, nbuser, n, nprow, npcol)
!
!     This subroutine try to choose an optimal block size
!     for the distributd matrix.
!
!
      implicit none
      integer :: N
      integer nb, nprow, npcol, nbuser

      nb = min ( n/nprow, n/npcol )
      if(nbuser.gt.0) then
        nb = min ( nb, nbuser )
      endif
      nb = max(nb,1)

      return
end subroutine blockset_ml
!---------------------------------

subroutine init_my_matrix_ml (context,iproc,nprocs,N,A,D_A, NPROW,NPCOL,MYROW,MYCOL,DESCA,DESCD_A)
      use mpi
      use ml_in_ndm_module, ONLY: iread_ml,debug
      use temporary_data_cov, ONLY: mlocal,d_mlocal, u_local,v_local
      implicit none
      integer      :: N
      integer      :: NPROW,NPCOL,MYROW,MYCOL
      integer      :: DESCA(:),DESCD_A(:)
      REAL(kind=8) :: A(:,:), D_A(:,:)
      integer :: context,iproc,nprocs
   real(kind=8) :: AIJ, D_AIJ!, MPI_Wtime

   integer :: I,J, jptest,ilocal

 integer :: prank,qrank
 real(8) :: temps_i1,temps_i2,temps_i3,temps_f,temps_b, &
            temps_i4
 real(kind=kind(1.d0)), dimension(:), allocatable :: mlocali,d_mlocali
 integer, dimension(:), allocatable :: u_locali,v_locali
 character*6 :: namefile
 integer :: imax, imaxi,inamefile
 integer :: simat(1,1),rimat(1,1)


 call BLACS_BARRIER( context, 'A' )


     temps_i1=MPI_Wtime()

   !the case when all the matrices are written of the HDD iread_ml==1:
   if (iread_ml==1) then
          inamefile=100000+iproc
          write(namefile,'(i6)') inamefile
          namefile=TRIM(namefile)
          open(71,file=namefile//".mml",status='unknown', access='sequential',form='unformatted')
          open(72,file=namefile//".uml",status='unknown', access='sequential',form='unformatted')
          open(73,file=namefile//".vml",status='unknown', access='sequential',form='unformatted')
          open(74,file=namefile//".dml",status='unknown', access='sequential',form='unformatted')

             read(71)imax
             read(74)imax
             write (*,*) 'reading ... iproc, imax,',iproc, imax
             do ilocal=1,imax
               !if (iproc==0) write(*,*) ilocal
               read (71) mlocal (ilocal)
               read (74) d_mlocal (ilocal)
               read (72) u_local(ilocal)
               read (73) v_local(ilocal)
            end do
           close (71,status='keep')
           close (72,status='keep')
           close (73,status='keep')
           close (74,status='keep')

     call BLACS_BARRIER( context, 'A' )
   end if !end iread_ml==1

    temps_i2=MPI_Wtime()
    if(debug) then
     if (iproc == 0 ) then
      write(*,*) 'ML: Reading time..............: ', temps_i2 - temps_i1
     end if
    end if

     temps_b=0.d0
     temps_f=0.d0
! End of the reading part ...................

   do jptest=1,nprocs

      temps_i3=MPI_Wtime()

      if (iproc == (jptest-1) ) then
         prank=int(iproc/npcol)
         qrank=iproc-prank*npcol
          if  (prank /= myrow)   then
           write(*,*) 'Problems in prank, myrow  ', prank,myrow
          end if
          if (qrank /= mycol )   then
           write(*,*) 'Problems in qrank,mycol  ', qrank,mycol
          end if
         simat(1,1)=size(mlocal)
         imaxi=size(mlocal)
         ! write(*,*) 'on iproc ', iproc, 'we send imax   ', simat(1,1), i_local_cov
         call IGEBS2D(context,'All',' ',1,1,simat,1)
       end if
       !
       if ( iproc /= (jptest-1) ) then
         prank=int((jptest-1)/npcol)
         qrank=(jptest-1)-prank*npcol
         call IGEBR2D(context,'All',' ',1,1,rimat,1,prank,qrank)
         !  write(*,*) 'on iproc ', iproc, 'we recieve imax', rimat(1,1)
         imaxi=rimat(1,1)
       end if

       if (allocated(mlocali))    deallocate (mlocali)  ; allocate(mlocali(imaxi))
       if (allocated(d_mlocali))  deallocate (d_mlocali); allocate(d_mlocali(imaxi))
       if (allocated(u_locali))   deallocate (u_locali) ; allocate(u_locali(imaxi))
       if (allocated(v_locali))   deallocate (v_locali) ; allocate(v_locali(imaxi))

       call BLACS_BARRIER( context, 'A' )
       if (debug) then
        if (iproc==0) write(*,*) 'ML:  The imax was set on all procs ....'
       end if
! the second broadcast with the matrices ...

       if (iproc == (jptest-1)) then
         prank=int(iproc/npcol)
         qrank=iproc-prank*npcol
         if (debug) then
           write(*,*) 'ML: proc, size mlocal mlocali', iproc, size(mlocali), size(mlocal)
         end if
         mlocali(:)=mlocal(:)
         d_mlocali(:)=d_mlocal(:)
         u_locali(:)=u_local(:)
         v_locali(:)=v_local(:)
         call DGEBS2D(context,'All',' ',imaxi,1,mlocal ,imaxi)
         call DGEBS2D(context,'All',' ',imaxi,1,d_mlocal ,imaxi)
         call IGEBS2D(context,'All',' ',imaxi,1,u_local,imaxi)
         call IGEBS2D(context,'All',' ',imaxi,1,v_local,imaxi)
       end if
       !
       if ( iproc /= (jptest-1) ) then
        prank=int((jptest-1)/npcol)
        qrank=(jptest-1)-prank*npcol
        call DGEBR2D(context,'All',' ',imaxi,1,mlocali ,imaxi,prank,qrank)
        call DGEBR2D(context,'All',' ',imaxi,1,d_mlocali ,imaxi,prank,qrank)
        call IGEBR2D(context,'All',' ',imaxi,1,u_locali,imaxi,prank,qrank)
        call IGEBR2D(context,'All',' ',imaxi,1,v_locali,imaxi,prank,qrank)
      end if

!     call MPI_BARRIER(MPI_COMM_WORLD,codeml)

     call BLACS_BARRIER( context, 'A' )


      temps_i4=MPI_Wtime()
      temps_b=temps_b + (temps_i4-temps_i3)
          do ilocal=1,imaxi
          I= u_locali(ilocal)
          J= v_locali(ilocal)
          AIJ=mlocali(ilocal)
          D_AIJ=d_mlocali(ilocal)
          if (I > 3*N) write(*,*) 'WARNING I',I
          CALL PDELSET(A,I,J,DESCA,AIJ)
          CALL PDELSET(D_A,I,J,DESCD_A,D_AIJ)

      end do

     end do

 return
end subroutine init_my_matrix_ml

!-----------------------------------------------------------------------
!
      subroutine get_solution_ml (n,i,z,descz,eigvect)
      implicit none
      integer :: n
      integer :: descz(:)
      real(kind=8) ::  z(:,:)
      real(kind=8) :: eigvect(n)
      integer :: i
      integer :: j

      do j= 1, n
        write(*,*) j
        call pdelget('A',' ',eigvect(j),z,j,i,descz)
      enddo

      return
      end subroutine get_solution_ml


!
!-----------------------------------------------------------------------
!

subroutine  diago_scalapack(N,W)
      use gen_com_m, ONLY: im,imm,angst,at
      !use tab_imm_m, ONLY :xp
      use mpi
      implicit none
!      include 'mpif.h'
       
      logical :: leigenvectors=.true.
      integer   :: istat,info,i
      real(kind=8),dimension(:,:),allocatable  :: a
      real(kind=8),dimension(:,:),allocatable  :: d_a
      real(kind=8),dimension(:,:),allocatable  :: c
      real(kind=8),dimension(:,:),allocatable  :: z
      integer :: N
      real(kind=8) :: W(N)
      !real(kind=8),dimension(:),allocatable  :: solution
      integer,dimension(:),allocatable :: ipiv
 
      real(kind=8) :: temps1,temps2
      integer   :: root
      real(kind=8),dimension(:),allocatable    :: memory_a
 
      integer   ::       context, iam
      integer   ::       mycol, myrow, nb
      integer   ::       npcol, nprocs, nprow
      integer   ::       l_nrowsa,l_ncolsa
      integer   ::       l_nrowsc,l_ncolsc
 !     PDSYEVX interface
      integer, allocatable, dimension(:,:)    :: ra, ca
      integer                                 :: imin,imax
     integer                                  :: nvectors, nvalues
      integer   :: nns,npp, npq,lwork,liwork !,numroc
      integer, allocatable, dimension(:)      :: ifail,iclustr,iwork
      real(kind=8), allocatable, dimension(:) :: work,gap
      real(kind=8)                            :: orfac, vl, vu, abstol
      real(kind=8), external                  :: pslamch

 
      integer,parameter :: descriptor_len=9
      integer   ::       desca( descriptor_len )
      integer   ::       descd_a( descriptor_len )
      integer   ::       descz( descriptor_len )
      integer   ::       descc( descriptor_len )
 
      integer   ::       numroc, ITEM 

      integer :: neig, nn,np0,mq0,anb,sqnpc,nps,nsytrd_lwopt
      integer , external :: pjlaenv, iceil
      integer :: itest

 
      !interface
      !  !
      !subroutine init_my_matrix_ml (context,iam,nprocs,n,a,nprow,npcol,myrow,mycol,desca)
      !implicit none
      !integer :: N
      !integer :: context,iam,nprocs
      !integer :: nprow,npcol,myrow,mycol
      !integer :: desca(:)
      !real(kind=8)    :: a(:,:)
            !    integer :: n
            !    integer :: context,iam,nprocs
            !    integer :: nprow,npcol,myrow,mycol
            !    integer :: desca(:), descd_a(:)
            !    real(kind=8)    :: a(:,:), d_a(:,:)
      !end subroutine init_my_matrix_ml


      !subroutine get_solution_ml (N,i,z,descz,eigvect)
      !implicit none
      !integer :: N
      !integer :: i 
      !integer :: descz(:)
      !real(kind=8) ::  z(:,:)
      !real(kind=8) :: eigvect(N)
      !end subroutine get_solution_ml

      !end interface
 
! ---------------------------------------------------------------------------
 
      root = 0
 
! -----    Initialize the blacs.  Note: processors are counted starting at 0.
!
    call blacs_pinfo( iam, nprocs )
    ITEM = 1000 + iam
    if  (iam == 0)    write (6,'(" nprocs = ",i4)') nprocs
    call gridsetup_ml (nprocs,nprow,npcol)
    call blacs_get( -1, 0, context )
    call BLACS_GRIDINIT( context, 'r', nprow, npcol )
    call BLACS_GRIDINFO( context, nprow, npcol, myrow, mycol )
    call blockset_ml( nb, 64, n, nprow, npcol)
    if (iam==0) write(*,*) 'The block was set to ....',nb
!
! -----    Distributed matrices: get num. local rows/cols. Create description.
! -----    INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
! -----    NUMROC computes the NUMber of Rows Or Columns of a distributed
!------    matrix owned by the process indicated by IPROC.

      l_nrowsa = numroc(n,nb,myrow,0,nprow)
      l_ncolsa = numroc(n,nb,mycol,0,npcol)
!      write(*,'("1 proc, n,nb,l_nrowsa,l_ncolsa  ",5i7)') iam, n,nb,l_nrowsa,l_ncolsa

      call descinit( desca, n, n, nb, nb, 0, 0, context, l_nrowsa, info )
      call descinit( descz, n, n, nb, nb, 0, 0, context, l_nrowsa, info )
 
      l_nrowsc = numroc(n,nb,myrow,0,nprow)
      l_ncolsc = numroc(1,nb,mycol,0,npcol)
      call descinit( descc, n, 1, nb, nb, 0, 0, context, l_nrowsc, info )
!
! -----   Allocate LHS, RHS, pivot, and solution -----
!
      allocate (a(l_nrowsa, l_ncolsa ), stat=istat)
      if (istat/=0) stop "ERR:ALLOCATE FAILS for A"
 
      allocate (z(l_nrowsa, l_ncolsa ), stat=istat)
      if (istat/=0) stop "ERR:ALLOCATE FAILS for Z"

 
      allocate (c(l_nrowsc,l_ncolsc), stat=istat)
      if (istat/=0) stop "ERR:ALLOCATE FAILS for C"
 
      allocate (ipiv (n), stat=istat)
      if (istat/=0) stop "ERR:ALLOCATE FAILS for IPIV"
 
 
      allocate (memory_a (nprocs), stat=istat)
      if (istat/=0) stop "ERR:ALLOCATE FAILS for memory_a"

!     PDSYEVX interface
      allocate (gap (nprocs), stat=istat)
      if (istat/=0) stop "ERR:ALLOCATE FAILS for gap"
      allocate (ifail (n), stat=istat)
      if (istat/=0) stop "ERR:ALLOCATE FAILS for ifail"
      allocate (iclustr (2*nprow*npcol), stat=istat)
      if (istat/=0) stop "ERR:ALLOCATE FAILS for iclustr"
      allocate (ra (nb,nb), stat=istat)
      if (istat/=0) stop "ERR:ALLOCATE FAILS for ra"
      allocate (ca (nb,nb), stat=istat)
      if (istat/=0) stop "ERR:ALLOCATE FAILS for ca"


 
!
! -----    Initialize LHS and RHS
!
      !call  read_my_matrix (context,iam,nprocs,n,a,nprow,npcol,myrow,mycol,desca)
      !call init_my_matrix_ml (context,iam,nprocs,N,A,nprow,npcol,myrow,mycol,desca)
      

      call init_my_matrix_ml (context,iam,nprocs,N,A,D_A,nprow,npcol,myrow,mycol,desca,descd_a)
      !symmetrization of the hessian matrix A = 1/2(A + A^T)
      Z(:,:) = A(:,:)
      call pdtran(N,N,0.5d0,A,1,1,desca,0.5d0,Z,1,1,descz)
      A(:,:) = Z(:,:)
!
! -----    Show how arrays distributed
!
 


 !     PDSYEVX interface

 vl    = -15.   ! eig min. if  'V' instead 'A'
 vu    =  15.   ! eig max  if 'V' instead 'A'
 imin    = -13  ! Ind  min.if  'I'  'A'
 imax    =  13  ! Ind max. if 'I instead 'A'
 abstol  =  -1. ! if  abstol=pslamch(context,'U') but dob't TOUCH
 abstol=pslamch(context,'U') ! If  JOBZ='V',  setting  ABSTOL to PDLAMCH( CONTEXT, 'U') yields
                             ! the most orthogonal eigenvectors 
 orfac   =   0.   ! Ortho

        temps1=MPI_Wtime()

 	  nns = MAX( n, nb, 2 )
 	  npp = numroc(nns, nb, 0, 0, nprow)
 	  npq = max(n, nprow*npcol+1, 4) 
      if (leigenvectors) then
          neig=n   !number of requested eigenvectors
          nn = max(n, nb, 2)
          np0 = numroc(nn, nb, 0, 0, nprow)
          mq0 = numroc(max(neig, nb, 2), nb, 0, 0, npcol)
          lwork= 5*n + max(5*nn, np0*mq0 + 2*NB*NB) + iceil(neig, nprow*npcol)*nn
           
          anb = pjlaenv(desca, 3, 'pdsyttrd', 'L', 0, 0, 0, 0)
          sqnpc = int(sqrt(dble(NPROW * NPCOL)))
          nps = max(numroc(n, 1, 0, 0, sqnpc), 2*anb)
          nsytrd_lwopt = n + 2*(anb+1)*(4*nps+2) + (nps + 3)*nps

          lwork= max(2*lwork, 5*n + nsytrd_lwopt)
        !  write(*,*) 'setup', lwork, 7*n + nsytrd_lwopt
      else    
 	  lwork=max(nb*(npp+1),5*nns)+5*n
      end if

 	  liwork=6*npq
          allocate (iwork (liwork), stat=istat)
          if (istat/=0) stop "ERR:ALLOCATE FAILS for iwork"
          allocate (work (lwork), stat=istat)
          if (istat/=0) stop "ERR:ALLOCATE FAILS for work"
  
   
        !V/N   = eigenvalue and eigenvetore/eigenvalues
        !A/V/I = all eignvalues / aleigenvalues in the interval [vl,vu] / eigenvalues index from  imin to imax
        !U/L = upper / lower triangular 
        ! n  = matrix dimension
        ! A  = the matrix
        ! 1,1 = the beginning of the matrix
        ! desca = descriptor
        ! vl,vu,imin,imax - described above
        ! abstol - if jobz (the first argument is set to V
        ! nvalues (output) - total number of eigenvalues 
        ! nvectore (output) - total number of eigenvectors which fill the matrix Z    



   if (leigenvectors) then
       if (iam==0) print*,'Eigenvalues and eigenvectors are computed'
        call PDSYEVX('V','A','L',n,A,1,1,desca,vl,vu,imin,imax,  &
                abstol,nvalues,nvectors,W,orfac,Z,1,1,descz,work, &
                lwork,iwork,liwork,ifail,iclustr,gap,info )
      
        if (iam==0) then 
         print*,'The number of eigenvalues  computed', nvalues
         print*,'The number of eigenvectors computed', nvectors
         itest=0
         do i=1,n
          if (ifail(i)/=0) then 
           itest=1
          end if 
         end do

         if (itest==1) then
          print*,'WARNING .... some eigenvalues failed to converge'
         end if  
        end if 
!         call get_solution (N,3*504,Z,DESCZ,eigvect)
!         do i= 1, n
!           call pdelget('A',' ',eigvect(i),z,i,3*504,descz)
!         enddo




!        if (iam==0) then
!         write(*,*)  eigvect (:)
!        end if 

        
    else 
        call PDSYEVX('N','A','L',n,A,1,1,desca,vl,vu,imin,imax,  &
                abstol,nvalues,nvectors,W,orfac,Z,1,1,descz,work, &
                lwork,iwork,liwork,ifail,iclustr,gap,info )
    end if

    if (info/=0) then
         write(*,*) 'Diago not OK'
         write(*,*) 'info',info
    end if 


  	temps2=MPI_Wtime()

   if (iam==0) then
    print*,'That diago time ',temps2-temps1,' seconds'
  end if
 if ( iam == 0 ) then
   do i = 1, n
!     print '(25X,"W(",I5,") : ",1PE12.4)', i, W(i)
!      write(10,'(i6,1D18.10)') i, W(i)
     write(10,*)  i,w(i)

   end do
 end if

    deallocate (a,ipiv,memory_a,Z)
!
! -----    Exit BLACS cleanly -----
!
      call blacs_gridexit( context )
!      call blacs_exit( 0 )
 
end subroutine  diago_scalapack

!endif of PARAPH
#endif
end module 
