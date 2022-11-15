
  subroutine diago_threading_real()
    USE T_kind_param_m, ONLY:  double
    use gen_com_m
    use var_pot
    use tab_imm_m
    use phondy_in_ndm_module
    use thermo_module

 implicit none

  !diago driver
  real(double), dimension(:), allocatable :: work
  character*1  :: uplo
  integer :: info, lwork
  !diago driver
  integer :: count1,count2,count_rate,count_max
  real(double) :: time
  integer :: i

  call system_clock (count1,count_rate,count_max)

 uplo='u'
 info=0
 lwork=max(1,3*nmat-1)+nmat
 allocate(work(lwork))
 write(*,*) 'PHONDY: Starting diago ....'
 if (leigenvectors) then
    call DSYEV( 'V', 'L', nmat, matfor, nmat, w, work, lwork, info )
    if (info/=0) then
       write(*,*) 'WARNING there are problems in the diago with eigenvalues and eigenvectors'
    end if
   else
    call DSYEV( 'N', 'L', nmat, matfor, nmat, w, work, lwork, info )
 end if

 write(*,*) 'PHONDY: ... the diago completed'
!debug do i=1,size(W, DIM=1)
!debug    write(774,*) i, sign(1.d0,W(i))*dsqrt(dabs(W(i)))*unit_nu*1.D3
!debug  end do

call system_clock (count2,count_rate,count_max)
 time=real((count2-count1))/real(count_rate)
 write(*,"(' PHONDY: DIAGO was  ................:  ',f16.8,' s')") time
#if (PARAPH)
 write(*,*) 'NO LDOS IN THIS COMPILATION'
#else
 if (lldos) call  write_ldos_threading (nmat,w,matfor)
#endif

end subroutine diago_threading_real



subroutine diago_threading_complex()
  USE T_kind_param_m, ONLY:  double
  use gen_com_m
  use var_pot
  use tab_imm_m
  use phondy_in_ndm_module
  use thermo_module

  implicit none

  !diago driver
  complex(kind(1.d0)), dimension(:), allocatable :: work
  real(double), dimension(:), allocatable :: rwork
  character*1  :: uplo
  integer :: info, lwork
  !diago driver
  integer :: count1,count2,count_rate,count_max
  real(double) :: time
  integer :: i

  call system_clock (count1,count_rate,count_max)

 uplo='u'
 info=0
 lwork=max(1,3*nmat-1)+nmat
 if (allocated(work)) deallocate(work) ; allocate(work(lwork))
 if (allocated(rwork)) deallocate(rwork) ; allocate(rwork(3*nmat-2))

!if (debug)

!end if
if (debug_ph) then
  if (rangph==0) write(6,*) 'PHONDY: Starting diago .... for matrix nmat', nmat
  if (rangph==0) write(6,*) 'nmat, lwork, size(rwork(:),DIM=1)',   nmat, lwork, size(rwork(:),DIM=1)
end if
if (leigenvectors) then
    call ZHEEV( 'V', 'L', nmat, dynmat, nmat, w, work, lwork, rwork, info )
    if (info/=0) then
       write(*,*) 'WARNING there are problems in the diago with eigenvalues and eigenvectors'
    end if
   else
    call ZHEEV( 'N', 'L', nmat, dynmat, nmat, w, work, lwork, rwork, info )
 end if

if (debug_ph)  write(*,*) 'PHONDY: ... the diago completed'


 call system_clock (count2,count_rate,count_max)
 time=real((count2-count1))/real(count_rate)
if (debug_ph)  write(*,"(' PHONDY: DIAGO was  ................:  ',f16.8,' s')") time

end subroutine diago_threading_complex




#if(PARAPH)



subroutine  diago_scalapack(N,W)
      use gen_com_m, ONLY: im,imm,angst,at, rangph
      !use mod_mpi_phondy
      !use tab_imm_m, ONLY :xp
      use phondy_in_ndm_module, ONLY: leigenvectors,lldos,lmodes,unit_nu,imodes,nmodes, l_thermo_atoms
      use mpi
      use thermo_module
      use diago_scalapack_real, only : Z, DESCZ
      implicit none
!      include 'mpif.h'

      integer   :: istat,info,i
      real(kind=8),dimension(:,:),allocatable  :: a
      real(kind=8),dimension(:,:),allocatable  :: c
      !debug real(kind=8),dimension(:,:),allocatable  :: z
      integer :: N
      real(kind=8) :: W(N)
      !real(kind=8),dimension(:),allocatable  :: solution
      integer,dimension(:),allocatable :: ipiv

      real(kind=8) :: memory_local, memory_sum,temps1,temps2,memory_min,memory_max
      integer   :: root, ierr
      real(kind=8),dimension(:),allocatable    :: memory_a
      logical :: PRT_MEM

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
      real(kind=8), external                  :: pslamch, pdlamch


      integer,parameter :: descriptor_len=9
      integer   ::       desca( descriptor_len )
      integer   ::       descc( descriptor_len )

      integer   ::       numroc, ITEM,  NH, NH1

      integer :: neig, nn,np0,mq0,anb,sqnpc,nps,nsytrd_lwopt
      integer , external :: pjlaenv, iceil
      integer :: itest

!------------------local dos related-----------------------------
!for future developement (Q-points)



! --------------------------------------------------------------------------
      interface
        !
        subroutine init_my_matrix (context,iam,nprocs,n,a,nprow,npcol,myrow,mycol,desca)
          implicit none
          integer :: N
          integer :: context,iam,nprocs
          integer :: nprow,npcol,myrow,mycol
          integer :: desca(:)
          real(kind=8)    :: a(:,:)
        end subroutine init_my_matrix

        !
!        end subroutine write_thermo_atoms


        subroutine get_solution (N,i,z,descz,eigvect)
          implicit none
          integer :: N
          integer :: i
          integer :: descz(:)
          real(kind=8) ::  z(:,:)
          real(kind=8) :: eigvect(N)
        end subroutine get_solution

      end interface

! ---------------------------------------------------------------------------

      if (allocated(DESCZ)) deallocate(DESCZ) ;    allocate(DESCZ(descriptor_len))
!      N=128000
    !  N=431*3
      PRT_MEM = .TRUE.


!
!     Trace intermediate steps only for small values of N

      root = 0

!     N must be EVEN

      NH    = N/2
      NH1   = NH + 1


!
! -----    Initialize the blacs.  Note: processors are counted starting at 0.
!
      call blacs_pinfo( iam, nprocs )


     if (iam==0) then
      IF (2*NH .NE. N) THEN
        WRITE(6,2001) N
 2001   FORMAT("N must be EVEN, got N = ",i9)
        print*, 'But this was in the past. Now we can compute even for EVEN. Sic!'
       ! STOP
      ENDIF
     end if



! ---------------------------------------------------------------------------

     ITEM = 1000 + iam


        if  (iam == 0)    write (6,'(" nprocs = ",i4)') nprocs
! -----    Set the dimension of the 2d processors grid.
!
!
! This subroutine factorizes the number of processors (nproc)
! into nprow and npcol,  that are the sizes of the 2d processors mesh.


      !call MPI_Dims_Create(nprocs,2,dims,ierr)
      !nprow=dims(1)
      !npcol=dims(2)
       call gridsetup(nprocs,nprow,npcol)
!
! -----    Initialize a single blacs context.  Determine which processor I
!          am in the 2D process or grid.
!
      !call blacs_get( -1, 0, context )
      call blacs_get( -1, 0, context )
      call BLACS_GRIDINIT( context, 'r', nprow, npcol )
      call BLACS_GRIDINFO( context, nprow, npcol, myrow, mycol )

!
! -----    Calculate the blocking factor for the matrix.
!
!     This subroutine try to choose an optimal block size
!     for the distributd matrix.

      call blockset( nb, 64, n, nprow, npcol)
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

! -----    Reassemble memory_a on all processors. Add and Print on 0 -----
!

       memory_local = l_nrowsa*l_ncolsa*8

       call MPI_Gather(memory_local, 1, MPI_REAL8,                      &
     &                memory_a, 1, MPI_REAL8,   root,                   &
     &                MPI_COMM_WORLD, ierr)


      IF (IAM.EQ.0) THEN
        IF ( PRT_MEM ) THEN
          ITEM = 6500
          memory_sum = 0
          write (6,550) ITEM,N
  550     FORMAT(I5," N = ",i12 /" Distributed memory for A: ")
          DO I=0,nprocs-1
            ITEM = 6500 + I
!            WRITE (6,551) ITEM,I,memory_a(I+1)
            memory_sum = memory_sum + memory_a(I+1)
!  551       FORMAT(I5,"       processor = ",I5," ",f14.1)
          ENDDO

          memory_min=MINVAL(memory_a(1:nprocs))
          memory_max=MAXVAL(memory_a(1:nprocs))

          ITEM = ITEM + 1
          WRITE (6,552) ITEM,memory_sum/1.0d+6, float(N)*float(N)*8/1.0d6
  552     FORMAT(I5," TOTAL MEMORY (Mo) for A = ",f14.1," ==  N*N*8 = ",f14.1)
          write(6,'("The distributed maximum memory on one proc......:",f14.1)')  memory_max/1.0d+6
          write(6,'("The distributed minimum memory on one proc......:",f14.1)')  memory_min/1.0d+6
        ENDIF
        call flush(6_4)

      ENDIF

!
! -----    Initialize LHS and RHS
!
      !call  read_my_matrix (context,iam,nprocs,n,a,nprow,npcol,myrow,mycol,desca)

      call init_my_matrix (context,iam,nprocs,N,A,nprow,npcol,myrow,mycol,desca)
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
 !old abstol=pslamch(context,'U') ! If  JOBZ='V',  setting  ABSTOL to PDLAMCH( CONTEXT, 'U') yields
 abstol=pdlamch(context,'U') ! If  JOBZ='V',  setting  ABSTOL to PDLAMCH( CONTEXT, 'U') yields
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



   if (allocated(work)) deallocate(work)
   if (allocated(iwork)) deallocate(iwork)
   !deallocate (a,c,z,ipiv,memory_a,gap,ifail,iclustr,ra,ca)
   deallocate (a,c,ipiv,memory_a,gap,ifail,iclustr,ra,ca)


  !debug call blacs_gridexit( context )
  !debug call blacs_exit( 0 )

end subroutine  diago_scalapack

subroutine gridsetup(nproc,nprow,npcol)
!
! This subroutine factorizes the number of processors (nproc)
! into nprow and npcol,  that are the sizes of the 2d processors mesh.
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
end subroutine gridsetup
!-----------------------------------------------------------------------
subroutine blockset( nb, nbuser, n, nprow, npcol)
!
!     This subroutine try to choose an optimal block size
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
end subroutine blockset
!---------------------------------

subroutine init_my_matrix (context,iproc,nprocs,N,A,NPROW,NPCOL,MYROW,MYCOL,DESCA)
      use mpi
      use mod_mpi_phondy
      use gen_com_m, ONLY : rangph
      use phondy_in_ndm_module, ONLY: codeph,iread,imax,mlocal,u_local,v_local
      implicit none
      integer :: N
      INTEGER      :: NPROW,NPCOL,MYROW,MYCOL
      INTEGER      :: DESCA(:)
      REAL(kind=8) :: A(:,:)
      integer :: context,iproc,nprocs
   real(kind=8) :: AIJ!, MPI_Wtime

   integer :: I,J, jptest,ilocal

 integer :: prank,qrank
 real(8) :: temps_i1,temps_i2,temps_i3,temps_f,temps_b, &
            temps_i4
 real(8), dimension(:), allocatable :: mlocali
 integer, dimension(:), allocatable :: u_locali,v_locali
 character*6 :: namefile
 integer :: imaxi,inamefile
 integer :: simat(1,1),rimat(1,1)


 call BLACS_BARRIER( context, 'A' )

     temps_i1=MPI_Wtime()

   !the case when all the matrices are written of the HDD iread==1:
   if (iread==1) then
          inamefile=100000+iproc
          write(namefile,'(i6)') inamefile
          namefile=TRIM(namefile)
          open(71,file=namefile//".m",status='unknown', access='sequential',form='unformatted')
          open(72,file=namefile//".u",status='unknown', access='sequential',form='unformatted')
          open(73,file=namefile//".v",status='unknown', access='sequential',form='unformatted')

             read(71)imax
           !  allocate(mlocal(imax),u_local(imax),v_local(imax))
            !debug    write (*,*) 'reading ... iproc, imax,',iproc, imax
             do ilocal=1,imax
               !if (iproc==0) write(*,*) ilocal
               read (71) mlocal (ilocal)
               read (72) u_local(ilocal)
               read (73) v_local(ilocal)
            end do
           close (71,status='keep')
           close (72,status='keep')
           close (73,status='keep')

     call BLACS_BARRIER( context, 'A' )
   end if !end iread==1

    temps_i2=MPI_Wtime()

    if (iproc == 0 ) then
      write(*,*) 'Reading time..............: ', temps_i2 - temps_i1
    end if

     temps_b=0.d0
     temps_f=0.d0
! End of the reading part ...................
!debug   write(*,*) 'test iproc', iproc, rangph, nprocs, imax
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
         simat(1,1)=imax
         imaxi=imax
!debug          write(*,*) 'on iproc ', iproc, 'we send imax   ', imax
         call IGEBS2D(context,'All',' ',1,1,simat,1)
       end if
       !
       if ( iproc /= (jptest-1) ) then
         prank=int((jptest-1)/npcol)
         qrank=(jptest-1)-prank*npcol
         call IGEBR2D(context,'All',' ',1,1,rimat,1,prank,qrank)
!debug          write(*,*) 'on iproc ', iproc, 'we recieve imax', rimat(1,1)
         imaxi=rimat(1,1)
       end if

       if (allocated(mlocali))  deallocate (mlocali)
       if (allocated(u_locali)) deallocate (u_locali)
       if (allocated(v_locali)) deallocate (v_locali)
       allocate (mlocali(imaxi),u_locali(imaxi),v_locali(imaxi))

       call BLACS_BARRIER( context, 'A' )
!debug        if (iproc==0) write(*,*) 'The imax was set on all procs ....'

! the second broadcast with the matrices ...

       if (iproc == (jptest-1)) then
         prank=int(iproc/npcol)
         qrank=iproc-prank*npcol
         mlocali(:)=mlocal(:)
         u_locali(:)=u_local(:)
         v_locali(:)=v_local(:)
         call DGEBS2D(context,'All',' ',imaxi,1,mlocal ,imaxi)
         call IGEBS2D(context,'All',' ',imaxi,1,u_local,imaxi)
         call IGEBS2D(context,'All',' ',imaxi,1,v_local,imaxi)
       end if
       !
       if ( iproc /= (jptest-1) ) then
        prank=int((jptest-1)/npcol)
        qrank=(jptest-1)-prank*npcol
        call DGEBR2D(context,'All',' ',imaxi,1,mlocali ,imaxi,prank,qrank)
        call IGEBR2D(context,'All',' ',imaxi,1,u_locali,imaxi,prank,qrank)
        call IGEBR2D(context,'All',' ',imaxi,1,v_locali,imaxi,prank,qrank)
      end if
     call MPI_BARRIER(MPI_COMM_WORLD,codeph)

     call BLACS_BARRIER( context, 'A' )
      temps_i4=MPI_Wtime()
      temps_b=temps_b + (temps_i4-temps_i3)

!debug write(*,*) 'imaxi iproc', iproc,imaxi
          do ilocal=1,imaxi

          I= u_locali(ilocal)
          J= v_locali(ilocal)
          AIJ=mlocali(ilocal)
          if (I > 3*N) write(*,*) 'WARNING I',I
          CALL PDELSET(A,I,J,DESCA,AIJ)
      end do

     end do

      return
      end

!-----------------------------------------------------------------------
!
      subroutine get_solution (n,i,z,descz,eigvect)
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
      end

!
!-----------------------------------------------------------------------
!

!endif of PARAPH
#endif



      subroutine timestamp( date_time, day,hour,minute,second,millisec)
      implicit none
      character*23 :: date_time
      integer      :: day, hour, minute, second, millisec
      integer      :: elements(8)
      character*3  :: months(12)

      data months /'Jan','Feb','Mar','Apr','May','Jun',                 &
     &             'Juk','Aug','Sep','Oct','Nov','Dec' /
      call date_and_time ( VALUES=elements )

      if ( elements(1) .ne. -HUGE(0) ) then
      century: if ( elements(1) .lt. 2000 ) then
                 elements(1) = elements(1) - 1900
               else
                 elements(1) = elements(1) - 2000
               endif century
               write(date_time,1010) elements(3),                       &
     &                       months ( elements(2) ),                    &
     &                                elements(1)  , elements(5),       &
     &                 elements(6)  , elements(7)  , elements(8)

 1010          format(i2.2,1x,a3,1x,i2.2,1x,                            &
     &                i2.2,':',i2.2,':',i2.2,'.',i4.4 )

               day      = elements(3)
               hour     = elements(5)
               minute   = elements(6)
               second   = elements(7)
               millisec = elements(8)

      else
               date_time=' '

               day      = 0
               hour     = 0
               minute   = 0
               second   = 0
               millisec = 0

      endif

      end subroutine timestamp
