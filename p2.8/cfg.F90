MODULE cfg_module


CONTAINS

  SUBROUTINE WriteCfg(xp, ityp, im, at, out, mask, nAux_int, aux_int, nAux_real, aux_real, aux_title)
    ! permet d'ecrire les configurations au format cfg directement compris par atomeye.
    ! 
    ! Parametres d'entree obligatoires:
    ! xp(1:3,:) : coordonnees cartesiennes des atomes 
    ! ityp(:) : type des atomes
    ! im : nombre d'atomes
    ! at(1:3,i) : vecteur de periodicite i
    ! out : numero de l'unite connecte au fichier de sortie
    !
    ! Parametres d'entree optionnels:
    ! mask(:) : .true.  => l'atome correspondant est inclu dans le fichier cfg
    !           .false. => ------------------------- ignore
    ! nAux_int : nombre de proprietes auxiliaires au format INTEGER
    ! aux_int(1:nAux_int,:) : proprietes auxiliaires au format INTEGER
    ! nAux_real : nombre de proprietes auxiliaires au format REAL
    ! aux_real(1:nAux_real,:) : proprietes auxiliaires au format REAL
    ! aux_title(1:nAux_int+nAux_real) : nom des proprietes correspondantes
  USE temp_com,only:imm ! A EFFACER
  USE T_kind_param_m, ONLY:  double
    USE gen_com_m, ONLY: umass
    USE var_pot, ONLY:ntyp,cm,ty
    USE Mat_utils_mod
    IMPLICIT NONE
    REAL(double),                 intent(in) :: xp(3,imm)
    INTEGER,                      intent(in) :: ityp(imm)
    INTEGER,                      intent(in) :: im
    REAL(double), dimension(3,3), intent(in) :: at
    INTEGER,                      intent(in) :: out

    LOGICAL, dimension(:),           intent(in), optional :: mask
    INTEGER,                         intent(in), optional :: nAux_int, nAux_real
    INTEGER, dimension(:,:),         intent(in), optional :: aux_int
    REAL(double), dimension(:,:),    intent(in), optional :: aux_real
    CHARACTER(len=20), dimension(:), intent(in), optional :: aux_title

    INTEGER :: n, i, ic, j
    CHARACTER(len=50) :: out_format
    LOGICAL, dimension(:), allocatable :: local_mask
    LOGICAL :: test_aux_int, test_aux_real
    REAL(double), dimension(3,3) :: inv_at

    !APARA
    IF (Allocated(local_mask)) Deallocate(local_mask)
    Allocate(local_mask(1:im))

    IF (Present(mask)) THEN
       local_mask(1:im)=mask(1:im)
    ELSE
       local_mask = .TRUE.
    END IF

    IF ( Present(nAux_int) ) THEN
       IF ( (.NOT.Present(aux_int)) .OR. (.NOT.Present(aux_title)) ) THEN
          WRITE(0,'(a)') 'All auxiliary properties have to be defined'
          STOP '< WriteCfg >'
       END IF
       IF (Size(aux_int,1).LT.nAux_int)  THEN
          WRITE(0,'(a)') 'Problem with size of auxiliary properties'
          STOP '< WriteCfg >'
       END IF
       test_aux_int=.TRUE.
    ELSE
       test_aux_int=.FALSE.
    END IF
    IF ( Present(nAux_real) ) THEN
       IF ( (.NOT.Present(aux_real)) .OR. (.NOT.Present(aux_title)) ) THEN
          WRITE(0,'(a)') 'All auxiliary properties have to be defined'
          STOP '< WriteCfg >'
       END IF
       IF (Size(aux_real,1).LT.nAux_real)  THEN
          WRITE(0,'(a)') 'Problem with size of auxiliary properties'
          STOP '< WriteCfg >'
       END IF
       test_aux_real=.TRUE.
    ELSE
       test_aux_real=.FALSE.
    END IF

    Call MatInv(at,inv_at)

    write(out,'(a,i0)')'Number of particles = ', Count(local_mask(1:im))
    write(out,'(a)')'A = 1.000 Angstrom (basic length-scale)'
    do j=1,3
       write(out,'(a,i0)') '# Unit cell vector #', j    
       do ic=1,3
          write(out,'(A,I1,A,I1,A,g16.8,A)')'H0(',j,',',ic,') = ',1e8*at(ic,j),' A'
       end do
    end do
    write(out,'(A)')'.NO_VELOCITY.'
    IF (test_aux_int.AND.test_aux_real) THEN
       write(out,'(A,I0)')'entry_count = ', 3+nAux_int+nAux_real
       DO i=1, nAux_int+nAux_real
          WRITE(out,'(a,i0,2a)') 'auxiliary[',i-1,'] = ', aux_title(i)
       END DO
       WRITE(out_format,'(a,3(i0,a))') &
            '(',3,'(g24.16,1x),',nAux_int,'(1x,i0),',nAux_real,'(1x,g24.16))'
    ELSEIF (test_aux_int) THEN
       write(out,'(A,I0)')'entry_count = ', 3+nAux_int
       DO i=1, nAux_int
          WRITE(out,'(a,i0,2a)') 'auxiliary[',i-1,'] = ', aux_title(i)
       END DO
       WRITE(out_format,'(a,2(i0,a))') &
            '(',3,'(g24.16,1x),',nAux_int,'(1x,i0)),'
    ELSEIF (test_aux_real) THEN
       write(out,'(A,I0)')'entry_count = ', 3+nAux_real
       DO i=1, nAux_real
          WRITE(out,'(a,i0,2a)') 'auxiliary[',i-1,'] = ', aux_title(i)
       END DO
       WRITE(out_format,'(a,2(i0,a))') &
            '(',3,'(g24.16,1x),',nAux_real,'(1x,g24.16))'
    ELSE
       write(out,'(A,I0)')'entry_count = ', 3
       WRITE(out_format,'(a,i0,a)') '(',3,'(g24.16,1x))'
    ENDIF

    DO n=1, ntyp       ! Loop on all types
       IF ( Count( iTyp(1:im).EQ.n ).LE.0) Cycle
       IF (test_aux_int.AND.test_aux_real) THEN
          WRITE(out,'(f0.3)') cm(n)/umass        ! Mass (g/mol)
          WRITE(out,'(a)') ty(n)                 ! Atom type
          do i=1, im
             IF ( (ityp(i).EQ.n).AND.(local_mask(i)) ) &
                  write(out,out_format) MatMul(inv_at(:,:),xp(:,i)), aux_int(1:nAux_int,i), &
                  aux_real(1:nAux_real,i)
          end do
       ELSEIF (test_aux_int) THEN
          WRITE(out,'(f0.3)') cm(n)/umass        ! Mass (g/mol)
          WRITE(out,'(a)') ty(n)                 ! Atom type
          do i=1, im
             IF ( (ityp(i).EQ.n).AND.(local_mask(i)) ) &
                  write(out,out_format) MatMul(inv_at(:,:),xp(:,i)), aux_int(1:nAux_int,i)
          end do
       ELSEIF (test_aux_real) THEN
          WRITE(out,'(f0.3)') cm(n)/umass        ! Mass (g/mol)
          WRITE(out,'(a)') ty(n)                 ! Atom type
          do i=1, im
             IF ( (ityp(i).EQ.n).AND.(local_mask(i)) ) &
                  write(out,out_format) MatMul(inv_at(:,:),xp(:,i)), aux_real(1:nAux_real,i)
          end do
       ELSE
          WRITE(out,'(f0.3)') cm(n)/umass        ! Mass (g/mol)
          WRITE(out,'(a)') ty(n)                 ! Atom type
          do i=1, im
             IF ( (ityp(i).EQ.n).AND.(local_mask(i)) ) &
                  write(out,out_format) MatMul(inv_at(:,:),xp(:,i))
          end do
       ENDIF
    END DO

    Deallocate(local_mask)

  END SUBROUTINE WriteCfg

END MODULE cfg_module
