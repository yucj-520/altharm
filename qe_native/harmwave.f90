!
! Copyright (C) 2026 the altharm project
! This file is distributed under the terms of the GNU General Public
! License. See the file `License' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM harmwave
  !-----------------------------------------------------------------------
  !
  ! Reads a finished 'bands' calculation (collinear LSDA, nosym=.true.,
  ! noinv=.true., an explicit k-point list forming a circle in the k_x k_y
  ! plane, read from the k-points rather than assumed) and computes,
  ! DIRECTLY from Quantum ESPRESSO's own
  ! in-memory eigenvalues and its own symmetry-finding routine, two things:
  !
  !   1. The magnetic character (chi = +1 preserves the two magnetic
  !      sublattices, chi = -1 exchanges them) of every operation of the
  !      NON-magnetic parent structure, found by QE's own find_sym() called
  !      on a copy of the structure where the two magnetic sublattices are
  !      collapsed back to one chemical species (their true element, read
  !      from the pseudopotential's own psd field -- NOT from ityp, since
  !      this group's convention is to give the two sublattices distinct
  !      ATOMIC_SPECIES entries so LSDA can start them with opposite
  !      starting_magnetization, and find_sym() would otherwise never let a
  !      symmetry operation map one onto the other).
  !
  !   2. A raw least-squares fit of the spin splitting
  !      dE(theta) = E_up(theta) - E_down(theta) onto cos(l*theta)/sin(l*theta)
  !      up to l = lmax, computed by LAPACK directly on the Rydberg
  !      eigenvalues QE holds in memory (wvfct's et array) -- no XML
  !      round-trip.
  !
  ! Both outputs are written to a plain text file for cross-validation
  ! against the independent Python implementation in the altharm package
  ! (altharm.symmetry, altharm.project): the (rotation, chi) list against
  ! altharm.symmetry.spin_operations(), and the fitted (a_l, b_l) against
  ! altharm.project.project(). Which harmonics that chi-list makes ALLOWED,
  ! and the resulting noise-floor separation, are deliberately left to that
  ! already-tested code rather than re-derived here -- see PLAN.md P9.
  !
  USE kinds,           ONLY : DP
  USE constants,       ONLY : RYTOEV, PI
  USE io_files,        ONLY : prefix, tmp_dir
  USE io_global,       ONLY : ionode, ionode_id, stdout
  USE mp,              ONLY : mp_bcast
  USE mp_global,       ONLY : mp_startup
  USE mp_images,       ONLY : intra_image_comm
  USE environment,     ONLY : environment_start, environment_end
  USE klist,           ONLY : nkstot, xk
  USE wvfct,           ONLY : nbnd, et
  USE lsda_mod,        ONLY : nspin
  USE noncollin_module,ONLY : noncolin
  USE ions_base,       ONLY : nat, ityp, tau, nsp
  USE uspp_param,      ONLY : upf
  USE cell_base,       ONLY : at, bg, alat
  USE symm_base,       ONLY : nsym, s, sr, irt, find_sym
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER, PARAMETER :: max_ntyp = 20, max_ops = 200
  !
  CHARACTER(LEN=256) :: outdir, report_file
  INTEGER :: lmax, band_index
  REAL(DP) :: magmom(max_ntyp), moment_tol
  !
  NAMELIST / inputpp / prefix, outdir, lmax, band_index, magmom, &
                        moment_tol, report_file
  !
  INTEGER :: ios, it, ia, isym, ik, n, l, j, ncoef, nkept, lwork, info
  INTEGER :: chem_ityp(1000)          ! collapsed chemical type per atom
  INTEGER :: n_chem                   ! number of distinct chemical species found
  CHARACTER(LEN=2) :: chem_label(max_ntyp)
  REAL(DP) :: dummy_mloc(3, 1000)
  REAL(DP), ALLOCATABLE :: op_cart(:,:,:)   ! (3,3,nkept) kept Cartesian rotations
  INTEGER,  ALLOCATABLE :: op_chi(:)        ! (nkept)
  LOGICAL :: ok_plus, ok_minus, is_dup
  INTEGER :: ja
  REAL(DP) :: si, sj
  !
  REAL(DP), ALLOCATABLE :: theta(:), dE(:), design(:,:), coef(:), work(:)
  REAL(DP) :: kz_ref, s_sv(200)
  REAL(DP) :: nplane(3), e1(3), e2(3), c_off, v1(3), v2(3), vn(3), nrm
  INTEGER  :: ik2, ik3
  INTEGER :: iunit
  !
#if defined(__MPI)
  CALL mp_startup( )
#endif
  CALL environment_start( 'HARMWAVE' )
  !
  ! ... defaults
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM(outdir) == ' ' ) outdir = './'
  lmax = 12
  band_index = 0
  magmom(:) = 0.0_DP
  moment_tol = 1.0d-3
  report_file = 'harmwave_report.txt'
  !
  ios = 0
  IF ( ionode ) THEN
     CALL input_from_file( )
     READ( 5, inputpp, iostat = ios )
     tmp_dir = trimcheck( outdir )
  END IF
  CALL mp_bcast( ios, ionode_id, intra_image_comm )
  IF ( ios /= 0 ) CALL errore( 'harmwave', 'reading inputpp namelist', ABS(ios) )
  !
  CALL mp_bcast( tmp_dir,     ionode_id, intra_image_comm )
  CALL mp_bcast( prefix,      ionode_id, intra_image_comm )
  CALL mp_bcast( lmax,        ionode_id, intra_image_comm )
  CALL mp_bcast( band_index,  ionode_id, intra_image_comm )
  CALL mp_bcast( magmom,      ionode_id, intra_image_comm )
  CALL mp_bcast( moment_tol,  ionode_id, intra_image_comm )
  CALL mp_bcast( report_file, ionode_id, intra_image_comm )
  !
  ! ... load the finished 'bands' run exactly as QE itself left it:
  ! ... et, xk, nkstot, nspin, tau, ityp, at, bg, alat all populated.
  !
  CALL read_file( )
  !
  IF ( nspin /= 2 .OR. noncolin ) CALL errore( 'harmwave', &
       'v1 handles collinear LSDA only (nspin=2, noncolin=.false.)', 1 )
  IF ( MOD(nkstot,2) /= 0 ) CALL errore( 'harmwave', &
       'expected an even nkstot (LSDA doubles the k-list); got an odd one', nkstot )
  IF ( band_index < 1 .OR. band_index > nbnd ) CALL errore( 'harmwave', &
       'band_index out of range 1..nbnd', band_index )
  IF ( nat > 1000 ) CALL errore( 'harmwave', 'raise the hardcoded atom limit', nat )
  !
  n = nkstot / 2       ! k-points per spin channel; ik=1..n up, ik=n+1..2n down
  !
  IF ( n < 3 ) CALL errore( 'harmwave', 'need at least 3 k-points to define a plane', n )
  !
  !=====================================================================
  ! Part 0: the circle's plane, taken FROM THE K-POINTS rather than assumed.
  !
  ! v1 was hardcoded to the k_x k_y plane, which is where this project's own
  ! circles happen to lie and is not where another material's will. The plane
  ! is over-determined by the k-list, so it is read off it: the normal from
  ! the first three non-collinear points, then every point checked against
  ! it. This has to run BEFORE Part 1, because the symmetry filter below
  ! tests operations against this normal.
  !
  ! The in-plane frame (e1,e2) sets where theta = 0 and which way it turns.
  ! Amplitudes do not care, individual a_l and b_l do, so the frame is
  ! written into the report rather than left implicit. When the normal is
  ! +z the frame is (x,y) exactly, so a k_x k_y circle reproduces v1's
  ! numbers bit for bit.
  !=====================================================================
  !
  v1(:) = xk(:,2) - xk(:,1)
  nplane(:) = 0.0_DP
  DO ik3 = 3, n
     v2(:) = xk(:,ik3) - xk(:,1)
     vn(1) = v1(2)*v2(3) - v1(3)*v2(2)
     vn(2) = v1(3)*v2(1) - v1(1)*v2(3)
     vn(3) = v1(1)*v2(2) - v1(2)*v2(1)
     nrm = SQRT( SUM( vn(:)**2 ) )
     IF ( nrm > 1.0d-8 ) THEN
        nplane(:) = vn(:) / nrm
        EXIT
     END IF
  END DO
  IF ( SUM( nplane(:)**2 ) < 0.5_DP ) CALL errore( 'harmwave', &
       'the k-points are collinear: no plane, so no angular harmonics', 1 )
  !
  ! Orient the normal so its largest component is positive: the sign of the
  ! cross product depends on the k-list's order, and theta's sense should not.
  IF ( nplane(MAXLOC(ABS(nplane),1)) < 0.0_DP ) nplane(:) = -nplane(:)
  !
  c_off = SUM( nplane(:) * xk(:,1) )
  DO ik = 2, n
     IF ( ABS( SUM( nplane(:) * xk(:,ik) ) - c_off ) > 1.0d-6 ) &
        CALL errore( 'harmwave', 'k-point off the plane the first three define', ik )
  END DO
  !
  ! The in-plane frame. For a +z normal this is exactly (x,y).
  IF ( ABS(nplane(3) - 1.0_DP) < 1.0d-8 ) THEN
     e1 = (/ 1.0_DP, 0.0_DP, 0.0_DP /)
     e2 = (/ 0.0_DP, 1.0_DP, 0.0_DP /)
  ELSE
     ! e1 from the least-aligned cartesian axis, projected into the plane
     j = MINLOC( ABS(nplane), 1 )
     e1(:) = 0.0_DP ; e1(j) = 1.0_DP
     e1(:) = e1(:) - SUM( e1(:)*nplane(:) ) * nplane(:)
     e1(:) = e1(:) / SQRT( SUM( e1(:)**2 ) )
     e2(1) = nplane(2)*e1(3) - nplane(3)*e1(2)
     e2(2) = nplane(3)*e1(1) - nplane(1)*e1(3)
     e2(3) = nplane(1)*e1(2) - nplane(2)*e1(1)
  END IF
  !
  !=====================================================================
  ! Part 1: the non-magnetic parent's symmetry, from QE's own find_sym(),
  ! with the two (or more) magnetic sublattices collapsed to their TRUE
  ! chemical species -- read from the pseudopotential, not from ityp.
  !=====================================================================
  !
  n_chem = 0
  DO it = 1, nsp
     is_dup = .FALSE.
     DO j = 1, n_chem
        IF ( chem_label(j) == upf(it)%psd ) THEN
           is_dup = .TRUE.
           EXIT
        END IF
     END DO
     IF ( .NOT. is_dup ) THEN
        n_chem = n_chem + 1
        IF ( n_chem > max_ntyp ) CALL errore( 'harmwave', 'raise max_ntyp', n_chem )
        chem_label(n_chem) = upf(it)%psd
     END IF
  END DO
  DO ia = 1, nat
     DO j = 1, n_chem
        IF ( chem_label(j) == upf(ityp(ia))%psd ) chem_ityp(ia) = j
     END DO
  END DO
  !
  dummy_mloc(:,1:nat) = 0.0_DP
  CALL find_sym( nat, tau, chem_ityp(1:nat), .FALSE., dummy_mloc(:,1:nat), .FALSE. )
  ! find_sym() has just populated (module symm_base): nsym, s(:,:,1:nsym),
  ! sr(:,:,1:nsym) (Cartesian, via its own internal s_axis_to_cart() call),
  ! and irt(1:nsym,1:nat) -- the atom each operation maps atom ia onto.
  !
  ! First pass: classify EVERY crystallographic operation's character,
  ! independent of whether it happens to preserve this particular circle's
  ! plane. Getting this order right matters: the operation that actually
  ! defines the altermagnet -- e.g. a horizontal mirror that exchanges the
  ! two sublattices -- typically flips k_z and so fails the plane filter on
  ! its own, but combined with time reversal (k -> -k) it does not, and that
  ! combination must not be missed. Building the full, time-reversal-closed
  ! set FIRST and filtering for the plane SECOND (mirroring
  ! altharm.symmetry's own order of operations) is what a first version of
  ! this routine got backwards: it time-reversed only the operations that
  ! had already survived the plane filter, silently dropping every
  ! kz-flipping operation's time-reversed partner and so over-predicting
  ! the allowed harmonics (found {0,3,6,9,12} on Cr2Se2As instead of the
  ! correct {3,9,15,...}, checked against altharm.symmetry.allowed_harmonics
  ! fed this routine's own operator dump).
  !
  ALLOCATE( op_cart(3,3,max_ops), op_chi(max_ops) )
  nkept = 0
  !
  DO isym = 1, nsym
     ok_plus  = .TRUE.
     ok_minus = .TRUE.
     DO ia = 1, nat
        IF ( ABS(magmom(ityp(ia))) < moment_tol ) CYCLE
        ja = irt(isym,ia)
        si = SIGN(1.0_DP, magmom(ityp(ia)))
        sj = SIGN(1.0_DP, magmom(ityp(ja)))
        IF ( sj /=  si ) ok_plus  = .FALSE.
        IF ( sj /= -si ) ok_minus = .FALSE.
     END DO
     IF ( ok_plus .EQV. ok_minus ) CYCLE        ! genuinely broken by the
                                                 ! magnetic order either way
     !
     ! chi for the bare operation, and for it combined with time reversal
     ! (k -> -k; character unchanged, since time reversal alone is chi=+1
     ! for a collinear, spin-orbit-free calculation -- each spin channel is
     ! separately time-reversal symmetric). Try BOTH against the plane
     ! filter: whichever of R or -R leaves k_z pointwise fixed is the one
     ! that constrains this circle.
     !
     CALL keep_if_plane_preserving( sr(:,:,isym),  MERGE(1,-1,ok_plus) )
     CALL keep_if_plane_preserving( -sr(:,:,isym), MERGE(1,-1,ok_plus) )
  END DO
  !
  !=====================================================================
  ! Part 2: the raw harmonic fit, on QE's own in-memory eigenvalues.
  !=====================================================================
  !
  ! the plane was established in Part 0; theta is measured in its own frame
  ALLOCATE( theta(n), dE(n) )
  DO ik = 1, n
     theta(ik) = ATAN2( SUM( e2(:)*xk(:,ik) ), SUM( e1(:)*xk(:,ik) ) )
     dE(ik) = ( et(band_index,ik) - et(band_index,ik+n) ) * RYTOEV
  END DO
  !
  ncoef = 2*lmax + 1
  ALLOCATE( design(n,ncoef), coef(MAX(n,ncoef)) )
  design(:,1) = 1.0_DP
  DO l = 1, lmax
     design(:,2*l)   = COS( l * theta(:) )
     design(:,2*l+1) = SIN( l * theta(:) )
  END DO
  coef(1:n) = dE(1:n)
  !
  lwork = MAX( 1, 3*MIN(n,ncoef) + MAX(2*MIN(n,ncoef), MAX(n,ncoef)) ) + 64*MAX(n,ncoef)
  ALLOCATE( work(lwork) )
  CALL DGELSS( n, ncoef, 1, design, n, coef, MAX(n,ncoef), s_sv, -1.0_DP, &
               j, work, lwork, info )
  IF ( info /= 0 ) CALL errore( 'harmwave', 'DGELSS failed', info )
  !
  !=====================================================================
  ! Report -- everything an independent Python check needs, nothing it
  ! has to trust blindly: raw (theta, dE) samples, the kept (R_cart, chi)
  ! list, and the raw fitted coefficients.
  !=====================================================================
  !
  IF ( ionode ) THEN
     iunit = 77
     OPEN( UNIT=iunit, FILE=TRIM(report_file), STATUS='unknown' )
     WRITE( iunit, '(A)' ) '# harmwave.x v1 -- direct-from-QE harmonic fit'
     WRITE( iunit, '(A,I0)' ) '# nat = ', nat
     WRITE( iunit, '(A,I0)' ) '# nsp (species as given to pw.x) = ', nsp
     WRITE( iunit, '(A,I0)' ) '# n_chem (true chemical species, collapsed) = ', n_chem
     WRITE( iunit, '(A,F18.12)' ) '# alat (bohr) = ', alat
     WRITE( iunit, '(A,F18.12)' ) '# RYTOEV = ', RYTOEV
     WRITE( iunit, '(A,I0)' ) '# n (k-points per spin channel) = ', n
     WRITE( iunit, '(A,I0)' ) '# band_index = ', band_index
     WRITE( iunit, '(A,I0)' ) '# lmax = ', lmax
     WRITE( iunit, '(A,3F16.9)' ) '# plane normal (cartesian, 2pi/alat) = ', nplane(:)
     WRITE( iunit, '(A,F16.9)' )  '# plane offset  n.k                  = ', c_off
     WRITE( iunit, '(A,3F16.9)' ) '# theta frame e1                     = ', e1(:)
     WRITE( iunit, '(A,3F16.9)' ) '# theta frame e2                     = ', e2(:)
     WRITE( iunit, '(A,I0)' ) '# nsym (crystallographic, before plane filter) = ', nsym
     WRITE( iunit, '(A,I0)' ) '# nkept (plane-preserving, chi-classified ops) = ', nkept
     WRITE( iunit, '(A)' ) '#'
     WRITE( iunit, '(A)' ) '# OPERATIONS: chi  R11 R12 R13 R21 R22 R23 R31 R32 R33'
     DO isym = 1, nkept
        WRITE( iunit, '(A,I3,9F14.9)' ) 'OP  ', op_chi(isym), &
             op_cart(1,1,isym), op_cart(1,2,isym), op_cart(1,3,isym), &
             op_cart(2,1,isym), op_cart(2,2,isym), op_cart(2,3,isym), &
             op_cart(3,1,isym), op_cart(3,2,isym), op_cart(3,3,isym)
     END DO
     WRITE( iunit, '(A)' ) '#'
     WRITE( iunit, '(A)' ) '# SAMPLES: theta(rad)  dE(eV) [E_up - E_down, band_index]'
     DO ik = 1, n
        WRITE( iunit, '(A,2F20.14)' ) 'PT  ', theta(ik), dE(ik)
     END DO
     WRITE( iunit, '(A)' ) '#'
     WRITE( iunit, '(A)' ) '# COEFFICIENTS: l  kind(c|s)  value(eV)'
     WRITE( iunit, '(A,I3,A,F20.14)' ) 'CF  ', 0, '  c  ', coef(1)
     DO l = 1, lmax
        WRITE( iunit, '(A,I3,A,F20.14)' ) 'CF  ', l, '  c  ', coef(2*l)
        WRITE( iunit, '(A,I3,A,F20.14)' ) 'CF  ', l, '  s  ', coef(2*l+1)
     END DO
     CLOSE( iunit )
     WRITE( stdout, '(5X,A,A)' ) 'harmwave: report written to ', TRIM(report_file)
  END IF
  !
  DEALLOCATE( op_cart, op_chi, theta, dE, design, coef, work )
  !
  CALL environment_end( 'HARMWAVE' )
  CALL stop_pp
  STOP
  !
CONTAINS
  !
  SUBROUTINE keep_if_plane_preserving( rcand, chi_cand )
    !! Append (rcand, chi_cand) to the kept list iff rcand leaves the
    !! component along the circle's normal fixed pointwise, and it is not
    !! already present with the same character.
    !!
    !! k' = R k, so n.k is preserved for every k iff R^T n = n. For n = +z
    !! that reduces to row 3 of R being (0,0,1), which is what v1 tested
    !! directly -- the general form is the same condition without the
    !! assumption. The plane's OFFSET needs no separate check: if R^T n = n
    !! then n.(Rk) = (R^T n).k = n.k = c_off automatically.
    REAL(DP), INTENT(IN) :: rcand(3,3)
    INTEGER,  INTENT(IN) :: chi_cand
    INTEGER :: jj
    LOGICAL :: dup
    IF ( ANY( ABS( MATMUL( TRANSPOSE(rcand), nplane ) - nplane ) > 1.0d-6 ) ) RETURN
    dup = .FALSE.
    DO jj = 1, nkept
       IF ( op_chi(jj) == chi_cand .AND. &
            ALL( ABS(op_cart(:,:,jj) - rcand) < 1.0d-8 ) ) dup = .TRUE.
    END DO
    IF ( dup ) RETURN
    nkept = nkept + 1
    IF ( nkept > max_ops ) CALL errore( 'harmwave', 'raise max_ops', nkept )
    op_cart(:,:,nkept) = rcand
    op_chi(nkept) = chi_cand
  END SUBROUTINE keep_if_plane_preserving
  !
END PROGRAM harmwave
