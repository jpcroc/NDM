module spline_mod
        !implicit none
        contains 

subroutine cspline(n, x, y, b, c, d)
  !-----------------------------------------------
  !   M o d u l e s
  !-----------------------------------------------
  USE T_kind_param_m
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer , intent(in) :: n
  real(double) , intent(in) :: x(n)
  real(double) , intent(in) :: y(n)
  real(double) , intent(inout) :: b(n)
  real(double) , intent(inout) :: c(n)
  real(double) , intent(inout) :: d(n)
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: nm1, ib, i
  real(double) :: t
  !-----------------------------------------------
  !
  !  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
  !  for a cubic interpolating spline
  !
  !    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
  !
  !    for  x(i) .le. x .le. x(i+1)
  !
  !  input..
  !
  !    n = the number of data points or knots (n.ge.2)
  !    x = the abscissas of the knots in strictly increasing order
  !    y = the ordinates of the knots
  !
  !  output..
  !
  !    b, c, d  = arrays of spline coefficients as defined above.
  !
  !  using  p  to denote differentiation,
  !
  !    y(i) = s(x(i))
  !    b(i) = sp(x(i))
  !    c(i) = spp(x(i))/2
  !    d(i) = sppp(x(i))/6  (derivative from the right)
  !
  !  the accompanying function subprogram  seval  can be USEd
  !  to evaluate the spline.
  !
  !
  !
! write(6,*)'spline'
  nm1 = n-1
  if (n<2) return
  if (n>=3) then
     !
     !  set up tridiagonal system
     !
     !  b = diagonal, d = offdiagonal, c = right hand side.
     !
     d(1) = x(2)-x(1)
     c(2) = (y(2)-y(1))/d(1)
     d(2:nm1) = x(3:nm1+1)-x(2:nm1)
     b(2:nm1) = 2.*(d(:nm1-1)+d(2:nm1))
     c(3:nm1+1) = (y(3:nm1+1)-y(2:nm1))/d(2:nm1)
     c(2:nm1) = c(3:nm1+1)-c(2:nm1)
     !
     !  end conditions.  third derivatives at  x(1)  and  x(n)
     !  obtained from divided differences
     !
     b(1) = -d(1)
     b(n) = -d(n-1)
     c(1) = 0.
     c(n) = 0.
     if (n/=3) then
        c(1) = c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
        c(n) = c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
        c(1) = c(1)*d(1)**2/(x(4)-x(1))
        c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
     endif
     !
     !  forward elimination
     !
     do i = 2, n
        t = d(i-1)/b(i-1)
        b(i) = b(i)-t*d(i-1)
        c(i) = c(i)-t*c(i-1)
     end do
     !
     !  back substitution
     !
     c(n) = c(n)/b(n)
     do ib = 1, nm1
        i = n-ib
        c(i) = (c(i)-d(i)*c(i+1))/b(i)
     end do
     !
     !  c(i) is now the sigma(i) of the text
     !
     !  compute polynomial coefficients
     !
     b(n) = (y(n)-y(nm1))/d(nm1)+d(nm1)*(c(nm1)+2.*c(n))
     b(:nm1) = (y(2:nm1+1)-y(:nm1))/d(:nm1)-d(:nm1)*(c(2:nm1+1)+2.*c(:nm1))
     d(:nm1) = (c(2:nm1+1)-c(:nm1))/d(:nm1)
     c(:nm1) = 3.*c(:nm1)
     c(n) = 3.*c(n)
     d(n) = d(n-1)
     return
  endif
  !
  b(1) = (y(2)-y(1))/(x(2)-x(1))
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end subroutine cspline


      SUBROUTINE DCSSMO(H, N, TNODE, G, WGS, RHO, GSMO, B, C, D)       !DCS   10
!
!  THIS SUBROUTINE COMPUTES THE DISCRETE NATURAL CUBIC
!  SPLINE DEFINED ON THE INTERVAL (TNODE(1),TNODE(N)) WHICH
!  SMOOTHS THROUGH THE DATA (TNODE(I),G(I)),I=1,2,...,N.
!  N MUST BE 2 OR GREATER. THE NODES MUST SATISFY TNODE(I)
!  .LT.TNODE(I+1). THE SOLUTION S(T) FOR T IN THE INTERVAL
!  (TNODE(I),TNODE(I+1)) IS GIVEN BY
!
!     S(T)=GSMO(I)+B(I)*(T-TNODE(I))+
!             C(I)*(T-TNODE(I))**2+D(I)*(T-TNODE(I))**3
!
      real(kind(1.d0))  :: TNODE(N), G(N), WGS(N), GSMO(N),& 
                           B(N), C(N), D(N)
      real(kind(1.0)) :: RHO
!
!  INPUT  PARAMETERS(NONE OF THE INPUT PARAMETERS ARE CHANGED
!         BY THIS SUBROUTINE)
!
!  H     - THE STEP SIZE USED FOR THE DISCRETE CUBIC SPLINE
!  N     - NUMBER OF NODES (TNODE) AND DATA VALUES(G)
!  TNODE - REAL ARRAY CONTAINING THE NODES (TNODE(I).LT.
!          TNODE(I+1)).
!  G     - REAL ARRAY CONTAINING THE DATA VALUES.
!  WGS   - REAL ARRAY CONTAINING THE WEIGHTS WGS(I)
!          CORRESPONDING TO THE DATA (TNODE(I),G(I)).
!  RHO   - SIMPLE REAL VARIABLE CONTAINING THE POSITIVE
!          PARAMETER FOR VARYING THE SMOOTHNESS OF THE FIT.
!          IF RHO IS SMALL SMOOTHNESS IS EMPHASIZED.
!          IF RHO IS LARGE DATA FITTING IS EMPHASIZED.
!
!  OUTPUT PARAMETERS
!
!  GSMO  - REAL ARRAY CONTAINING THE SMOOTHED VALUES OF
!          THE DATA G(I),I=1,2,....,N.
!  B     - REAL ARRAY CONTAINING THE COEFFICIENTS B(I) FOR
!          THE TERMS (T-TNODE(I)).
!  C     - REAL ARRAY CONTAINING THE COEFFICIENTS C(I) FOR
!          THE TERMS (T-TNODE(I))**2.
!  D     - REAL ARRAY CONTAINING THE COEFFICIENTS D(I) FOR
!          THE TERMS (T-TNODE(I))**3.
!
      IF (N.EQ.2) GO TO 180
      N1 = N - 1
      N2 = N1 - 1
      N3 = N2 - 1
!  THE RIGHT HAND SIDE OF THE LINEAR SYSTEM FOR THE
!  C(I)'S WILL NOW BE CONSTRUCTED.
      DO 10 I=1,N
        C(I) = G(I)
   10 CONTINUE
      DO 20 I=1,N1
        C(I) = (C(I+1)-C(I))/(TNODE(I+1)-TNODE(I))
   20 CONTINUE
      DO 30 I=1,N2
        C(I) = 3.0*(C(I+1)-C(I))
   30 CONTINUE
!  THE RIGHT HAND SIDE IS NOW IN ARRAY C.
!
!  THE P.D. 5 BANDED SYMMETRIC MATRIX WILL NOW BE CONSTRUCTED.
!  THE THREE NEEDED DIAGONALS WILL BE STORED IN ARRAYS
!  GSMO,B,D.
      H2 = H*H
      H3 = H2*H
      R6 = 6.0*H3/RHO
      HI3 = TNODE(2) - TNODE(1)
      HI4 = TNODE(3) - TNODE(2)
      ETA3 = HI3 + HI3 + H2/HI3
      BETA3 = R6/(WGS(1)*HI3)
      BETA4 = R6/(WGS(2)*HI3*HI4)
      EPS3 = (BETA3+BETA4*HI4)/HI3
      H2DHI = H2/HI4
      ETA4 = HI4 + HI4 + H2DHI
      IF (N.EQ.3) GO TO 60
      HI5 = TNODE(4) - TNODE(3)
      BETA5 = R6/(WGS(3)*HI4*HI5)
      EPS4 = (BETA4*HI3+BETA5*HI5)/HI4
      GSMO(1) = ETA3 + ETA4 + BETA4 + BETA4 + EPS3 + EPS4
      P = H2DHI + BETA4 + BETA5 + EPS4
      B(1) = HI4 - P
      IF (N.EQ.4) GO TO 50
      DO 40 I=2,N3
        HI3 = HI4
        HI4 = HI5
        HI5 = TNODE(I+3) - TNODE(I+2)
        ETA3 = ETA4
        H2DHI = H2/HI4
        ETA4 = HI4 + HI4 + H2DHI
        BETA3 = BETA4
        BETA4 = BETA5
        BETA5 = R6/(WGS(I+2)*HI4*HI5)
        EPS3 = EPS4
        EPS4 = (BETA4*HI3+BETA5*HI5)/HI4
        D(I-1) = BETA4
        P = H2DHI + BETA4 + BETA5 + EPS4
        B(I) = HI4 - P
        GSMO(I) = ETA3 + ETA4 + BETA4 + BETA4 + EPS3 + EPS4
   40 CONTINUE
   50 HI3 = HI4
      HI4 = HI5
      ETA3 = ETA4
      ETA4 = HI4 + HI4 + H2/HI4
      BETA4 = BETA5
      EPS3 = EPS4
   60 BETA5 = R6/(WGS(N)*HI4)
      EPS4 = (BETA4*HI3+BETA5)/HI4
      GSMO(N2) = ETA3 + ETA4 + BETA4 + BETA4 + EPS3 + EPS4
!  THE P.D. 5 BANDED SYMMETRIC MATRIX IS COMPLETE.
!  THE SYSTEM OF LINEAR EQUATION WILL NOW BE SOLVED FOR THE
!  C(I)'S.
      IF (N.GT.3) GO TO 70
      C(1) = C(1)/GSMO(1)
      GO TO 150
   70 IF (N.GT.4) GO TO 80
      C(1) = (C(1)*GSMO(2)-C(2)*B(1))/(GSMO(1)*GSMO(2)-B(1)**2)
      C(2) = (C(2)-C(1)*B(1))/GSMO(2)
      GO TO 150
!  THIS SOLVE THE 5 BANDED SYSTEM WHEN K=N-2.GT.3.
   80 K = N2
      K1 = K - 1
      K2 = K1 - 1
      K3 = K2 - 1
!  THE 5 BANDED MATRIX WILL NOW BE FACTORED.
      B(1) = B(1)/GSMO(1)
      D(1) = D(1)/GSMO(1)
      P = GSMO(1)*B(1)
      GSMO(2) = GSMO(2) - P*B(1)
      B(2) = (B(2)-P*D(1))/GSMO(2)
      IF (K.EQ.3) GO TO 110
      D(2) = D(2)/GSMO(2)
      IF (K.EQ.4) GO TO 100
      DO 90 I=3,K2
        I1 = I - 1
        I2 = I1 - 1
        P = GSMO(I1)*B(I1)
        GSMO(I) = GSMO(I) - GSMO(I2)*(D(I2)**2) - P*B(I1)
        B(I) = (B(I)-P*D(I1))/GSMO(I)
        D(I) = D(I)/GSMO(I)
   90 CONTINUE
  100 P = GSMO(K2)*B(K2)
      GSMO(K1) = GSMO(K1) - GSMO(K3)*(D(K3)**2) - P*B(K2)
      B(K1) = (B(K1)-P*D(K2))/GSMO(K1)
  110 GSMO(K) = GSMO(K) - GSMO(K2)*(D(K2)**2) - GSMO(K1)*(B(K1)**2)
!  FACTORIZATION COMPLETE.
!  CARRY OUT FORWARD  AND BACKWARD SUBSTITUTION.
      C(2) = C(2) - B(1)*C(1)
      DO 120 I=3,K
        I1 = I - 1
        I2 = I - 2
        C(I) = C(I) - B(I1)*C(I1) - D(I2)*C(I2)
  120 CONTINUE
      DO 130 I=1,K
        C(I) = C(I)/GSMO(I)
  130 CONTINUE
      C(K1) = C(K1) - B(K1)*C(K)
      DO 140 I=2,K1
        J = K - I
        C(J) = C(J) - B(J)*C(J+1) - D(J)*C(J+2)
  140 CONTINUE
!  THE 5 BANDED SYSTEM HAS BEEN SOLVED.THE SOLUTION IS IN
!  ARRAY C. THE COEFFICIENTS GSMO, B, C, AND D WILL NOW BE
!  SET UP.
  150 C(N) = 0.0
      D(N) = 0.0
      C(N1) = C(N2)
      HK1 = TNODE(N) - TNODE(N1)
      D(N1) = -C(N1)/(3.0*HK1)
      GSMO(N) = G(N) + R6*D(N1)/WGS(N)
      IF (N.EQ.3) GO TO 170
      DO 160 I=2,N2
        K = N - I
        K1 = K + 1
        HK2 = HK1
        HK1 = TNODE(K1) - TNODE(K)
        C(K) = C(K-1)
        D(K) = (C(K1)-C(K))/(3.0*HK1)
        GSMO(K1) = G(K1) - R6*(D(K1)-D(K))/WGS(K1)
        B(K1) = (GSMO(K1+1)-GSMO(K1))/HK2 - HK2*(C(K1)+C(K1)+C(K1+1))/ &
         3.0
  160 CONTINUE
  170 C(1) = 0.0
      HK2 = HK1
      HK1 = TNODE(2) - TNODE(1)
      D(1) = (C(2)-C(1))/(3.0*HK1)
      GSMO(2) = G(2) - R6*(D(2)-D(1))/WGS(2)
      GSMO(1) = G(1) - R6*D(1)/WGS(1)
      B(2) = (GSMO(3)-GSMO(2))/HK2 - HK2*(C(2)+C(2)+C(3))/3.0
      B(1) = (GSMO(2)-GSMO(1))/HK1 - HK1*(C(1)+C(1)+C(2))/3.0
!  THE DISCRETE CUBIC SMOOTHING SPLINE IS NOW COMPLETE.
      RETURN
!  THE TRIVIAL CASE WHEN N=2 IS HANDLED HERE.
  180 GSMO(1) = G(1)
      GSMO(2) = G(2)
      B(1) = (G(2)-G(1))/(TNODE(2)-TNODE(1))
      C(1) = 0.0
      D(1) = 0.0
      RETURN
      END
      SUBROUTINE DCSINT(IENT, H, N, TNODE, G, END1, ENDN, B, C, D)     
!
!  THIS SUBROUTINE COMPUTES THE DISCRETE CUBIC SPLINE
!  DEFINED ON THE INTERVAL (TNODE(1),TNODE(N)),WHICH INTER-
!  POLATES THE DATA (TNODE(I),G(I)),I=1,2,...,N. WE REQUIRE
!  THAT TNODE(I).LT.TNODE(I+1). END1 AND ENDN CONTAIN THE
!  VALUES OF THE END CONDITIONS BEING USED.
!
!  IF IENT=1,THE FIRST CENTRAL DIVIDED DIFFERENCE END
!  CONDITIONS ARE BEING USED.
!
!  IF IENT=2,THE SECOND CENTRAL DIVIDED DIFFERENCE END
!  CONDITIONS ARE BEING USED.
!
!  IF IENT=3,THE PERODIC END CONDITIONS ARE BEING USED.
!  FOR THIS CASE THE CONTENTS OF G(N), END1,AND ENDN ARE
!  IGNORED.
!
!  FOR ALL THREE END CONDITIONS N MUST BE GREATER THAN OR
!  EQUAL TO 2.
!
!  THE DISCRETE CUBIC SPLINE IS REPRESENTED BY PIECEWISE
!  CUBIC POLYNOMIALS. FOR T IN THE INTERNAL (TNODE(I),TNODE
!  (I+1)) THE CUBIC SPLINE IS
!
!     S(T)=G(I)+B(I)*(T-TNODE(I))
!                          +C(I)*(T-TNODE(I))**2
!                                   +D(I)*(T-TNODE(I))**3
!
      USE  T_kind_param_m, ONLY : double
      real(double) :: TNODE(N), G(N), B(N), C(N), D(N)
      real(double) :: H,END1,ENDN
!
!  INPUT PARAMETERS (NONE OF THESE PARAMETERS
!                          ARE CHANGED BY THIS SUBROUTINE.)
!
!  IENT  - SPECIFIES END CONDITIONS WHICH ARE IN EFFECT.
!  H     - THE STEP SIZE USED FOR THE DISCRETE CUBIC SPLINE.
!  N     - NUMBER OF NODES (TNODE) AND DATA VALUES (G).
!          (N.GE.2)
!  TNODE - REAL ARRAY CONTAINING THE NODES (TNODE(I).LT.
!          TNODE(I+1)).
!  G     - REAL ARRAY CONTAINING THE INTERPOLATING DATA.
!  END1  - END CONDITION VALUE AT TNODE(1).
!  ENDN  - END CONDITION VALUE AT TNODE(N).
!
!  OUTPUT PARAMETERS
!
!  B     - REAL ARRAY CONTAINING COEFFICIENTS OF
!          (T-TNODE(I)),I=1,2,....,N-1.
!  C     - REAL ARRAY CONTAINING COEFFICIENTS OF
!          (T-TNODE(I))**2,I=1,2,...,N-1.
!  D     - REAL ARRAY CONTAINING COEFFICIENTS OF
!          (T-TNODE(I))**3,I=1,2,...,N-1.
!
!  SPECIAL CASES ARE ACCOUNTED FOR HERE.
      IF (N.EQ.2 .AND. IENT.EQ.3) GO TO 220
      H2 = H*H
      N1 = N - 1
      IF (N.EQ.2 .AND. IENT.EQ.2) GO TO 180
      N2 = N1 - 1
!  THE SYMMETRIC TRIDIAGONAL(OR NEAR TRIDIAGONAL) LINEAR
!  SYSTEM WILL NOW BE SET UP FOR THE APPROPRIATE END
!  CONDITIONS
      HI = TNODE(2) - TNODE(1)
      H2DHI = H2/HI
      ETA2 = HI + HI + H2DHI
      GO TO (10, 40, 60), IENT
!  IENT=1 - FIRST CENTRAL DIVIDED DIFFERENCE END CONDITONS
!  SET UP.
   10 B(1) = ETA2
      D(1) = HI - H2DHI
      G2 = (G(2)-G(1))/HI
      C(1) = 3.0*(G2-END1)
      IF (N.EQ.2) GO TO 30
      DO 20 I=2,N1
        ETA1 = ETA2
        G1 = G2
        HI = TNODE(I+1) - TNODE(I)
        H2DHI = H2/HI
        ETA2 = HI + HI + H2DHI
        B(I) = ETA1 + ETA2
        D(I) = HI - H2DHI
        G2 = (G(I+1)-G(I))/HI
        C(I) = 3.0*(G2-G1)
   20 CONTINUE
   30 B(N) = ETA2
      C(N) = 3.0*(ENDN-G2)
      L = N
!  SET UP FOR (1) FIRST CENTRAL DIVIDED DIFFERENCE END
!  CONDITING COMPLETE.THE LINEAR EQUATIONS WILL BE NOW
!  SOLVED.
      GO TO 110
!  IENT=2 - SECOND CENTRAL DIVIDED DIFFERENCE END CONDITIONS
!  SET UP.
   40 GAMMA = HI - H2DHI
      G2 = (G(2)-G(1))/HI
      DO 50 I=1,N2
        ETA1 = ETA2
        G1 = G2
        HI = TNODE(I+2) - TNODE(I+1)
        H2DHI = H2/HI
        ETA2 = HI + HI + H2DHI
        B(I) = ETA1 + ETA2
        D(I) = HI - H2DHI
        G2 = (G(I+2)-G(I+1))/HI
        C(I) = 3.0*(G2-G1)
   50 CONTINUE
      C(1) = C(1) - GAMMA*END1/2.0
      HI = TNODE(N) - TNODE(N1)
      GAMMA = HI - H2/HI
      C(N2) = C(N2) - GAMMA*ENDN/2.0
!  STEP UP FOR (2) SECOND CENTRAL DIVIDED DIFFERENCE
!  END CONDITIONS COMPLETE. THE LINEAR EQUATIONS WILL NOW
!  BE SOLVED.
      IF (N2.EQ.1) GO TO 100
      L = N2
      GO TO 110
!  IENT=3 - PERIODIC END CONDITIONS SET UP.
   60 B(1) = ETA2
      D(1) = HI - H2DHI
      DO 70 I=2,N1
        ETA1 = ETA2
        HI = TNODE(I+1) - TNODE(I)
        H2DHI = H2/HI
        ETA2 = HI + HI + H2DHI
        B(I) = ETA1 + ETA2
        D(I) = HI - H2DHI
        C(I) = 0.0
   70 CONTINUE
      CT = D(N1)
      C(1) = CT
      C(N1) = CT
      B(1) = B(1) + ETA2 - CT
      B(N1) = B(N1) - CT
      L = N1
      ITRANS = 1
      GO TO 120
   80 G1 = (G(N1)-G(N2))/(TNODE(N1)-TNODE(N2))
      G2 = (G(1)-G(N1))/(TNODE(N)-TNODE(N1))
      DEN = (1.0+C(1)+C(N1))
      CT = 3.0*(G2-G1)
      BS1 = CT*C(N1)
      C(N1) = CT
      DO 90 I=1,N2
        HI = TNODE(I+1) - TNODE(I)
        G1 = G2
        G2 = (G(I+1)-G(I))/HI
        CT = 3.0*(G2-G1)
        BS1 = BS1 + CT*C(I)
        C(I) = CT
   90 CONTINUE
      BS1 = BS1/DEN
      C(1) = C(1) - BS1
      C(N1) = C(N1) - BS1
      ITRANS = 0
      GO TO 140
!  THE SET UP AND MOST OF THE DETAILS FOR SOLVING THE
!  LINEAR EQUQTION FOR (3)  THE PERIODIC END CONDITIONS
!  ARE COMPLETED.
!
!  THE LINEAR EQUATION ARE SOLVED HERE.
  100 C(2) = C(1)/B(1)
      GO TO 180
  110 ITRANS = 0
  120 L1 = L - 1
      DO 130 I=1,L1
        T = D(I)/B(I)
        B(I+1) = B(I+1) - T*D(I)
        D(I) = T
  130 CONTINUE
!  THE LINEAR EQUATION SOLVER IS ENTERED AT THIS POINT
!  IF THE LU FACTORIZATION HAS ALREADY BEEN DONE.
  140 DO 150 I=1,L1
        C(I+1) = C(I+1) - D(I)*C(I)
  150 CONTINUE
      C(L) = C(L)/B(L)
      DO 160 I=1,L1
        LI = L - I
        C(LI) = C(LI)/B(LI) - D(LI)*C(LI+1)
  160 CONTINUE
      IF (ITRANS.GE.1) GO TO 80
!  THE LINEAR SYSTEM HAS BEEN SOLVE FOR THE C-VECTOR
      IF (IENT.EQ.3) C(N) = C(1)
      IF (IENT.NE.2) GO TO 190
      DO 170 I=1,N2
        LI = N - I
        C(LI) = C(LI-1)
  170 CONTINUE
  180 C(1) = END1/2.0
      C(N) = ENDN/2.0
  190 C1 = C(1)
      C2 = C(2)
      HI = TNODE(2) - TNODE(1)
      IF (N.EQ.2) GO TO 210
      DO 200 I=1,N2
        B(I) = (G(I+1)-G(I))/HI - HI*(C1+C1+C2)/3.0
        D(I) = (C2-C1)/(3.0*HI)
        HI = TNODE(I+2) - TNODE(I+1)
        C1 = C2
        C2 = C(I+2)
  200 CONTINUE
  210 GN = G(1)
      IF (IENT.NE.3) GN = G(N)
      B(N1) = (GN-G(N1))/HI - HI*(C1+C1+C2)/3.0
      D(N1) = (C2-C1)/(3.0*HI)
!  THE INTERPOLATING DISCRETE CUBIC SPLINE HAS BEEN
!  CONSTRUCTED.
      RETURN
!  THE FOLLOWING HANDLES THE TRIVIAL PERIODIC CASE
!  (IENT.EQ.3) WHEN N.EQ.2.
  220 B(1) = 0.0
      C(1) = 0.0
      D(1) = 0.0
      RETURN
      END

end module
