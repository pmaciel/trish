LOGICAL FUNCTION INSIDE(P, LV, XV, YV, ZV, NV, LISTV, IER)
    INTEGER LV, NV, LISTV(NV), IER
    REAL P(3), XV(LV), YV(LV), ZV(LV)
!
!***********************************************************
!   This function locates a point P relative to a polygonal
! region R on the surface of the unit sphere, returning
! INSIDE = TRUE if and only if P is contained in R.  R is
! defined by a cyclically ordered sequence of vertices which
! form a positively-oriented simple closed curve.  Adjacent
! vertices need not be distinct but the curve must not be
! self-intersecting.  Also, while polygon edges are by defi-
! nition restricted to a single hemisphere, R is not so
! restricted.  Its interior is the region to the left as the
! vertices are traversed in order.
!
!   The algorithm consists of selecting a point Q in R and
! then finding all points at which the great circle defined
! by P and Q intersects the boundary of R.  P lies inside R
! if and only if there is an even number of intersection
! points between Q and P.  Q is taken to be a point immedi-
! ately to the left of a directed boundary edge -- the first
! one that results in no consistency-check failures.
!
!   If P is close to the polygon boundary, the problem is
! ill-conditioned and the decision may be incorrect.  Also,
! an incorrect decision may result from a poor choice of Q
! (if, for example, a boundary edge lies on the great cir-
! cle defined by P and Q).  A more reliable result could be
! obtained by a sequence of calls to INSIDE with the ver-
! tices cyclically permuted before each call (to alter the
! choice of Q).
!
!
! On input:
!
!       P = Array of length 3 containing the Cartesian
!           coordinates of the point (unit vector) to be
!           located.
!
!       LV = Length of arrays XV, YV, and ZV.
!
!       XV,YV,ZV = Arrays of length LV containing the Carte-
!                  sian coordinates of unit vectors (points
!                  on the unit sphere).  These values are
!                  not tested for validity.
!
!       NV = Number of vertices in the polygon.  3 .LE. NV
!            .LE. LV.
!
!       LISTV = Array of length NV containing the indexes
!               (for XV, YV, and ZV) of a cyclically-ordered
!               (and CCW-ordered) sequence of vertices that
!               define R.  The last vertex (indexed by
!               LISTV(NV)) is followed by the first (indexed
!               by LISTV(1)).  LISTV entries must be in the
!               range 1 to LV.
!
! Input parameters are not altered by this function.
!
! On output:
!
!       INSIDE = TRUE if and only if P lies inside R unless
!                IER .NE. 0, in which case the value is not
!                altered.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if LV or NV is outside its valid
!                     range.
!             IER = 2 if a LISTV entry is outside its valid
!                     range.
!             IER = 3 if the polygon boundary was found to
!                     be self-intersecting.  This error will
!                     not necessarily be detected.
!             IER = 4 if every choice of Q (one for each
!                     boundary edge) led to failure of some
!                     internal consistency check.  The most
!                     likely cause of this error is invalid
!                     input:  P = (0,0,0), a null or self-
!                     intersecting polygon, etc.
!
!***********************************************************
!
    INTEGER I1, I2, IERR, IMX, K, K0, N, NI
    LOGICAL EVEN, LFT1, LFT2, PINR, QINR
    REAL B(3), BP, BQ, CN(3), D, EPS, PN(3), Q(3), QN(3), QNRM, V1(3), V2(3), VN(3), VNRM
!
! Local parameters:
!
! B =         Intersection point between the boundary and
!               the great circle defined by P and Q
! BP,BQ =     <B,P> and <B,Q>, respectively, maximized over
!               intersection points B that lie between P and
!               Q (on the shorter arc) -- used to find the
!               closest intersection points to P and Q
! CN =        Q X P = normal to the plane of P and Q
! D =         Dot product <B,P> or <B,Q>
! EPS =       Parameter used to define Q as the point whose
!               orthogonal distance to (the midpoint of)
!               boundary edge V1->V2 is approximately EPS/
!               (2*Cos(A/2)), where <V1,V2> = Cos(A).
! EVEN =      TRUE iff an even number of intersection points
!               lie between P and Q (on the shorter arc)
! I1,I2 =     Indexes (LISTV elements) of a pair of adjacent
!               boundary vertices (endpoints of a boundary
!               edge)
! IERR =      Error flag for calls to INTRSC (not tested)
! IMX =       Local copy of LV and maximum value of I1 and
!               I2
! K =         DO-loop index and LISTV index
! K0 =        LISTV index of the first endpoint of the
!               boundary edge used to compute Q
! LFT1,LFT2 = Logical variables associated with I1 and I2 in
!               the boundary traversal:  TRUE iff the vertex
!               is strictly to the left of Q->P (<V,CN> > 0)
! N =         Local copy of NV
! NI =        Number of intersections (between the boundary
!               curve and the great circle P-Q) encountered
! PINR =      TRUE iff P is to the left of the directed
!               boundary edge associated with the closest
!               intersection point to P that lies between P
!               and Q (a left-to-right intersection as
!               viewed from Q), or there is no intersection
!               between P and Q (on the shorter arc)
! PN,QN =     P X CN and CN X Q, respectively:  used to
!               locate intersections B relative to arc Q->P
! Q =         (V1 + V2 + EPS*VN/VNRM)/QNRM, where V1->V2 is
!               the boundary edge indexed by LISTV(K0) ->
!               LISTV(K0+1)
! QINR =      TRUE iff Q is to the left of the directed
!               boundary edge associated with the closest
!               intersection point to Q that lies between P
!               and Q (a right-to-left intersection as
!               viewed from Q), or there is no intersection
!               between P and Q (on the shorter arc)
! QNRM =      Euclidean norm of V1+V2+EPS*VN/VNRM used to
!               compute (normalize) Q
! V1,V2 =     Vertices indexed by I1 and I2 in the boundary
!               traversal
! VN =        V1 X V2, where V1->V2 is the boundary edge
!               indexed by LISTV(K0) -> LISTV(K0+1)
! VNRM =      Euclidean norm of VN
!
    DATA EPS/1.E-3/
!
! Store local parameters, test for error 1, and initialize
!   K0.
!
    IMX = LV
    N = NV
    IF (N .LT. 3 .OR. N .GT. IMX) GO TO 11
    K0 = 0
    I1 = LISTV(1)
    IF (I1 .LT. 1 .OR. I1 .GT. IMX) GO TO 12
!
! Increment K0 and set Q to a point immediately to the left
!   of the midpoint of edge V1->V2 = LISTV(K0)->LISTV(K0+1):
!   Q = (V1 + V2 + EPS*VN/VNRM)/QNRM, where VN = V1 X V2.
!
1   K0 = K0 + 1
    IF (K0 .GT. N) GO TO 14
    I1 = LISTV(K0)
    IF (K0 .LT. N) THEN
        I2 = LISTV(K0 + 1)
    ELSE
        I2 = LISTV(1)
    END IF
    IF (I2 .LT. 1 .OR. I2 .GT. IMX) GO TO 12
    VN(1) = YV(I1) * ZV(I2) - ZV(I1) * YV(I2)
    VN(2) = ZV(I1) * XV(I2) - XV(I1) * ZV(I2)
    VN(3) = XV(I1) * YV(I2) - YV(I1) * XV(I2)
    VNRM = SQRT(VN(1) * VN(1) + VN(2) * VN(2) + VN(3) * VN(3))
    IF (VNRM .EQ. 0.) GO TO 1
    Q(1) = XV(I1) + XV(I2) + EPS * VN(1) / VNRM
    Q(2) = YV(I1) + YV(I2) + EPS * VN(2) / VNRM
    Q(3) = ZV(I1) + ZV(I2) + EPS * VN(3) / VNRM
    QNRM = SQRT(Q(1) * Q(1) + Q(2) * Q(2) + Q(3) * Q(3))
    Q(1) = Q(1) / QNRM
    Q(2) = Q(2) / QNRM
    Q(3) = Q(3) / QNRM
!
! Compute CN = Q X P, PN = P X CN, and QN = CN X Q.
!
    CN(1) = Q(2) * P(3) - Q(3) * P(2)
    CN(2) = Q(3) * P(1) - Q(1) * P(3)
    CN(3) = Q(1) * P(2) - Q(2) * P(1)
    IF (CN(1) .EQ. 0. .AND. CN(2) .EQ. 0. .AND. CN(3) .EQ. 0.) GO TO 1
    PN(1) = P(2) * CN(3) - P(3) * CN(2)
    PN(2) = P(3) * CN(1) - P(1) * CN(3)
    PN(3) = P(1) * CN(2) - P(2) * CN(1)
    QN(1) = CN(2) * Q(3) - CN(3) * Q(2)
    QN(2) = CN(3) * Q(1) - CN(1) * Q(3)
    QN(3) = CN(1) * Q(2) - CN(2) * Q(1)
!
! Initialize parameters for the boundary traversal.
!
    NI = 0
    EVEN = .TRUE.
    BP = -2.
    BQ = -2.
    PINR = .TRUE.
    QINR = .TRUE.
    I2 = LISTV(N)
    IF (I2 .LT. 1 .OR. I2 .GT. IMX) GO TO 12
    LFT2 = CN(1) * XV(I2) + CN(2) * YV(I2) + CN(3) * ZV(I2) .GT. 0.
!
! Loop on boundary arcs I1->I2.
!
    DO 2 K = 1, N
        I1 = I2
        LFT1 = LFT2
        I2 = LISTV(K)
        IF (I2 .LT. 1 .OR. I2 .GT. IMX) GO TO 12
        LFT2 = CN(1) * XV(I2) + CN(2) * YV(I2) + CN(3) * ZV(I2) .GT. 0.
        IF (LFT1 .EQV. LFT2) GO TO 2
!
!   I1 and I2 are on opposite sides of Q->P.  Compute the
!     point of intersection B.
!
        NI = NI + 1
        V1(1) = XV(I1)
        V1(2) = YV(I1)
        V1(3) = ZV(I1)
        V2(1) = XV(I2)
        V2(2) = YV(I2)
        V2(3) = ZV(I2)
        CALL INTRSC(V1, V2, CN, B, IERR)
!
!   B is between Q and P (on the shorter arc) iff
!     B Forward Q->P and B Forward P->Q       iff
!     <B,QN> > 0 and <B,PN> > 0.
!
        IF (B(1) * QN(1) + B(2) * QN(2) + B(3) * QN(3) .GT. 0. .AND. B(1) * PN(1) + B(2) * PN(2) + B(3) * PN(3) .GT. 0.) THEN
!
!   Update EVEN, BQ, QINR, BP, and PINR.
!
            EVEN = .NOT. EVEN
            D = B(1) * Q(1) + B(2) * Q(2) + B(3) * Q(3)
            IF (D .GT. BQ) THEN
                BQ = D
                QINR = LFT2
            END IF
            D = B(1) * P(1) + B(2) * P(2) + B(3) * P(3)
            IF (D .GT. BP) THEN
                BP = D
                PINR = LFT1
            END IF
        END IF
2       CONTINUE
!
! Test for consistency:  NI must be even and QINR must be
!   TRUE.
!
        IF (NI .NE. 2 * (NI / 2) .OR. .NOT. QINR) GO TO 1
!
! Test for error 3:  different values of PINR and EVEN.
!
        IF (PINR .NEQV. EVEN) GO TO 13
!
! No error encountered.
!
        IER = 0
        INSIDE = EVEN
        RETURN
!
! LV or NV is outside its valid range.
!
11      IER = 1
        RETURN
!
! A LISTV entry is outside its valid range.
!
12      IER = 2
        RETURN
!
! The polygon boundary is self-intersecting.
!
13      IER = 3
        RETURN
!
! Consistency tests failed for all values of Q.
!
14      IER = 4
        RETURN
    END
