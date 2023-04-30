SUBROUTINE TRFIND(NST, P, N, X, Y, Z, LIST, LPTR, LEND, B1, B2, B3, I1, I2, I3)
    INTEGER NST, N, LIST(*), LPTR(*), LEND(N), I1, I2, I3
    REAL P(3), X(N), Y(N), Z(N), B1, B2, B3
!
!***********************************************************
!   This subroutine locates a point P relative to a triangu-
! lation created by Subroutine TRMESH.  If P is contained in
! a triangle, the three vertex indexes and barycentric coor-
! dinates are returned.  Otherwise, the indexes of the
! visible boundary nodes are returned.
!
!
! On input:
!
!       NST = Index of a node at which TRFIND begins its
!             search.  Search time depends on the proximity
!             of this node to P.
!
!       P = Array of length 3 containing the x, y, and z
!           coordinates (in that order) of the point P to be
!           located.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y,Z = Arrays of length N containing the Cartesian
!               coordinates of the triangulation nodes (unit
!               vectors).  (X(I),Y(I),Z(I)) defines node I
!               for I = 1 to N.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to Subroutine
!                        TRMESH.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       B1,B2,B3 = Unnormalized barycentric coordinates of
!                  the central projection of P onto the un-
!                  derlying planar triangle if P is in the
!                  convex hull of the nodes.  These parame-
!                  ters are not altered if I1 = 0.
!
!       I1,I2,I3 = Counterclockwise-ordered vertex indexes
!                  of a triangle containing P if P is con-
!                  tained in a triangle.  If P is not in the
!                  convex hull of the nodes, I1 and I2 are
!                  the rightmost and leftmost (boundary)
!                  nodes that are visible from P, and
!                  I3 = 0.  (If all boundary nodes are vis-
!                  ible from P, then I1 and I2 coincide.)
!                  I1 = I2 = I3 = 0 if P and all of the
!                  nodes are coplanar (lie on a common great
!                  circle.
!
!***********************************************************
!
    INTEGER JRAND, LSTPTR
    INTEGER IX, IY, IZ, LP, N0, N1, N1S, N2, N2S, N3, N4, NEXT, NF, NL
    REAL STORE
    REAL DET, EPS, PTN1, PTN2, Q(3), S12, TOL, XP, YP, ZP
    REAL X0, X1, X2, Y0, Y1, Y2, Z0, Z1, Z2
!
    SAVE IX, IY, IZ
    DATA IX/1/, IY/2/, IZ/3/
!
! Local parameters:
!
! EPS =      Machine precision
! IX,IY,IZ = Integer seeds for JRAND
! LP =       LIST pointer
! N0,N1,N2 = Nodes in counterclockwise order defining a
!              cone (with vertex N0) containing P, or end-
!              points of a boundary edge such that P Right
!              N1->N2
! N1S,N2S =  Initially-determined values of N1 and N2
! N3,N4 =    Nodes opposite N1->N2 and N2->N1, respectively
! NEXT =     Candidate for I1 or I2 when P is exterior
! NF,NL =    First and last neighbors of N0, or first
!              (rightmost) and last (leftmost) nodes
!              visible from P when P is exterior to the
!              triangulation
! PTN1 =     Scalar product <P,N1>
! PTN2 =     Scalar product <P,N2>
! Q =        (N2 X N1) X N2  or  N1 X (N2 X N1) -- used in
!              the boundary traversal when P is exterior
! S12 =      Scalar product <N1,N2>
! TOL =      Tolerance (multiple of EPS) defining an upper
!              bound on the magnitude of a negative bary-
!              centric coordinate (B1 or B2) for P in a
!              triangle -- used to avoid an infinite number
!              of restarts with 0 <= B3 < EPS and B1 < 0 or
!              B2 < 0 but small in magnitude
! XP,YP,ZP = Local variables containing P(1), P(2), and P(3)
! X0,Y0,Z0 = Dummy arguments for DET
! X1,Y1,Z1 = Dummy arguments for DET
! X2,Y2,Z2 = Dummy arguments for DET
!
! Statement function:
!
! DET(X1,...,Z0) .GE. 0 if and only if (X0,Y0,Z0) is in the
!                       (closed) left hemisphere defined by
!                       the plane containing (0,0,0),
!                       (X1,Y1,Z1), and (X2,Y2,Z2), where
!                       left is defined relative to an ob-
!                       server at (X1,Y1,Z1) facing
!                       (X2,Y2,Z2).
!
    DET(X1, Y1, Z1, X2, Y2, Z2, X0, Y0, Z0) = X0 * (Y1 * Z2 - Y2 * Z1) - Y0 * (X1 * Z2 - X2 * Z1) + Z0 * (X1 * Y2 - X2 * Y1)
!
! Initialize variables.
!
    XP = P(1)
    YP = P(2)
    ZP = P(3)
    N0 = NST
    IF (N0 .LT. 1 .OR. N0 .GT. N) N0 = JRAND(N, IX, IY, IZ)
!
! Compute the relative machine precision EPS and TOL.
!
    EPS = 1.E0
1   EPS = EPS / 2.E0
    IF (STORE(EPS + 1.E0) .GT. 1.E0) GO TO 1
    EPS = 2.E0 * EPS
    TOL = 100.E0 * EPS
!
! Set NF and NL to the first and last neighbors of N0, and
!   initialize N1 = NF.
!
2   LP = LEND(N0)
    NL = LIST(LP)
    LP = LPTR(LP)
    NF = LIST(LP)
    N1 = NF
!
! Find a pair of adjacent neighbors N1,N2 of N0 that define
!   a wedge containing P:  P LEFT N0->N1 and P RIGHT N0->N2.
!
    IF (NL .GT. 0) THEN
!
!   N0 is an interior node.  Find N1.
!
3       IF (DET(X(N0), Y(N0), Z(N0), X(N1), Y(N1), Z(N1), XP, YP, ZP) .LT. 0.) THEN
            LP = LPTR(LP)
            N1 = LIST(LP)
            IF (N1 .EQ. NL) GO TO 6
            GO TO 3
        END IF
    ELSE
!
!   N0 is a boundary node.  Test for P exterior.
!
        NL = -NL
        IF (DET(X(N0), Y(N0), Z(N0), X(NF), Y(NF), Z(NF), XP, YP, ZP) .LT. 0.) THEN
!
!   P is to the right of the boundary edge N0->NF.
!
            N1 = N0
            N2 = NF
            GO TO 9
        END IF
        IF (DET(X(NL), Y(NL), Z(NL), X(N0), Y(N0), Z(N0), XP, YP, ZP) .LT. 0.) THEN
!
!   P is to the right of the boundary edge NL->N0.
!
            N1 = NL
            N2 = N0
            GO TO 9
        END IF
    END IF
!
! P is to the left of arcs N0->N1 and NL->N0.  Set N2 to the
!   next neighbor of N0 (following N1).
!
4   LP = LPTR(LP)
    N2 = ABS(LIST(LP))
    IF (DET(X(N0), Y(N0), Z(N0), X(N2), Y(N2), Z(N2), XP, YP, ZP) .LT. 0.) GO TO 7
    N1 = N2
    IF (N1 .NE. NL) GO TO 4
    IF (DET(X(N0), Y(N0), Z(N0), X(NF), Y(NF), Z(NF), XP, YP, ZP) .LT. 0.) GO TO 6
!
! P is left of or on arcs N0->NB for all neighbors NB
!   of N0.  Test for P = +/-N0.
!
    IF (STORE(ABS(X(N0) * XP + Y(N0) * YP + Z(N0) * ZP)) .LT. 1.0 - 4.0 * EPS) THEN
!
!   All points are collinear iff P Left NB->N0 for all
!     neighbors NB of N0.  Search the neighbors of N0.
!     Note:  N1 = NL and LP points to NL.
!
5       IF (DET(X(N1), Y(N1), Z(N1), X(N0), Y(N0), Z(N0), XP, YP, ZP) .GE. 0.) THEN
            LP = LPTR(LP)
            N1 = ABS(LIST(LP))
            IF (N1 .EQ. NL) GO TO 14
            GO TO 5
        END IF
    END IF
!
! P is to the right of N1->N0, or P = +/-N0.  Set N0 to N1
!   and start over.
!
    N0 = N1
    GO TO 2
!
! P is between arcs N0->N1 and N0->NF.
!
6   N2 = NF
!
! P is contained in a wedge defined by geodesics N0-N1 and
!   N0-N2, where N1 is adjacent to N2.  Save N1 and N2 to
!   test for cycling.
!
7   N3 = N0
    N1S = N1
    N2S = N2
!
! Top of edge-hopping loop:
!
8   B3 = DET(X(N1), Y(N1), Z(N1), X(N2), Y(N2), Z(N2), XP, YP, ZP)
    IF (B3 .LT. 0.) THEN
!
!   Set N4 to the first neighbor of N2 following N1 (the
!     node opposite N2->N1) unless N1->N2 is a boundary arc.
!
        LP = LSTPTR(LEND(N2), N1, LIST, LPTR)
        IF (LIST(LP) .LT. 0) GO TO 9
        LP = LPTR(LP)
        N4 = ABS(LIST(LP))
!
!   Define a new arc N1->N2 which intersects the geodesic
!     N0-P.
!
        IF (DET(X(N0), Y(N0), Z(N0), X(N4), Y(N4), Z(N4), XP, YP, ZP) .LT. 0.) THEN
            N3 = N2
            N2 = N4
            N1S = N1
            IF (N2 .NE. N2S .AND. N2 .NE. N0) GO TO 8
        ELSE
            N3 = N1
            N1 = N4
            N2S = N2
            IF (N1 .NE. N1S .AND. N1 .NE. N0) GO TO 8
        END IF
!
!   The starting node N0 or edge N1-N2 was encountered
!     again, implying a cycle (infinite loop).  Restart
!     with N0 randomly selected.
!
        N0 = JRAND(N, IX, IY, IZ)
        GO TO 2
    END IF
!
! P is in (N1,N2,N3) unless N0, N1, N2, and P are collinear
!   or P is close to -N0.
!
    IF (B3 .GE. EPS) THEN
!
!   B3 .NE. 0.
!
        B1 = DET(X(N2), Y(N2), Z(N2), X(N3), Y(N3), Z(N3), XP, YP, ZP)
        B2 = DET(X(N3), Y(N3), Z(N3), X(N1), Y(N1), Z(N1), XP, YP, ZP)
        IF (B1 .LT. -TOL .OR. B2 .LT. -TOL) THEN
!
!   Restart with N0 randomly selected.
!
            N0 = JRAND(N, IX, IY, IZ)
            GO TO 2
        END IF
    ELSE
!
!   B3 = 0 and thus P lies on N1->N2. Compute
!     B1 = Det(P,N2 X N1,N2) and B2 = Det(P,N1,N2 X N1).
!
        B3 = 0.
        S12 = X(N1) * X(N2) + Y(N1) * Y(N2) + Z(N1) * Z(N2)
        PTN1 = XP * X(N1) + YP * Y(N1) + ZP * Z(N1)
        PTN2 = XP * X(N2) + YP * Y(N2) + ZP * Z(N2)
        B1 = PTN1 - S12 * PTN2
        B2 = PTN2 - S12 * PTN1
        IF (B1 .LT. -TOL .OR. B2 .LT. -TOL) THEN
!
!   Restart with N0 randomly selected.
!
            N0 = JRAND(N, IX, IY, IZ)
            GO TO 2
        END IF
    END IF
!
! P is in (N1,N2,N3).
!
    I1 = N1
    I2 = N2
    I3 = N3
    IF (B1 .LT. 0.0) B1 = 0.0
    IF (B2 .LT. 0.0) B2 = 0.0
    RETURN
!
! P Right N1->N2, where N1->N2 is a boundary edge.
!   Save N1 and N2, and set NL = 0 to indicate that
!   NL has not yet been found.
!
9   N1S = N1
    N2S = N2
    NL = 0
!
!           Counterclockwise Boundary Traversal:
!
10  LP = LEND(N2)
    LP = LPTR(LP)
    NEXT = LIST(LP)
    IF (DET(X(N2), Y(N2), Z(N2), X(NEXT), Y(NEXT), Z(NEXT), XP, YP, ZP) .GE. 0.) THEN
!
!   N2 is the rightmost visible node if P Forward N2->N1
!     or NEXT Forward N2->N1.  Set Q to (N2 X N1) X N2.
!
        S12 = X(N1) * X(N2) + Y(N1) * Y(N2) + Z(N1) * Z(N2)
        Q(1) = X(N1) - S12 * X(N2)
        Q(2) = Y(N1) - S12 * Y(N2)
        Q(3) = Z(N1) - S12 * Z(N2)
        IF (XP * Q(1) + YP * Q(2) + ZP * Q(3) .GE. 0.) GO TO 11
        IF (X(NEXT) * Q(1) + Y(NEXT) * Q(2) + Z(NEXT) * Q(3) .GE. 0.) GO TO 11
!
!   N1, N2, NEXT, and P are nearly collinear, and N2 is
!     the leftmost visible node.
!
        NL = N2
    END IF
!
! Bottom of counterclockwise loop:
!
    N1 = N2
    N2 = NEXT
    IF (N2 .NE. N1S) GO TO 10
!
! All boundary nodes are visible from P.
!
    I1 = N1S
    I2 = N1S
    I3 = 0
    RETURN
!
! N2 is the rightmost visible node.
!
11  NF = N2
    IF (NL .EQ. 0) THEN
!
! Restore initial values of N1 and N2, and begin the search
!   for the leftmost visible node.
!
        N2 = N2S
        N1 = N1S
!
!           Clockwise Boundary Traversal:
!
12      LP = LEND(N1)
        NEXT = -LIST(LP)
        IF (DET(X(NEXT), Y(NEXT), Z(NEXT), X(N1), Y(N1), Z(N1), XP, YP, ZP) .GE. 0.) THEN
!
!   N1 is the leftmost visible node if P or NEXT is
!     forward of N1->N2.  Compute Q = N1 X (N2 X N1).
!
            S12 = X(N1) * X(N2) + Y(N1) * Y(N2) + Z(N1) * Z(N2)
            Q(1) = X(N2) - S12 * X(N1)
            Q(2) = Y(N2) - S12 * Y(N1)
            Q(3) = Z(N2) - S12 * Z(N1)
            IF (XP * Q(1) + YP * Q(2) + ZP * Q(3) .GE. 0.) GO TO 13
            IF (X(NEXT) * Q(1) + Y(NEXT) * Q(2) + Z(NEXT) * Q(3) .GE. 0.) GO TO 13
!
!   P, NEXT, N1, and N2 are nearly collinear and N1 is the
!     rightmost visible node.
!
            NF = N1
        END IF
!
! Bottom of clockwise loop:
!
        N2 = N1
        N1 = NEXT
        IF (N1 .NE. N1S) GO TO 12
!
! All boundary nodes are visible from P.
!
        I1 = N1
        I2 = N1
        I3 = 0
        RETURN
!
! N1 is the leftmost visible node.
!
13      NL = N1
    END IF
!
! NF and NL have been found.
!
    I1 = NF
    I2 = NL
    I3 = 0
    RETURN
!
! All points are collinear (coplanar).
!
14  I1 = 0
    I2 = 0
    I3 = 0
    RETURN
END
