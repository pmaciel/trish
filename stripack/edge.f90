SUBROUTINE EDGE(IN1, IN2, X, Y, Z, LWK, IWK, LIST, LPTR, LEND, IER)
    INTEGER IN1, IN2, LWK, IWK(2, *), LIST(*), LPTR(*), LEND(*), IER
    REAL X(*), Y(*), Z(*)
!
!***********************************************************
!   Given a triangulation of N nodes and a pair of nodal
! indexes IN1 and IN2, this routine swaps arcs as necessary
! to force IN1 and IN2 to be adjacent.  Only arcs which
! intersect IN1-IN2 are swapped out.  If a Delaunay triangu-
! lation is input, the resulting triangulation is as close
! as possible to a Delaunay triangulation in the sense that
! all arcs other than IN1-IN2 are locally optimal.
!
!   A sequence of calls to EDGE may be used to force the
! presence of a set of edges defining the boundary of a non-
! convex and/or multiply connected region, or to introduce
! barriers into the triangulation.  Note that Subroutine
! GETNP will not necessarily return closest nodes if the
! triangulation has been constrained by a call to EDGE.
! However, this is appropriate in some applications, such
! as triangle-based interpolation on a nonconvex domain.
!
!
! On input:
!
!       IN1,IN2 = Indexes (of X, Y, and Z) in the range 1 to
!                 N defining a pair of nodes to be connected
!                 by an arc.
!
!       X,Y,Z = Arrays of length N containing the Cartesian
!               coordinates of the nodes.
!
! The above parameters are not altered by this routine.
!
!       LWK = Number of columns reserved for IWK.  This must
!             be at least NI -- the number of arcs that
!             intersect IN1-IN2.  (NI is bounded by N-3.)
!
!       IWK = Integer work array of length at least 2*LWK.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to Subroutine
!                        TRMESH.
!
! On output:
!
!       LWK = Number of arcs which intersect IN1-IN2 (but
!             not more than the input value of LWK) unless
!             IER = 1 or IER = 3.  LWK = 0 if and only if
!             IN1 and IN2 were adjacent (or LWK=0) on input.
!
!       IWK = Array containing the indexes of the endpoints
!             of the new arcs other than IN1-IN2 unless
!             IER > 0 or LWK = 0.  New arcs to the left of
!             IN1->IN2 are stored in the first K-1 columns
!             (left portion of IWK), column K contains
!             zeros, and new arcs to the right of IN1->IN2
!             occupy columns K+1,...,LWK.  (K can be deter-
!             mined by searching IWK for the zeros.)
!
!       LIST,LPTR,LEND = Data structure updated if necessary
!                        to reflect the presence of an arc
!                        connecting IN1 and IN2 unless IER >
!                        0.  The data structure has been
!                        altered if IER >= 4.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if IN1 < 1, IN2 < 1, IN1 = IN2,
!                     or LWK < 0 on input.
!             IER = 2 if more space is required in IWK.
!                     Refer to LWK.
!             IER = 3 if IN1 and IN2 could not be connected
!                     due to either an invalid data struc-
!                     ture or collinear nodes (and floating
!                     point error).
!             IER = 4 if an error flag other than IER = 1
!                     was returned by OPTIM.
!             IER = 5 if error flag 1 was returned by OPTIM.
!                     This is not necessarily an error, but
!                     the arcs other than IN1-IN2 may not
!                     be optimal.
!
!   An error message is written to the standard output unit
! in the case of IER = 3 or IER = 4.
!
!***********************************************************
!
    LOGICAL LEFT
    INTEGER I, IERR, IWC, IWCP1, IWEND, IWF, IWL, LFT, LP, LP21, LPL, N0, N1, N1FRST, N1LST, N2, NEXT, NIT, NL, NR
    REAL DP12, DP1L, DP1R, DP2L, DP2R, X0, X1, X2, Y0, Y1, Y2, Z0, Z1, Z2
!
! Local parameters:
!
! DPij =     Dot product <Ni,Nj>
! I =        DO-loop index and column index for IWK
! IERR =     Error flag returned by Subroutine OPTIM
! IWC =      IWK index between IWF and IWL -- NL->NR is
!              stored in IWK(1,IWC)->IWK(2,IWC)
! IWCP1 =    IWC + 1
! IWEND =    Input or output value of LWK
! IWF =      IWK (column) index of the first (leftmost) arc
!              which intersects IN1->IN2
! IWL =      IWK (column) index of the last (rightmost) are
!              which intersects IN1->IN2
! LFT =      Flag used to determine if a swap results in the
!              new arc intersecting IN1-IN2 -- LFT = 0 iff
!              N0 = IN1, LFT = -1 implies N0 LEFT IN1->IN2,
!              and LFT = 1 implies N0 LEFT IN2->IN1
! LP =       List pointer (index for LIST and LPTR)
! LP21 =     Unused parameter returned by SWAP
! LPL =      Pointer to the last neighbor of IN1 or NL
! N0 =       Neighbor of N1 or node opposite NR->NL
! N1,N2 =    Local copies of IN1 and IN2
! N1FRST =   First neighbor of IN1
! N1LST =    (Signed) last neighbor of IN1
! NEXT =     Node opposite NL->NR
! NIT =      Flag or number of iterations employed by OPTIM
! NL,NR =    Endpoints of an arc which intersects IN1-IN2
!              with NL LEFT IN1->IN2
! X0,Y0,Z0 = Coordinates of N0
! X1,Y1,Z1 = Coordinates of IN1
! X2,Y2,Z2 = Coordinates of IN2
!
!
! Store IN1, IN2, and LWK in local variables and test for
!   errors.
!
    N1 = IN1
    N2 = IN2
    IWEND = LWK
    IF (N1 .LT. 1 .OR. N2 .LT. 1 .OR. N1 .EQ. N2 .OR. IWEND .LT. 0) GO TO 31
!
! Test for N2 as a neighbor of N1.  LPL points to the last
!   neighbor of N1.
!
    LPL = LEND(N1)
    N0 = ABS(LIST(LPL))
    LP = LPL
1   IF (N0 .EQ. N2) GO TO 30
    LP = LPTR(LP)
    N0 = LIST(LP)
    IF (LP .NE. LPL) GO TO 1
!
! Initialize parameters.
!
    IWL = 0
    NIT = 0
!
! Store the coordinates of N1 and N2.
!
2   X1 = X(N1)
    Y1 = Y(N1)
    Z1 = Z(N1)
    X2 = X(N2)
    Y2 = Y(N2)
    Z2 = Z(N2)
!
! Set NR and NL to adjacent neighbors of N1 such that
!   NR LEFT N2->N1 and NL LEFT N1->N2,
!   (NR Forward N1->N2 or NL Forward N1->N2), and
!   (NR Forward N2->N1 or NL Forward N2->N1).
!
!   Initialization:  Set N1FRST and N1LST to the first and
!     (signed) last neighbors of N1, respectively, and
!     initialize NL to N1FRST.
!
    LPL = LEND(N1)
    N1LST = LIST(LPL)
    LP = LPTR(LPL)
    N1FRST = LIST(LP)
    NL = N1FRST
    IF (N1LST .LT. 0) GO TO 4
!
!   N1 is an interior node.  Set NL to the first candidate
!     for NR (NL LEFT N2->N1).
!
3   IF (LEFT(X2, Y2, Z2, X1, Y1, Z1, X(NL), Y(NL), Z(NL))) GO TO 4
    LP = LPTR(LP)
    NL = LIST(LP)
    IF (NL .NE. N1FRST) GO TO 3
!
!   All neighbors of N1 are strictly left of N1->N2.
!
    GO TO 5
!
!   NL = LIST(LP) LEFT N2->N1.  Set NR to NL and NL to the
!     following neighbor of N1.
!
4   NR = NL
    LP = LPTR(LP)
    NL = ABS(LIST(LP))
    IF (LEFT(X1, Y1, Z1, X2, Y2, Z2, X(NL), Y(NL), Z(NL))) THEN
!
!   NL LEFT N1->N2 and NR LEFT N2->N1.  The Forward tests
!     are employed to avoid an error associated with
!     collinear nodes.
!
        DP12 = X1 * X2 + Y1 * Y2 + Z1 * Z2
        DP1L = X1 * X(NL) + Y1 * Y(NL) + Z1 * Z(NL)
        DP2L = X2 * X(NL) + Y2 * Y(NL) + Z2 * Z(NL)
        DP1R = X1 * X(NR) + Y1 * Y(NR) + Z1 * Z(NR)
        DP2R = X2 * X(NR) + Y2 * Y(NR) + Z2 * Z(NR)
        IF ((DP2L - DP12 * DP1L .GE. 0. .OR. DP2R - DP12 * DP1R .GE. 0.) &
            .AND. (DP1L - DP12 * DP2L .GE. 0. .OR. DP1R - DP12 * DP2R .GE. 0.)) GO TO 6
!
!   NL-NR does not intersect N1-N2.  However, there is
!     another candidate for the first arc if NL lies on
!     the line N1-N2.
!
        IF (.NOT. LEFT(X2, Y2, Z2, X1, Y1, Z1, X(NL), Y(NL), Z(NL))) GO TO 5
    END IF
!
!   Bottom of loop.
!
    IF (NL .NE. N1FRST) GO TO 4
!
! Either the triangulation is invalid or N1-N2 lies on the
!   convex hull boundary and an edge NR->NL (opposite N1 and
!   intersecting N1-N2) was not found due to floating point
!   error.  Try interchanging N1 and N2 -- NIT > 0 iff this
!   has already been done.
!
5   IF (NIT .GT. 0) GO TO 33
    NIT = 1
    N1 = N2
    N2 = IN1
    GO TO 2
!
! Store the ordered sequence of intersecting edges NL->NR in
!   IWK(1,IWL)->IWK(2,IWL).
!
6   IWL = IWL + 1
    IF (IWL .GT. IWEND) GO TO 32
    IWK(1, IWL) = NL
    IWK(2, IWL) = NR
!
!   Set NEXT to the neighbor of NL which follows NR.
!
    LPL = LEND(NL)
    LP = LPTR(LPL)
!
!   Find NR as a neighbor of NL.  The search begins with
!     the first neighbor.
!
7   IF (LIST(LP) .EQ. NR) GO TO 8
    LP = LPTR(LP)
    IF (LP .NE. LPL) GO TO 7
!
!   NR must be the last neighbor, and NL->NR cannot be a
!     boundary edge.
!
    IF (LIST(LP) .NE. NR) GO TO 33
!
!   Set NEXT to the neighbor following NR, and test for
!     termination of the store loop.
!
8   LP = LPTR(LP)
    NEXT = ABS(LIST(LP))
    IF (NEXT .EQ. N2) GO TO 9
!
!   Set NL or NR to NEXT.
!
    IF (LEFT(X1, Y1, Z1, X2, Y2, Z2, X(NEXT), Y(NEXT), Z(NEXT))) THEN
        NL = NEXT
    ELSE
        NR = NEXT
    END IF
    GO TO 6
!
! IWL is the number of arcs which intersect N1-N2.
!   Store LWK.
!
9   LWK = IWL
    IWEND = IWL
!
! Initialize for edge swapping loop -- all possible swaps
!   are applied (even if the new arc again intersects
!   N1-N2), arcs to the left of N1->N2 are stored in the
!   left portion of IWK, and arcs to the right are stored in
!   the right portion.  IWF and IWL index the first and last
!   intersecting arcs.
!
    IWF = 1
!
! Top of loop -- set N0 to N1 and NL->NR to the first edge.
!   IWC points to the arc currently being processed.  LFT
!   .LE. 0 iff N0 LEFT N1->N2.
!
10  LFT = 0
    N0 = N1
    X0 = X1
    Y0 = Y1
    Z0 = Z1
    NL = IWK(1, IWF)
    NR = IWK(2, IWF)
    IWC = IWF
!
!   Set NEXT to the node opposite NL->NR unless IWC is the
!     last arc.
!
11  IF (IWC .EQ. IWL) GO TO 21
    IWCP1 = IWC + 1
    NEXT = IWK(1, IWCP1)
    IF (NEXT .NE. NL) GO TO 16
    NEXT = IWK(2, IWCP1)
!
!   NEXT RIGHT N1->N2 and IWC .LT. IWL.  Test for a possible
!     swap.
!
    IF (.NOT. LEFT(X0, Y0, Z0, X(NR), Y(NR), Z(NR), X(NEXT), Y(NEXT), Z(NEXT))) GO TO 14
    IF (LFT .GE. 0) GO TO 12
    IF (.NOT. LEFT(X(NL), Y(NL), Z(NL), X0, Y0, Z0, X(NEXT), Y(NEXT), Z(NEXT))) GO TO 14
!
!   Replace NL->NR with N0->NEXT.
!
    CALL SWAP(NEXT, N0, NL, NR, LIST, LPTR, LEND, LP21)
    IWK(1, IWC) = N0
    IWK(2, IWC) = NEXT
    GO TO 15
!
!   Swap NL-NR for N0-NEXT, shift columns IWC+1,...,IWL to
!     the left, and store N0-NEXT in the right portion of
!     IWK.
!
12  CALL SWAP(NEXT, N0, NL, NR, LIST, LPTR, LEND, LP21)
    DO 13 I = IWCP1, IWL
        IWK(1, I - 1) = IWK(1, I)
        IWK(2, I - 1) = IWK(2, I)
13      CONTINUE
        IWK(1, IWL) = N0
        IWK(2, IWL) = NEXT
        IWL = IWL - 1
        NR = NEXT
        GO TO 11
!
!   A swap is not possible.  Set N0 to NR.
!
14      N0 = NR
        X0 = X(N0)
        Y0 = Y(N0)
        Z0 = Z(N0)
        LFT = 1
!
!   Advance to the next arc.
!
15      NR = NEXT
        IWC = IWC + 1
        GO TO 11
!
!   NEXT LEFT N1->N2, NEXT .NE. N2, and IWC .LT. IWL.
!     Test for a possible swap.
!
16      IF (.NOT. LEFT(X(NL), Y(NL), Z(NL), X0, Y0, Z0, X(NEXT), Y(NEXT), Z(NEXT))) GO TO 19
        IF (LFT .LE. 0) GO TO 17
        IF (.NOT. LEFT(X0, Y0, Z0, X(NR), Y(NR), Z(NR), X(NEXT), Y(NEXT), Z(NEXT))) GO TO 19
!
!   Replace NL->NR with NEXT->N0.
!
        CALL SWAP(NEXT, N0, NL, NR, LIST, LPTR, LEND, LP21)
        IWK(1, IWC) = NEXT
        IWK(2, IWC) = N0
        GO TO 20
!
!   Swap NL-NR for N0-NEXT, shift columns IWF,...,IWC-1 to
!     the right, and store N0-NEXT in the left portion of
!     IWK.
!
17      CALL SWAP(NEXT, N0, NL, NR, LIST, LPTR, LEND, LP21)
        DO 18 I = IWC - 1, IWF, -1
            IWK(1, I + 1) = IWK(1, I)
            IWK(2, I + 1) = IWK(2, I)
18          CONTINUE
            IWK(1, IWF) = N0
            IWK(2, IWF) = NEXT
            IWF = IWF + 1
            GO TO 20
!
!   A swap is not possible.  Set N0 to NL.
!
19          N0 = NL
            X0 = X(N0)
            Y0 = Y(N0)
            Z0 = Z(N0)
            LFT = -1
!
!   Advance to the next arc.
!
20          NL = NEXT
            IWC = IWC + 1
            GO TO 11
!
!   N2 is opposite NL->NR (IWC = IWL).
!
21          IF (N0 .EQ. N1) GO TO 24
            IF (LFT .LT. 0) GO TO 22
!
!   N0 RIGHT N1->N2.  Test for a possible swap.
!
            IF (.NOT. LEFT(X0, Y0, Z0, X(NR), Y(NR), Z(NR), X2, Y2, Z2)) GO TO 10
!
!   Swap NL-NR for N0-N2 and store N0-N2 in the right
!     portion of IWK.
!
            CALL SWAP(N2, N0, NL, NR, LIST, LPTR, LEND, LP21)
            IWK(1, IWL) = N0
            IWK(2, IWL) = N2
            IWL = IWL - 1
            GO TO 10
!
!   N0 LEFT N1->N2.  Test for a possible swap.
!
22          IF (.NOT. LEFT(X(NL), Y(NL), Z(NL), X0, Y0, Z0, X2, Y2, Z2)) GO TO 10
!
!   Swap NL-NR for N0-N2, shift columns IWF,...,IWL-1 to the
!     right, and store N0-N2 in the left portion of IWK.
!
            CALL SWAP(N2, N0, NL, NR, LIST, LPTR, LEND, LP21)
            I = IWL
23          IWK(1, I) = IWK(1, I - 1)
            IWK(2, I) = IWK(2, I - 1)
            I = I - 1
            IF (I .GT. IWF) GO TO 23
            IWK(1, IWF) = N0
            IWK(2, IWF) = N2
            IWF = IWF + 1
            GO TO 10
!
! IWF = IWC = IWL.  Swap out the last arc for N1-N2 and
!   store zeros in IWK.
!
24          CALL SWAP(N2, N1, NL, NR, LIST, LPTR, LEND, LP21)
            IWK(1, IWC) = 0
            IWK(2, IWC) = 0
!
! Optimization procedure --
!
            IER = 0
            IF (IWC .GT. 1) THEN
!
!   Optimize the set of new arcs to the left of IN1->IN2.
!
                NIT = 4 * (IWC - 1)
                CALL OPTIM(X, Y, Z, IWC - 1, LIST, LPTR, LEND, NIT, IWK, IERR)
                IF (IERR .NE. 0 .AND. IERR .NE. 1) GO TO 34
                IF (IERR .EQ. 1) IER = 5
            END IF
            IF (IWC .LT. IWEND) THEN
!
!   Optimize the set of new arcs to the right of IN1->IN2.
!
                NIT = 4 * (IWEND - IWC)
                CALL OPTIM(X, Y, Z, IWEND - IWC, LIST, LPTR, LEND, NIT, IWK(1, IWC + 1), IERR)
                IF (IERR .NE. 0 .AND. IERR .NE. 1) GO TO 34
                IF (IERR .EQ. 1) GO TO 35
            END IF
            IF (IER .EQ. 5) GO TO 35
!
! Successful termination (IER = 0).
!
            RETURN
!
! IN1 and IN2 were adjacent on input.
!
30          IER = 0
            RETURN
!
! Invalid input parameter.
!
31          IER = 1
            RETURN
!
! Insufficient space reserved for IWK.
!
32          IER = 2
            RETURN
!
! Invalid triangulation data structure or collinear nodes
!   on convex hull boundary.
!
33          IER = 3
            WRITE (*, 130) IN1, IN2
130         FORMAT(//5X, '*** Error in EDGE:  Invalid triangula', 'tion or null triangles on boundary'/9X, &
                    'IN1 =', I4, ', IN2=', I4/)
            RETURN
!
! Error flag (other than 1) returned by OPTIM.
!
34          IER = 4
            WRITE (*, 140) NIT, IERR
140         FORMAT(//5X, '*** Error in OPTIM (called from EDGE):', '  NIT = ', I4, ', IER = ', I1, ' ***'/)
            RETURN
!
! Error flag 1 returned by OPTIM.
!
35          IER = 5
            RETURN
        END
