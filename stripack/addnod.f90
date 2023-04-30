SUBROUTINE ADDNOD(NST, K, X, Y, Z, LIST, LPTR, LEND, LNEW, IER)
    INTEGER NST, K, LIST(*), LPTR(*), LEND(K), LNEW, IER
    REAL X(K), Y(K), Z(K)
!
!***********************************************************
!   This subroutine adds node K to a triangulation of the
! convex hull of nodes 1,...,K-1, producing a triangulation
! of the convex hull of nodes 1,...,K.
!
!   The algorithm consists of the following steps:  node K
! is located relative to the triangulation (TRFIND), its
! index is added to the data structure (INTADD or BDYADD),
! and a sequence of swaps (SWPTST and SWAP) are applied to
! the arcs opposite K so that all arcs incident on node K
! and opposite node K are locally optimal (satisfy the cir-
! cumcircle test).  Thus, if a Delaunay triangulation is
! input, a Delaunay triangulation will result.
!
!
! On input:
!
!       NST = Index of a node at which TRFIND begins its
!             search.  Search time depends on the proximity
!             of this node to K.  If NST < 1, the search is
!             begun at node K-1.
!
!       K = Nodal index (index for X, Y, Z, and LEND) of the
!           new node to be added.  K .GE. 4.
!
!       X,Y,Z = Arrays of length .GE. K containing Car-
!               tesian coordinates of the nodes.
!               (X(I),Y(I),Z(I)) defines node I for
!               I = 1,...,K.
!
! The above parameters are not altered by this routine.
!
!       LIST,LPTR,LEND,LNEW = Data structure associated with
!                             the triangulation of nodes 1
!                             to K-1.  The array lengths are
!                             assumed to be large enough to
!                             add node K.  Refer to Subrou-
!                             tine TRMESH.
!
! On output:
!
!       LIST,LPTR,LEND,LNEW = Data structure updated with
!                             the addition of node K as the
!                             last entry unless IER .NE. 0
!                             and IER .NE. -3, in which case
!                             the arrays are not altered.
!
!       IER = Error indicator:
!             IER =  0 if no errors were encountered.
!             IER = -1 if K is outside its valid range
!                      on input.
!             IER = -2 if all nodes (including K) are col-
!                      linear (lie on a common geodesic).
!             IER =  L if nodes L and K coincide for some
!                      L < K.
!
!***********************************************************
!
    INTEGER LSTPTR
    INTEGER I1, I2, I3, IO1, IO2, IN1, IST, KK, KM1, L, LP, LPF, LPO1, LPO1S
    LOGICAL SWPTST
    REAL B1, B2, B3, P(3)
!
! Local parameters:
!
! B1,B2,B3 = Unnormalized barycentric coordinates returned
!              by TRFIND.
! I1,I2,I3 = Vertex indexes of a triangle containing K
! IN1 =      Vertex opposite K:  first neighbor of IO2
!              that precedes IO1.  IN1,IO1,IO2 are in
!              counterclockwise order.
! IO1,IO2 =  Adjacent neighbors of K defining an arc to
!              be tested for a swap
! IST =      Index of node at which TRFIND begins its search
! KK =       Local copy of K
! KM1 =      K-1
! L =        Vertex index (I1, I2, or I3) returned in IER
!              if node K coincides with a vertex
! LP =       LIST pointer
! LPF =      LIST pointer to the first neighbor of K
! LPO1 =     LIST pointer to IO1
! LPO1S =    Saved value of LPO1
! P =        Cartesian coordinates of node K
!
    KK = K
    IF (KK .LT. 4) GO TO 3
!
! Initialization:
!
    KM1 = KK - 1
    IST = NST
    IF (IST .LT. 1) IST = KM1
    P(1) = X(KK)
    P(2) = Y(KK)
    P(3) = Z(KK)
!
! Find a triangle (I1,I2,I3) containing K or the rightmost
!   (I1) and leftmost (I2) visible boundary nodes as viewed
!   from node K.
!
    CALL TRFIND(IST, P, KM1, X, Y, Z, LIST, LPTR, LEND, B1, B2, B3, I1, I2, I3)
!
!   Test for collinear or duplicate nodes.
!
    IF (I1 .EQ. 0) GO TO 4
    IF (I3 .NE. 0) THEN
        L = I1
        IF (P(1) .EQ. X(L) .AND. P(2) .EQ. Y(L) .AND. P(3) .EQ. Z(L)) GO TO 5
        L = I2
        IF (P(1) .EQ. X(L) .AND. P(2) .EQ. Y(L) .AND. P(3) .EQ. Z(L)) GO TO 5
        L = I3
        IF (P(1) .EQ. X(L) .AND. P(2) .EQ. Y(L) .AND. P(3) .EQ. Z(L)) GO TO 5
        CALL INTADD(KK, I1, I2, I3, LIST, LPTR, LEND, LNEW)
    ELSE
        IF (I1 .NE. I2) THEN
            CALL BDYADD(KK, I1, I2, LIST, LPTR, LEND, LNEW)
        ELSE
            CALL COVSPH(KK, I1, LIST, LPTR, LEND, LNEW)
        END IF
    END IF
    IER = 0
!
! Initialize variables for optimization of the
!   triangulation.
!
    LP = LEND(KK)
    LPF = LPTR(LP)
    IO2 = LIST(LPF)
    LPO1 = LPTR(LPF)
    IO1 = ABS(LIST(LPO1))
!
! Begin loop:  find the node opposite K.
!
1   LP = LSTPTR(LEND(IO1), IO2, LIST, LPTR)
    IF (LIST(LP) .LT. 0) GO TO 2
    LP = LPTR(LP)
    IN1 = ABS(LIST(LP))
!
! Swap test:  if a swap occurs, two new arcs are
!             opposite K and must be tested.
!
    LPO1S = LPO1
    IF (.NOT. SWPTST(IN1, KK, IO1, IO2, X, Y, Z)) GO TO 2
    CALL SWAP(IN1, KK, IO1, IO2, LIST, LPTR, LEND, LPO1)
    IF (LPO1 .EQ. 0) THEN
!
!   A swap is not possible because KK and IN1 are already
!     adjacent.  This error in SWPTST only occurs in the
!     neutral case and when there are nearly duplicate
!     nodes.
!
        LPO1 = LPO1S
        GO TO 2
    END IF
    IO1 = IN1
    GO TO 1
!
! No swap occurred.  Test for termination and reset
!   IO2 and IO1.
!
2   IF (LPO1 .EQ. LPF .OR. LIST(LPO1) .LT. 0) RETURN
    IO2 = IO1
    LPO1 = LPTR(LPO1)
    IO1 = ABS(LIST(LPO1))
    GO TO 1
!
! KK < 4.
!
3   IER = -1
    RETURN
!
! All nodes are collinear.
!
4   IER = -2
    RETURN
!
! Nodes L and K coincide.
!
5   IER = L
    RETURN
END
