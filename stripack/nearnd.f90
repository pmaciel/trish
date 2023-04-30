INTEGER FUNCTION NEARND(P, IST, N, X, Y, Z, LIST, LPTR, LEND, AL)
    INTEGER IST, N, LIST(*), LPTR(*), LEND(N)
    REAL P(3), X(N), Y(N), Z(N), AL
!
!***********************************************************
!   Given a point P on the surface of the unit sphere and a
! Delaunay triangulation created by Subroutine TRMESH, this
! function returns the index of the nearest triangulation
! node to P.
!
!   The algorithm consists of implicitly adding P to the
! triangulation, finding the nearest neighbor to P, and
! implicitly deleting P from the triangulation.  Thus, it
! is based on the fact that, if P is a node in a Delaunay
! triangulation, the nearest node to P is a neighbor of P.
!
!
! On input:
!
!       P = Array of length 3 containing the Cartesian coor-
!           dinates of the point P to be located relative to
!           the triangulation.  It is assumed without a test
!           that P(1)**2 + P(2)**2 + P(3)**2 = 1.
!
!       IST = Index of a node at which TRFIND begins the
!             search.  Search time depends on the proximity
!             of this node to P.
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y,Z = Arrays of length N containing the Cartesian
!               coordinates of the nodes.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to TRMESH.
!
! Input parameters are not altered by this function.
!
! On output:
!
!       NEARND = Nodal index of the nearest node to P, or 0
!                if N < 3 or the triangulation data struc-
!                ture is invalid.
!
!       AL = Arc length (angular distance in radians) be-
!            tween P and NEARND unless NEARND = 0.
!
!       Note that the number of candidates for NEARND
!       (neighbors of P) is limited to LMAX defined in
!       the PARAMETER statement below.
!
!***********************************************************
!
    INTEGER LSTPTR
    INTEGER LMAX
    PARAMETER(LMAX=25)
    INTEGER I1, I2, I3, L, LISTP(LMAX), LP, LP1, LP2, LPL, LPTRP(LMAX), N1, N2, N3, NN, NR, NST
    REAL B1, B2, B3, DS1, DSR, DX1, DX2, DX3, DY1, DY2, DY3, DZ1, DZ2, DZ3
!
! Local parameters:
!
! B1,B2,B3 =  Unnormalized barycentric coordinates returned
!               by TRFIND
! DS1 =       (Negative cosine of the) distance from P to N1
! DSR =       (Negative cosine of the) distance from P to NR
! DX1,..DZ3 = Components of vectors used by the swap test
! I1,I2,I3 =  Nodal indexes of a triangle containing P, or
!               the rightmost (I1) and leftmost (I2) visible
!               boundary nodes as viewed from P
! L =         Length of LISTP/LPTRP and number of neighbors
!               of P
! LMAX =      Maximum value of L
! LISTP =     Indexes of the neighbors of P
! LPTRP =     Array of pointers in 1-1 correspondence with
!               LISTP elements
! LP =        LIST pointer to a neighbor of N1 and LISTP
!               pointer
! LP1,LP2 =   LISTP indexes (pointers)
! LPL =       Pointer to the last neighbor of N1
! N1 =        Index of a node visible from P
! N2 =        Index of an endpoint of an arc opposite P
! N3 =        Index of the node opposite N1->N2
! NN =        Local copy of N
! NR =        Index of a candidate for the nearest node to P
! NST =       Index of the node at which TRFIND begins the
!               search
!
!
! Store local parameters and test for N invalid.
!
    NN = N
    IF (NN .LT. 3) GO TO 6
    NST = IST
    IF (NST .LT. 1 .OR. NST .GT. NN) NST = 1
!
! Find a triangle (I1,I2,I3) containing P, or the rightmost
!   (I1) and leftmost (I2) visible boundary nodes as viewed
!   from P.
!
    CALL TRFIND(NST, P, N, X, Y, Z, LIST, LPTR, LEND, B1, B2, B3, I1, I2, I3)
!
! Test for collinear nodes.
!
    IF (I1 .EQ. 0) GO TO 6
!
! Store the linked list of 'neighbors' of P in LISTP and
!   LPTRP.  I1 is the first neighbor, and 0 is stored as
!   the last neighbor if P is not contained in a triangle.
!   L is the length of LISTP and LPTRP, and is limited to
!   LMAX.
!
    IF (I3 .NE. 0) THEN
        LISTP(1) = I1
        LPTRP(1) = 2
        LISTP(2) = I2
        LPTRP(2) = 3
        LISTP(3) = I3
        LPTRP(3) = 1
        L = 3
    ELSE
        N1 = I1
        L = 1
        LP1 = 2
        LISTP(L) = N1
        LPTRP(L) = LP1
!
!   Loop on the ordered sequence of visible boundary nodes
!     N1 from I1 to I2.
!
1       LPL = LEND(N1)
        N1 = -LIST(LPL)
        L = LP1
        LP1 = L + 1
        LISTP(L) = N1
        LPTRP(L) = LP1
        IF (N1 .NE. I2 .AND. LP1 .LT. LMAX) GO TO 1
        L = LP1
        LISTP(L) = 0
        LPTRP(L) = 1
    END IF
!
! Initialize variables for a loop on arcs N1-N2 opposite P
!   in which new 'neighbors' are 'swapped' in.  N1 follows
!   N2 as a neighbor of P, and LP1 and LP2 are the LISTP
!   indexes of N1 and N2.
!
    LP2 = 1
    N2 = I1
    LP1 = LPTRP(1)
    N1 = LISTP(LP1)
!
! Begin loop:  find the node N3 opposite N1->N2.
!
2   LP = LSTPTR(LEND(N1), N2, LIST, LPTR)
    IF (LIST(LP) .LT. 0) GO TO 3
    LP = LPTR(LP)
    N3 = ABS(LIST(LP))
!
! Swap test:  Exit the loop if L = LMAX.
!
    IF (L .EQ. LMAX) GO TO 4
    DX1 = X(N1) - P(1)
    DY1 = Y(N1) - P(2)
    DZ1 = Z(N1) - P(3)
!
    DX2 = X(N2) - P(1)
    DY2 = Y(N2) - P(2)
    DZ2 = Z(N2) - P(3)
!
    DX3 = X(N3) - P(1)
    DY3 = Y(N3) - P(2)
    DZ3 = Z(N3) - P(3)
    IF (DX3 * (DY2 * DZ1 - DY1 * DZ2) - DY3 * (DX2 * DZ1 - DX1 * DZ2) + DZ3 * (DX2 * DY1 - DX1 * DY2) .LE. 0.) GO TO 3
!
! Swap:  Insert N3 following N2 in the adjacency list for P.
!        The two new arcs opposite P must be tested.
!
    L = L + 1
    LPTRP(LP2) = L
    LISTP(L) = N3
    LPTRP(L) = LP1
    LP1 = L
    N1 = N3
    GO TO 2
!
! No swap:  Advance to the next arc and test for termination
!           on N1 = I1 (LP1 = 1) or N1 followed by 0.
!
3   IF (LP1 .EQ. 1) GO TO 4
    LP2 = LP1
    N2 = N1
    LP1 = LPTRP(LP1)
    N1 = LISTP(LP1)
    IF (N1 .EQ. 0) GO TO 4
    GO TO 2
!
! Set NR and DSR to the index of the nearest node to P and
!   an increasing function (negative cosine) of its distance
!   from P, respectively.
!
4   NR = I1
    DSR = -(X(NR) * P(1) + Y(NR) * P(2) + Z(NR) * P(3))
    DO 5 LP = 2, L
        N1 = LISTP(LP)
        IF (N1 .EQ. 0) GO TO 5
        DS1 = -(X(N1) * P(1) + Y(N1) * P(2) + Z(N1) * P(3))
        IF (DS1 .LT. DSR) THEN
            NR = N1
            DSR = DS1
        END IF
5       CONTINUE
        DSR = -DSR
        IF (DSR .GT. 1.0) DSR = 1.0
        AL = ACOS(DSR)
        NEARND = NR
        RETURN
!
! Invalid input.
!
6       NEARND = 0
        RETURN
    END
