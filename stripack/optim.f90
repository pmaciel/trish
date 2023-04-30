SUBROUTINE OPTIM(X, Y, Z, NA, LIST, LPTR, LEND, NIT, IWK, IER)
    INTEGER NA, LIST(*), LPTR(*), LEND(*), NIT, IWK(2, NA), IER
    REAL X(*), Y(*), Z(*)
!
!***********************************************************
!   Given a set of NA triangulation arcs, this subroutine
! optimizes the portion of the triangulation consisting of
! the quadrilaterals (pairs of adjacent triangles) which
! have the arcs as diagonals by applying the circumcircle
! test and appropriate swaps to the arcs.
!
!   An iteration consists of applying the swap test and
! swaps to all NA arcs in the order in which they are
! stored.  The iteration is repeated until no swap occurs
! or NIT iterations have been performed.  The bound on the
! number of iterations may be necessary to prevent an
! infinite loop caused by cycling (reversing the effect of a
! previous swap) due to floating point inaccuracy when four
! or more nodes are nearly cocircular.
!
!
! On input:
!
!       X,Y,Z = Arrays containing the nodal coordinates.
!
!       NA = Number of arcs in the set.  NA .GE. 0.
!
! The above parameters are not altered by this routine.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to Subroutine
!                        TRMESH.
!
!       NIT = Maximum number of iterations to be performed.
!             NIT = 4*NA should be sufficient.  NIT .GE. 1.
!
!       IWK = Integer array dimensioned 2 by NA containing
!             the nodal indexes of the arc endpoints (pairs
!             of endpoints are stored in columns).
!
! On output:
!
!       LIST,LPTR,LEND = Updated triangulation data struc-
!                        ture reflecting the swaps.
!
!       NIT = Number of iterations performed.
!
!       IWK = Endpoint indexes of the new set of arcs
!             reflecting the swaps.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if a swap occurred on the last of
!                     MAXIT iterations, where MAXIT is the
!                     value of NIT on input.  The new set
!                     of arcs is not necessarily optimal
!                     in this case.
!             IER = 2 if NA < 0 or NIT < 1 on input.
!             IER = 3 if IWK(2,I) is not a neighbor of
!                     IWK(1,I) for some I in the range 1
!                     to NA.  A swap may have occurred in
!                     this case.
!             IER = 4 if a zero pointer was returned by
!                     Subroutine SWAP.
!
!***********************************************************
!
    INTEGER I, IO1, IO2, ITER, LP, LP21, LPL, LPP, MAXIT, N1, N2, NNA
    LOGICAL SWPTST
    LOGICAL SWP
!
! Local parameters:
!
! I =       Column index for IWK
! IO1,IO2 = Nodal indexes of the endpoints of an arc in IWK
! ITER =    Iteration count
! LP =      LIST pointer
! LP21 =    Parameter returned by SWAP (not used)
! LPL =     Pointer to the last neighbor of IO1
! LPP =     Pointer to the node preceding IO2 as a neighbor
!             of IO1
! MAXIT =   Input value of NIT
! N1,N2 =   Nodes opposite IO1->IO2 and IO2->IO1,
!             respectively
! NNA =     Local copy of NA
! SWP =     Flag set to TRUE iff a swap occurs in the
!             optimization loop
!
    NNA = NA
    MAXIT = NIT
    IF (NNA .LT. 0 .OR. MAXIT .LT. 1) GO TO 7
!
! Initialize iteration count ITER and test for NA = 0.
!
    ITER = 0
    IF (NNA .EQ. 0) GO TO 5
!
! Top of loop --
!   SWP = TRUE iff a swap occurred in the current iteration.
!
1   IF (ITER .EQ. MAXIT) GO TO 6
    ITER = ITER + 1
    SWP = .FALSE.
!
!   Inner loop on arcs IO1-IO2 --
!
    DO 4 I = 1, NNA
        IO1 = IWK(1, I)
        IO2 = IWK(2, I)
!
!   Set N1 and N2 to the nodes opposite IO1->IO2 and
!     IO2->IO1, respectively.  Determine the following:
!
!     LPL = pointer to the last neighbor of IO1,
!     LP = pointer to IO2 as a neighbor of IO1, and
!     LPP = pointer to the node N2 preceding IO2.
!
        LPL = LEND(IO1)
        LPP = LPL
        LP = LPTR(LPP)
2       IF (LIST(LP) .EQ. IO2) GO TO 3
        LPP = LP
        LP = LPTR(LPP)
        IF (LP .NE. LPL) GO TO 2
!
!   IO2 should be the last neighbor of IO1.  Test for no
!     arc and bypass the swap test if IO1 is a boundary
!     node.
!
        IF (ABS(LIST(LP)) .NE. IO2) GO TO 8
        IF (LIST(LP) .LT. 0) GO TO 4
!
!   Store N1 and N2, or bypass the swap test if IO1 is a
!     boundary node and IO2 is its first neighbor.
!
3       N2 = LIST(LPP)
        IF (N2 .LT. 0) GO TO 4
        LP = LPTR(LP)
        N1 = ABS(LIST(LP))
!
!   Test IO1-IO2 for a swap, and update IWK if necessary.
!
        IF (.NOT. SWPTST(N1, N2, IO1, IO2, X, Y, Z)) GO TO 4
        CALL SWAP(N1, N2, IO1, IO2, LIST, LPTR, LEND, LP21)
        IF (LP21 .EQ. 0) GO TO 9
        SWP = .TRUE.
        IWK(1, I) = N1
        IWK(2, I) = N2
4       CONTINUE
        IF (SWP) GO TO 1
!
! Successful termination.
!
5       NIT = ITER
        IER = 0
        RETURN
!
! MAXIT iterations performed without convergence.
!
6       NIT = MAXIT
        IER = 1
        RETURN
!
! Invalid input parameter.
!
7       NIT = 0
        IER = 2
        RETURN
!
! IO2 is not a neighbor of IO1.
!
8       NIT = ITER
        IER = 3
        RETURN
!
! Zero pointer returned by SWAP.
!
9       NIT = ITER
        IER = 4
        RETURN
    END
