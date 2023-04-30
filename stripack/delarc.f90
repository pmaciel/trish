SUBROUTINE DELARC(N, IO1, IO2, LIST, LPTR, LEND, LNEW, IER)
    INTEGER N, IO1, IO2, LIST(*), LPTR(*), LEND(N), LNEW, IER
!
!***********************************************************
!   This subroutine deletes a boundary arc from a triangula-
! tion.  It may be used to remove a null triangle from the
! convex hull boundary.  Note, however, that if the union of
! triangles is rendered nonconvex, Subroutines DELNOD, EDGE,
! and TRFIND (and hence ADDNOD) may fail.  Also, Function
! NEARND should not be called following an arc deletion.
!
!   This routine is identical to the similarly named routine
! in TRIPACK.
!
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 4.
!
!       IO1,IO2 = Indexes (in the range 1 to N) of a pair of
!                 adjacent boundary nodes defining the arc
!                 to be removed.
!
! The above parameters are not altered by this routine.
!
!       LIST,LPTR,LEND,LNEW = Triangulation data structure
!                             created by Subroutine TRMESH.
!
! On output:
!
!       LIST,LPTR,LEND,LNEW = Data structure updated with
!                             the removal of arc IO1-IO2
!                             unless IER > 0.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, IO1, or IO2 is outside its valid
!                     range, or IO1 = IO2.
!             IER = 2 if IO1-IO2 is not a boundary arc.
!             IER = 3 if the node opposite IO1-IO2 is al-
!                     ready a boundary node, and thus IO1
!                     or IO2 has only two neighbors or a
!                     deletion would result in two triangu-
!                     lations sharing a single node.
!             IER = 4 if one of the nodes is a neighbor of
!                     the other, but not vice versa, imply-
!                     ing an invalid triangulation data
!                     structure.
!
!***********************************************************
!
    INTEGER LSTPTR
    INTEGER LP, LPH, LPL, N1, N2, N3
!
! Local parameters:
!
! LP =       LIST pointer
! LPH =      LIST pointer or flag returned by DELNB
! LPL =      Pointer to the last neighbor of N1, N2, or N3
! N1,N2,N3 = Nodal indexes of a triangle such that N1->N2
!              is the directed boundary edge associated
!              with IO1-IO2
!
    N1 = IO1
    N2 = IO2
!
! Test for errors, and set N1->N2 to the directed boundary
!   edge associated with IO1-IO2:  (N1,N2,N3) is a triangle
!   for some N3.
!
    IF (N .LT. 4 .OR. N1 .LT. 1 .OR. N1 .GT. N .OR. N2 .LT. 1 .OR. N2 .GT. N .OR. N1 .EQ. N2) THEN
        IER = 1
        RETURN
    END IF
!
    LPL = LEND(N2)
    IF (-LIST(LPL) .NE. N1) THEN
        N1 = N2
        N2 = IO1
        LPL = LEND(N2)
        IF (-LIST(LPL) .NE. N1) THEN
            IER = 2
            RETURN
        END IF
    END IF
!
! Set N3 to the node opposite N1->N2 (the second neighbor
!   of N1), and test for error 3 (N3 already a boundary
!   node).
!
    LPL = LEND(N1)
    LP = LPTR(LPL)
    LP = LPTR(LP)
    N3 = ABS(LIST(LP))
    LPL = LEND(N3)
    IF (LIST(LPL) .LE. 0) THEN
        IER = 3
        RETURN
    END IF
!
! Delete N2 as a neighbor of N1, making N3 the first
!   neighbor, and test for error 4 (N2 not a neighbor
!   of N1).  Note that previously computed pointers may
!   no longer be valid following the call to DELNB.
!
    CALL DELNB(N1, N2, N, LIST, LPTR, LEND, LNEW, LPH)
    IF (LPH .LT. 0) THEN
        IER = 4
        RETURN
    END IF
!
! Delete N1 as a neighbor of N2, making N3 the new last
!   neighbor.
!
    CALL DELNB(N2, N1, N, LIST, LPTR, LEND, LNEW, LPH)
!
! Make N3 a boundary node with first neighbor N2 and last
!   neighbor N1.
!
    LP = LSTPTR(LEND(N3), N1, LIST, LPTR)
    LEND(N3) = LP
    LIST(LP) = -N1
!
! No errors encountered.
!
    IER = 0
    RETURN
END
