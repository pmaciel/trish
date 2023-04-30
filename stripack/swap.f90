SUBROUTINE SWAP(IN1, IN2, IO1, IO2, LIST, LPTR, LEND, LP21)
    INTEGER IN1, IN2, IO1, IO2, LIST(*), LPTR(*), LEND(*), LP21
!
!***********************************************************
!   Given a triangulation of a set of points on the unit
! sphere, this subroutine replaces a diagonal arc in a
! strictly convex quadrilateral (defined by a pair of adja-
! cent triangles) with the other diagonal.  Equivalently, a
! pair of adjacent triangles is replaced by another pair
! having the same union.
!
!
! On input:
!
!       IN1,IN2,IO1,IO2 = Nodal indexes of the vertices of
!                         the quadrilateral.  IO1-IO2 is re-
!                         placed by IN1-IN2.  (IO1,IO2,IN1)
!                         and (IO2,IO1,IN2) must be trian-
!                         gles on input.
!
! The above parameters are not altered by this routine.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to Subroutine
!                        TRMESH.
!
! On output:
!
!       LIST,LPTR,LEND = Data structure updated with the
!                        swap -- triangles (IO1,IO2,IN1) and
!                        (IO2,IO1,IN2) are replaced by
!                        (IN1,IN2,IO2) and (IN2,IN1,IO1)
!                        unless LP21 = 0.
!
!       LP21 = Index of IN1 as a neighbor of IN2 after the
!              swap is performed unless IN1 and IN2 are
!              adjacent on input, in which case LP21 = 0.
!
!***********************************************************
!
    INTEGER LSTPTR
    INTEGER LP, LPH, LPSAV
!
! Local parameters:
!
! LP,LPH,LPSAV = LIST pointers
!
!
! Test for IN1 and IN2 adjacent.
!
    LP = LSTPTR(LEND(IN1), IN2, LIST, LPTR)
    IF (ABS(LIST(LP)) .EQ. IN2) THEN
        LP21 = 0
        RETURN
    END IF
!
! Delete IO2 as a neighbor of IO1.
!
    LP = LSTPTR(LEND(IO1), IN2, LIST, LPTR)
    LPH = LPTR(LP)
    LPTR(LP) = LPTR(LPH)
!
! If IO2 is the last neighbor of IO1, make IN2 the
!   last neighbor.
!
    IF (LEND(IO1) .EQ. LPH) LEND(IO1) = LP
!
! Insert IN2 as a neighbor of IN1 following IO1
!   using the hole created above.
!
    LP = LSTPTR(LEND(IN1), IO1, LIST, LPTR)
    LPSAV = LPTR(LP)
    LPTR(LP) = LPH
    LIST(LPH) = IN2
    LPTR(LPH) = LPSAV
!
! Delete IO1 as a neighbor of IO2.
!
    LP = LSTPTR(LEND(IO2), IN1, LIST, LPTR)
    LPH = LPTR(LP)
    LPTR(LP) = LPTR(LPH)
!
! If IO1 is the last neighbor of IO2, make IN1 the
!   last neighbor.
!
    IF (LEND(IO2) .EQ. LPH) LEND(IO2) = LP
!
! Insert IN1 as a neighbor of IN2 following IO2.
!
    LP = LSTPTR(LEND(IN2), IO2, LIST, LPTR)
    LPSAV = LPTR(LP)
    LPTR(LP) = LPH
    LIST(LPH) = IN1
    LPTR(LPH) = LPSAV
    LP21 = LPH
    RETURN
END
