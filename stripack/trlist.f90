SUBROUTINE TRLIST(N, LIST, LPTR, LEND, NROW, NT, LTRI, IER)
    INTEGER N, LIST(*), LPTR(*), LEND(N), NROW, NT, LTRI(NROW, *), IER
!
!***********************************************************
!   This subroutine converts a triangulation data structure
! from the linked list created by Subroutine TRMESH to a
! triangle list.
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       LIST,LPTR,LEND = Linked list data structure defin-
!                        ing the triangulation.  Refer to
!                        Subroutine TRMESH.
!
!       NROW = Number of rows (entries per triangle) re-
!              served for the triangle list LTRI.  The value
!              must be 6 if only the vertex indexes and
!              neighboring triangle indexes are to be
!              stored, or 9 if arc indexes are also to be
!              assigned and stored.  Refer to LTRI.
!
! The above parameters are not altered by this routine.
!
!       LTRI = Integer array of length at least NROW*NT,
!              where NT is at most 2N-4.  (A sufficient
!              length is 12N if NROW=6 or 18N if NROW=9.)
!
! On output:
!
!       NT = Number of triangles in the triangulation unless
!            IER .NE. 0, in which case NT = 0.  NT = 2N-NB-2
!            if NB .GE. 3 or 2N-4 if NB = 0, where NB is the
!            number of boundary nodes.
!
!       LTRI = NROW by NT array whose J-th column contains
!              the vertex nodal indexes (first three rows),
!              neighboring triangle indexes (second three
!              rows), and, if NROW = 9, arc indexes (last
!              three rows) associated with triangle J for
!              J = 1,...,NT.  The vertices are ordered
!              counterclockwise with the first vertex taken
!              to be the one with smallest index.  Thus,
!              LTRI(2,J) and LTRI(3,J) are larger than
!              LTRI(1,J) and index adjacent neighbors of
!              node LTRI(1,J).  For I = 1,2,3, LTRI(I+3,J)
!              and LTRI(I+6,J) index the triangle and arc,
!              respectively, which are opposite (not shared
!              by) node LTRI(I,J), with LTRI(I+3,J) = 0 if
!              LTRI(I+6,J) indexes a boundary arc.  Vertex
!              indexes range from 1 to N, triangle indexes
!              from 0 to NT, and, if included, arc indexes
!              from 1 to NA, where NA = 3N-NB-3 if NB .GE. 3
!              or 3N-6 if NB = 0.  The triangles are or-
!              dered on first (smallest) vertex indexes.
!
!       IER = Error indicator.
!             IER = 0 if no errors were encountered.
!             IER = 1 if N or NROW is outside its valid
!                     range on input.
!             IER = 2 if the triangulation data structure
!                     (LIST,LPTR,LEND) is invalid.  Note,
!                     however, that these arrays are not
!                     completely tested for validity.
!
!***********************************************************
!
    INTEGER I, I1, I2, I3, ISV, J, KA, KN, KT, LP, LP2, LPL, LPLN1, N1, N2, N3, NM2
    LOGICAL ARCS
!
! Local parameters:
!
! ARCS =     Logical variable with value TRUE iff are
!              indexes are to be stored
! I,J =      LTRI row indexes (1 to 3) associated with
!              triangles KT and KN, respectively
! I1,I2,I3 = Nodal indexes of triangle KN
! ISV =      Variable used to permute indexes I1,I2,I3
! KA =       Arc index and number of currently stored arcs
! KN =       Index of the triangle that shares arc I1-I2
!              with KT
! KT =       Triangle index and number of currently stored
!              triangles
! LP =       LIST pointer
! LP2 =      Pointer to N2 as a neighbor of N1
! LPL =      Pointer to the last neighbor of I1
! LPLN1 =    Pointer to the last neighbor of N1
! N1,N2,N3 = Nodal indexes of triangle KT
! NM2 =      N-2
!
!
! Test for invalid input parameters.
!
    IF (N .LT. 3 .OR. (NROW .NE. 6 .AND. NROW .NE. 9)) GO TO 11
!
! Initialize parameters for loop on triangles KT = (N1,N2,
!   N3), where N1 < N2 and N1 < N3.
!
!   ARCS = TRUE iff arc indexes are to be stored.
!   KA,KT = Numbers of currently stored arcs and triangles.
!   NM2 = Upper bound on candidates for N1.
!
    ARCS = NROW .EQ. 9
    KA = 0
    KT = 0
    NM2 = N - 2
!
! Loop on nodes N1.
!
    DO 9 N1 = 1, NM2
!
! Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
!   to the last neighbor of N1, and LP2 points to N2.
!
        LPLN1 = LEND(N1)
        LP2 = LPLN1
1       LP2 = LPTR(LP2)
        N2 = LIST(LP2)
        LP = LPTR(LP2)
        N3 = ABS(LIST(LP))
        IF (N2 .LT. N1 .OR. N3 .LT. N1) GO TO 8
!
! Add a new triangle KT = (N1,N2,N3).
!
        KT = KT + 1
        LTRI(1, KT) = N1
        LTRI(2, KT) = N2
        LTRI(3, KT) = N3
!
! Loop on triangle sides (I2,I1) with neighboring triangles
!   KN = (I1,I2,I3).
!
        DO 7 I = 1, 3
        IF (I .EQ. 1) THEN
            I1 = N3
            I2 = N2
        ELSEIF (I .EQ. 2) THEN
            I1 = N1
            I2 = N3
        ELSE
            I1 = N2
            I2 = N1
        END IF
!
! Set I3 to the neighbor of I1 that follows I2 unless
!   I2->I1 is a boundary arc.
!
        LPL = LEND(I1)
        LP = LPTR(LPL)
2       IF (LIST(LP) .EQ. I2) GO TO 3
        LP = LPTR(LP)
        IF (LP .NE. LPL) GO TO 2
!
!   I2 is the last neighbor of I1 unless the data structure
!     is invalid.  Bypass the search for a neighboring
!     triangle if I2->I1 is a boundary arc.
!
        IF (ABS(LIST(LP)) .NE. I2) GO TO 12
        KN = 0
        IF (LIST(LP) .LT. 0) GO TO 6
!
!   I2->I1 is not a boundary arc, and LP points to I2 as
!     a neighbor of I1.
!
3       LP = LPTR(LP)
        I3 = ABS(LIST(LP))
!
! Find J such that LTRI(J,KN) = I3 (not used if KN > KT),
!   and permute the vertex indexes of KN so that I1 is
!   smallest.
!
        IF (I1 .LT. I2 .AND. I1 .LT. I3) THEN
            J = 3
        ELSEIF (I2 .LT. I3) THEN
            J = 2
            ISV = I1
            I1 = I2
            I2 = I3
            I3 = ISV
        ELSE
            J = 1
            ISV = I1
            I1 = I3
            I3 = I2
            I2 = ISV
        END IF
!
! Test for KN > KT (triangle index not yet assigned).
!
        IF (I1 .GT. N1) GO TO 7
!
! Find KN, if it exists, by searching the triangle list in
!   reverse order.
!
        DO 4 KN = KT - 1, 1, -1
            IF (LTRI(1, KN) .EQ. I1 .AND. LTRI(2, KN) .EQ. I2 .AND. LTRI(3, KN) .EQ. I3) GO TO 5
4           CONTINUE
            GO TO 7
!
! Store KT as a neighbor of KN.
!
5           LTRI(J + 3, KN) = KT
!
! Store KN as a neighbor of KT, and add a new arc KA.
!
6           LTRI(I + 3, KT) = KN
            IF (ARCS) THEN
                KA = KA + 1
                LTRI(I + 6, KT) = KA
                IF (KN .NE. 0) LTRI(J + 6, KN) = KA
            END IF
7           CONTINUE
!
! Bottom of loop on triangles.
!
8           IF (LP2 .NE. LPLN1) GO TO 1
9           CONTINUE
!
! No errors encountered.
!
            NT = KT
            IER = 0
            RETURN
!
! Invalid input parameter.
!
11          NT = 0
            IER = 1
            RETURN
!
! Invalid triangulation data structure:  I1 is a neighbor of
!   I2, but I2 is not a neighbor of I1.
!
12          NT = 0
            IER = 2
            RETURN
        END
