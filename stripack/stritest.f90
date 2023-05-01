! By default, a triangulation is created from a set of N nodes consisting of the north pole and N-1 points uniformly distributed
! around the 60-degree parallel (with constant longitudinal separation). However, by enabling the READ statements below (C# in the
! first two columns), testing may be performed on an arbitrary set of up to NMAX nodes (refer to the PARAMETER statement below). A
! data set consists of the following sequence of records:
!
! N = Number of nodes (3 to NMAX)
! (RLAT(I),RLON(I), I = 1,N) = Nodal coordinates in degrees latitude (-90 to 90) and degrees longitude (-180 to 180)

PROGRAM MAIN
    IMPLICIT NONE

    INTEGER NMAX
    PARAMETER(NMAX=100)

    ! Number of nodes in the triangulation.  N .GE. 3.
    INTEGER N
    PARAMETER(N=9)

    ! Number of triangles (currently stored)
    ! NT = 2N-NB-2 if NB .GE. 3 or
    ! NT = 2N-4 if NB = 0, where NB is the number of boundary nodes.
    ! NT = (N1,N2,N3), where N1 < N2 and N1 < N3.
    INTEGER NT

    ! LIST,LPTR,LEND = Linked list data structure defining the triangulation
    ! Set of nodal indexes which, along with LPTR,
    ! LEND, and LNEW, define the triangulation as a
    ! set of N adjacency lists -- counterclockwise-
    ! ordered sequences of neighboring nodes such
    ! that the first and last neighbors of a bound-
    ! ary node are boundary nodes (the first neigh-
    ! bor of an interior node is arbitrary).  In
    ! order to distinguish between interior and
    ! boundary nodes, the last neighbor of each
    ! boundary node is represented by the negative
    ! of its index.
    ! Array of length at least 6N-12
    INTEGER LIST(6 * NMAX)

    ! Set of pointers (LIST indexes) in one-to-one
    ! correspondence with the elements of LIST.
    ! LIST(LPTR(I)) indexes the node which follows
    ! LIST(I) in cyclical counterclockwise order
    ! (the first neighbor follows the last neigh-
    ! bor).
    ! Array of length at least 6N-12
    INTEGER LPTR(6 * NMAX)

    ! Set of pointers to adjacency lists.  LEND(K)
    ! points to the last neighbor of node K for
    ! K = 1,...,N.  Thus, LIST(LEND(K)) < 0 if and
    ! only if K is a boundary node.
    ! Array of length at least N.
    INTEGER LEND(NMAX)

    ! Pointer to the first empty location in LIST
    ! and LPTR (list length plus one).  LIST, LPTR,
    ! LEND, and LNEW are not altered if IER < 0,
    ! and are incomplete if IER > 0.
    INTEGER LNEW

    ! Storage for the triangulation and nodal coordinates
    REAL RLAT(NMAX), RLON(NMAX)

    ! Arrays of length N containing the Cartesian
    ! coordinates of distinct nodes.  (X(K),Y(K),
    ! Z(K)) is referred to as node K, and K is re-
    ! ferred to as a nodal index.  It is required
    ! that X(K)**2 + Y(K)**2 + Z(K)**2 = 1 for all
    ! K.  The first three nodes must not be col-
    ! linear (lie on a common great circle).
    REAL X(NMAX), Y(NMAX), Z(NMAX)

    ! Storage for the triangle list
    !
    ! Integer array of length at least NROW*NT, where NT is at most 2N-4.  (A sufficient length is 12N if NROW=6)
    !
    ! NROW by NT array whose J-th column contains the vertex nodal indexes (first three rows),
    ! neighboring triangle indexes (second three rows), and, if NROW = 9, arc indexes (last
    ! three rows) associated with triangle J for J = 1,...,NT.  The vertices are ordered
    ! counterclockwise with the first vertex taken to be the one with smallest index.  Thus,
    ! LTRI(2,J) and LTRI(3,J) are larger than LTRI(1,J) and index adjacent neighbors of
    ! node LTRI(1,J).  For I = 1,2,3, LTRI(I+3,J) and LTRI(I+6,J) index the triangle and arc,
    ! respectively, which are opposite (not shared by) node LTRI(I,J), with LTRI(I+3,J) = 0 if
    ! LTRI(I+6,J) indexes a boundary arc.  Vertex indexes range from 1 to N, triangle indexes
    ! from 0 to NT, and, if included, arc indexes from 1 to NA, where NA = 3N-NB-3 if NB .GE. 3
    ! or 3N-6 if NB = 0.  The triangles are ordered on first (smallest) vertex indexes.
    INTEGER LTRI(6, 2 * NMAX)

    ! Work space arrays of length at
    ! least N.  The space is used to
    ! efficiently determine the nearest
    ! triangulation node to each un-
    ! processed node for use by ADDNOD.
    INTEGER NEAR(N)
    INTEGER NEXT(N)
    REAL DIST(N)

    ! Local (1):
    !
    ! D =        (Negative cosine of) distance from node K to node I
    ! D1,D2,D3 = Distances from node K to nodes 1, 2, and 3, respectively
    ! I,J =      Nodal indexes
    ! I0 =       Index of the node preceding I in a sequence of unprocessed nodes:  I = NEXT(I0)
    ! K =        Index of node to be added and DO-loop index: K > 3
    ! LP =       LIST index (pointer) of a neighbor of K
    ! LPL =      Pointer to the last neighbor of K
    ! NEXTI =    NEXT(I)
    ! NN =       Local copy of N
    ! I,J =      LTRI row indexes (1 to 3) associated with triangles NT and KN, respectively
    ! I1,I2,I3 = Nodal indexes of triangle KN
    ! ISV =      Variable used to permute indexes I1,I2,I3
    ! KN =       Index of the triangle that shares arc I1-I2 with NT
    ! NT =       Triangle index and number of currently stored triangles
    ! LP =       LIST pointer
    ! LP2 =      Pointer to N2 as a neighbor of N1
    ! LPL =      Pointer to the last neighbor of I1
    ! LPLN1 =    Pointer to the last neighbor of N1
    ! N1,N2,N3 = Nodal indexes of triangle NT
    INTEGER I, I1, I2, I3, ISV, J, K, KN, LP, LP2, LPL, LPLN1, N1, N2, N3, IER, I0, NEXTI, NN
    REAL D, D1, D2, D3, DEL, PHI, THETA
    LOGICAL LEFT

    ! Generate latitudinal and longitudinal coordinates. DEL is the separation in degrees between the nodes on the 60-degree line of latitude.
    RLAT(1) = 90.
    RLON(1) = 0.
    DEL = 360./REAL(N - 1)
    DO K = 2, N
        RLAT(K) = 60.
        RLON(K) = REAL(K - 2) * DEL
    END DO

    IF (N .LT. 3 .OR. N .GT. NMAX) THEN
        PRINT *, 'N is outside its valid range:  N =' , N
        STOP
    END IF

    ! Set X,Y,Z on the unit sphere from RLON,RLAT [radian]
    DO K = 1, N
        PHI = ACOS(-1.) / 180. * RLAT(K)
        THETA = ACOS(-1.) / 180. * RLON(K)
        X(K) = COS(PHI) * COS(THETA)
        Y(K) = COS(PHI) * SIN(THETA)
        Z(K) = SIN(PHI)
    END DO



    ! Create the triangulation

    ! Store the first triangle in the linked list.
    IF (.NOT. LEFT(X(1), Y(1), Z(1), X(2), Y(2), Z(2), X(3), Y(3), Z(3))) THEN
    !
    !   The first triangle is (3,2,1) = (2,1,3) = (1,3,2).
    !
        LIST(1) = 3
        LPTR(1) = 2
        LIST(2) = -2
        LPTR(2) = 1
        LEND(1) = 2
    !
        LIST(3) = 1
        LPTR(3) = 4
        LIST(4) = -3
        LPTR(4) = 3
        LEND(2) = 4
    !
        LIST(5) = 2
        LPTR(5) = 6
        LIST(6) = -1
        LPTR(6) = 5
        LEND(3) = 6
    !
    ELSEIF (.NOT. LEFT(X(2), Y(2), Z(2), X(1), Y(1), Z(1), X(3), Y(3), Z(3))) THEN
    !
    !   The first triangle is (1,2,3):  3 Strictly Left 1->2,
    !     i.e., node 3 lies in the left hemisphere defined by
    !     arc 1->2.
    !
        LIST(1) = 2
        LPTR(1) = 2
        LIST(2) = -3
        LPTR(2) = 1
        LEND(1) = 2
    !
        LIST(3) = 3
        LPTR(3) = 4
        LIST(4) = -1
        LPTR(4) = 3
        LEND(2) = 4
    !
        LIST(5) = 1
        LPTR(5) = 6
        LIST(6) = -2
        LPTR(6) = 5
        LEND(3) = 6
    !
    ELSE
    !
    !   The first three nodes are collinear.
        PRINT *, "the first three nodes are collinear"
        STOP
    END IF
    !
    ! Initialize LNEW and test for N = 3.
    !
    LNEW = 7
    IF (NN .GT. 3) THEN
        ! A nearest-node data structure (NEAR, NEXT, and DIST) is
        !   used to obtain an expected-time (N*log(N)) incremental
        !   algorithm by enabling constant search time for locating
        !   each new node in the triangulation.
        !
        ! For each unprocessed node K, NEAR(K) is the index of the
        !   triangulation node closest to K (used as the starting
        !   point for the search in Subroutine TRFIND) and DIST(K)
        !   is an increasing function of the arc length (angular
        !   distance) between nodes K and NEAR(K):  -Cos(a) for arc
        !   length a.
        !
        ! Since it is necessary to efficiently find the subset of
        !   unprocessed nodes associated with each triangulation
        !   node J (those that have J as their NEAR entries), the
        !   subsets are stored in NEAR and NEXT as follows:  for
        !   each node J in the triangulation, I = NEAR(J) is the
        !   first unprocessed node in J's set (with I = 0 if the
        !   set is empty), L = NEXT(I) (if I > 0) is the second,
        !   NEXT(L) (if L > 0) is the third, etc.  The nodes in each
        !   set are initially ordered by increasing indexes (which
        !   maximizes efficiency) but that ordering is not main-
        !   tained as the data structure is updated.
        !
        ! Initialize the data structure for the single triangle.
        !
        NEAR(1) = 0
        NEAR(2) = 0
        NEAR(3) = 0
        DO 1 K = NN, 4, -1
            D1 = -(X(K) * X(1) + Y(K) * Y(1) + Z(K) * Z(1))
            D2 = -(X(K) * X(2) + Y(K) * Y(2) + Z(K) * Z(2))
            D3 = -(X(K) * X(3) + Y(K) * Y(3) + Z(K) * Z(3))
            IF (D1 .LE. D2 .AND. D1 .LE. D3) THEN
                NEAR(K) = 1
                DIST(K) = D1
                NEXT(K) = NEAR(1)
                NEAR(1) = K
            ELSEIF (D2 .LE. D1 .AND. D2 .LE. D3) THEN
                NEAR(K) = 2
                DIST(K) = D2
                NEXT(K) = NEAR(2)
                NEAR(2) = K
            ELSE
                NEAR(K) = 3
                DIST(K) = D3
                NEXT(K) = NEAR(3)
                NEAR(3) = K
            END IF
        1       CONTINUE
        !
        ! Add the remaining nodes
        !
            DO 6 K = 4, NN
                CALL ADDNOD(NEAR(K), K, X, Y, Z, LIST, LPTR, LEND, LNEW, IER)
                IF (IER .NE. 0) THEN
                    PRINT *, "Error: ADDNOD"
                    STOP
                END IF

        ! Remove K from the set of unprocessed nodes associated
        !   with NEAR(K).
        !
                I = NEAR(K)
                IF (NEAR(I) .EQ. K) THEN
                    NEAR(I) = NEXT(K)
                ELSE
                    I = NEAR(I)
        2               I0 = I
                    I = NEXT(I0)
                    IF (I .NE. K) GO TO 2
                    NEXT(I0) = NEXT(K)
                END IF
                NEAR(K) = 0
        !
        ! Loop on neighbors J of node K.
        !
                LPL = LEND(K)
                LP = LPL
        3           LP = LPTR(LP)
                J = ABS(LIST(LP))
        !
        ! Loop on elements I in the sequence of unprocessed nodes
        !   associated with J:  K is a candidate for replacing J
        !   as the nearest triangulation node to I.  The next value
        !   of I in the sequence, NEXT(I), must be saved before I
        !   is moved because it is altered by adding I to K's set.

        I = NEAR(J)
        4   IF (I .EQ. 0) GO TO 5
        NEXTI = NEXT(I)

        ! Test for the distance from I to K less than the distance from I to J.
        D = -(X(I) * X(K) + Y(I) * Y(K) + Z(I) * Z(K))
        IF (D .LT. DIST(I)) THEN
            ! Replace J by K as the nearest triangulation node to I:
            !   update NEAR(I) and DIST(I), and remove I from J's set
            !   of unprocessed nodes and add it to K's set.
            NEAR(I) = K
            DIST(I) = D
            IF (I .EQ. NEAR(J)) THEN
                NEAR(J) = NEXTI
            ELSE
                NEXT(I0) = NEXTI
            END IF
            NEXT(I) = NEAR(K)
            NEAR(K) = I
        ELSE
            I0 = I
        END IF

        ! Bottom of loop on I.
        I = NEXTI
        GO TO 4

        ! Bottom of loop on neighbors J.
        5   IF (LP .NE. LPL) GO TO 3
        6   CONTINUE
    END IF





    ! Converts a triangulation data structure from linked list to a triangle list
    ! Loop on nodes N1.
    NT = 0
    DO N1 = 1, N-2
        ! Loop on pairs of adjacent neighbors (N2,N3).  LPLN1 points
        !   to the last neighbor of N1, and LP2 points to N2.
        LPLN1 = LEND(N1)
        LP2 = LPLN1
21      LP2 = LPTR(LP2)
        N2 = LIST(LP2)
        LP = LPTR(LP2)
        N3 = ABS(LIST(LP))
        IF (N2 .LT. N1 .OR. N3 .LT. N1) GO TO 28

        ! Add a new triangle NT = (N1,N2,N3).
        NT = NT + 1
        LTRI(1, NT) = N1
        LTRI(2, NT) = N2
        LTRI(3, NT) = N3

        ! Loop on triangle sides (I2,I1) with neighboring triangles KN = (I1,I2,I3).
        DO 27 I = 1, 3
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

        ! Set I3 to the neighbor of I1 that follows I2 unless
        !   I2->I1 is a boundary arc.
        LPL = LEND(I1)
        LP = LPTR(LPL)
22       IF (LIST(LP) .EQ. I2) GO TO 23
        LP = LPTR(LP)
        IF (LP .NE. LPL) GO TO 22

        ! I2 is the last neighbor of I1 unless the data structure is invalid.  Bypass the search for a neighboring
        !   triangle if I2->I1 is a boundary arc.
        IF (ABS(LIST(LP)) .NE. I2) THEN
            PRINT *, "Invalid:  I1 is a neighbor of I2, but I2 is not a neighbor of I1."
            PRINT *, "triangulation data structure (LIST,LPTR,LEND) is invalid"
            STOP
        END IF

        KN = 0
        IF (LIST(LP) .LT. 0) GO TO 26

        ! I2->I1 is not a boundary arc, and LP points to I2 as a neighbor of I1.
23      LP = LPTR(LP)
        I3 = ABS(LIST(LP))

        ! Find J such that LTRI(J,KN) = I3 (not used if KN > NT), and permute the vertex indexes of KN so that I1 is
        ! smallest.
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

        ! Test for KN > NT (triangle index not yet assigned).
        IF (I1 .GT. N1) GO TO 27

        ! Find KN, if it exists, by searching the triangle list in reverse order.
        DO KN = NT - 1, 1, -1
            IF (LTRI(1, KN) .EQ. I1 .AND. LTRI(2, KN) .EQ. I2 .AND. LTRI(3, KN) .EQ. I3) GO TO 25
        END DO
        GO TO 27

        ! Store NT as a neighbor of KN.
25      LTRI(J + 3, KN) = NT

        ! Store KN as a neighbor of NT
26      LTRI(I + 3, NT) = KN
27      CONTINUE

        ! Bottom of loop on triangles.
28      IF (LP2 .NE. LPLN1) GO TO 21
    END DO
END
