SUBROUTINE TRMESH(N, X, Y, Z, LIST, LPTR, LEND, LNEW, NEAR, NEXT, DIST, IER)
    INTEGER N, LIST(*), LPTR(*), LEND(N), LNEW, NEAR(N), NEXT(N), IER
    REAL X(N), Y(N), Z(N), DIST(N)
!
!***********************************************************
!   This subroutine creates a Delaunay triangulation of a
! set of N arbitrarily distributed points, referred to as
! nodes, on the surface of the unit sphere.  The Delaunay
! triangulation is defined as a set of (spherical) triangles
! with the following five properties:
!
!  1)  The triangle vertices are nodes.
!  2)  No triangle contains a node other than its vertices.
!  3)  The interiors of the triangles are pairwise disjoint.
!  4)  The union of triangles is the convex hull of the set
!        of nodes (the smallest convex set that contains
!        the nodes).  If the nodes are not contained in a
!        single hemisphere, their convex hull is the en-
!        tire sphere and there are no boundary nodes.
!        Otherwise, there are at least three boundary nodes.
!  5)  The interior of the circumcircle of each triangle
!        contains no node.
!
! The first four properties define a triangulation, and the
! last property results in a triangulation which is as close
! as possible to equiangular in a certain sense and which is
! uniquely defined unless four or more nodes lie in a common
! plane.  This property makes the triangulation well-suited
! for solving closest-point problems and for triangle-based
! interpolation.
!
!   Provided the nodes are randomly ordered, the algorithm
! has expected time complexity O(N*log(N)) for most nodal
! distributions.  Note, however, that the complexity may be
! as high as O(N**2) if, for example, the nodes are ordered
! on increasing latitude.
!
!   Spherical coordinates (latitude and longitude) may be
! converted to Cartesian coordinates by Subroutine TRANS.
!
!   The following is a list of the software package modules
! which a user may wish to call directly:
!
!  ADDNOD - Updates the triangulation by appending a new
!             node.
!
!  AREAS  - Returns the area of a spherical triangle.
!
!  BNODES - Returns an array containing the indexes of the
!             boundary nodes (if any) in counterclockwise
!             order.  Counts of boundary nodes, triangles,
!             and arcs are also returned.
!
!  CIRCUM - Returns the circumcenter of a spherical trian-
!             gle.
!
!  CRLIST - Returns the set of triangle circumcenters
!             (Voronoi vertices) and circumradii associated
!             with a triangulation.
!
!  DELARC - Deletes a boundary arc from a triangulation.
!
!  DELNOD - Updates the triangulation with a nodal deletion.
!
!  EDGE   - Forces an arbitrary pair of nodes to be connec-
!             ted by an arc in the triangulation.
!
!  GETNP  - Determines the ordered sequence of L closest
!             nodes to a given node, along with the associ-
!             ated distances.
!
!  INSIDE - Locates a point relative to a polygon on the
!             surface of the sphere.
!
!  INTRSC - Returns the point of intersection between a
!             pair of great circle arcs.
!
!  JRAND  - Generates a uniformly distributed pseudo-random
!             integer.
!
!  LEFT   - Locates a point relative to a great circle.
!
!  NEARND - Returns the index of the nearest node to an
!             arbitrary point, along with its squared
!             distance.
!
!  SCOORD - Converts a point from Cartesian coordinates to
!             spherical coordinates.
!
!  STORE  - Forces a value to be stored in main memory so
!             that the precision of floating point numbers
!             in memory locations rather than registers is
!             computed.
!
!  TRANS  - Transforms spherical coordinates into Cartesian
!             coordinates on the unit sphere for input to
!             Subroutine TRMESH.
!
!  TRLIST - Converts the triangulation data structure to a
!             triangle list more suitable for use in a fin-
!             ite element code.
!
!  TRLPRT - Prints the triangle list created by Subroutine
!             TRLIST.
!
!  TRMESH - Creates a Delaunay triangulation of a set of
!             nodes.
!
!  TRPLOT - Creates a level-2 Encapsulated Postscript (EPS)
!             file containing a triangulation plot.
!
!  TRPRNT - Prints the triangulation data structure and,
!             optionally, the nodal coordinates.
!
!  VRPLOT - Creates a level-2 Encapsulated Postscript (EPS)
!             file containing a Voronoi diagram plot.
!
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       X,Y,Z = Arrays of length N containing the Cartesian
!               coordinates of distinct nodes.  (X(K),Y(K),
!               Z(K)) is referred to as node K, and K is re-
!               ferred to as a nodal index.  It is required
!               that X(K)**2 + Y(K)**2 + Z(K)**2 = 1 for all
!               K.  The first three nodes must not be col-
!               linear (lie on a common great circle).
!
! The above parameters are not altered by this routine.
!
!       LIST,LPTR = Arrays of length at least 6N-12.
!
!       LEND = Array of length at least N.
!
!       NEAR,NEXT,DIST = Work space arrays of length at
!                        least N.  The space is used to
!                        efficiently determine the nearest
!                        triangulation node to each un-
!                        processed node for use by ADDNOD.
!
! On output:
!
!       LIST = Set of nodal indexes which, along with LPTR,
!              LEND, and LNEW, define the triangulation as a
!              set of N adjacency lists -- counterclockwise-
!              ordered sequences of neighboring nodes such
!              that the first and last neighbors of a bound-
!              ary node are boundary nodes (the first neigh-
!              bor of an interior node is arbitrary).  In
!              order to distinguish between interior and
!              boundary nodes, the last neighbor of each
!              boundary node is represented by the negative
!              of its index.
!
!       LPTR = Set of pointers (LIST indexes) in one-to-one
!              correspondence with the elements of LIST.
!              LIST(LPTR(I)) indexes the node which follows
!              LIST(I) in cyclical counterclockwise order
!              (the first neighbor follows the last neigh-
!              bor).
!
!       LEND = Set of pointers to adjacency lists.  LEND(K)
!              points to the last neighbor of node K for
!              K = 1,...,N.  Thus, LIST(LEND(K)) < 0 if and
!              only if K is a boundary node.
!
!       LNEW = Pointer to the first empty location in LIST
!              and LPTR (list length plus one).  LIST, LPTR,
!              LEND, and LNEW are not altered if IER < 0,
!              and are incomplete if IER > 0.
!
!       NEAR,NEXT,DIST = Garbage.
!
!       IER = Error indicator:
!             IER =  0 if no errors were encountered.
!             IER = -1 if N < 3 on input.
!             IER = -2 if the first three nodes are
!                      collinear.
!             IER =  L if nodes L and M coincide for some
!                      M > L.  The data structure represents
!                      a triangulation of nodes 1 to M-1 in
!                      this case.
!
!***********************************************************
!
    INTEGER I, I0, J, K, LP, LPL, NEXTI, NN
    LOGICAL LEFT
    REAL D, D1, D2, D3
!
! Local parameters:
!
! D =        (Negative cosine of) distance from node K to
!              node I
! D1,D2,D3 = Distances from node K to nodes 1, 2, and 3,
!              respectively
! I,J =      Nodal indexes
! I0 =       Index of the node preceding I in a sequence of
!              unprocessed nodes:  I = NEXT(I0)
! K =        Index of node to be added and DO-loop index:
!              K > 3
! LP =       LIST index (pointer) of a neighbor of K
! LPL =      Pointer to the last neighbor of K
! NEXTI =    NEXT(I)
! NN =       Local copy of N
!
    NN = N
    IF (NN .LT. 3) THEN
        IER = -1
        RETURN
    END IF
!
! Store the first triangle in the linked list.
!
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
!
        IER = -2
        RETURN
    END IF
!
! Initialize LNEW and test for N = 3.
!
    LNEW = 7
    IF (NN .EQ. 3) THEN
        IER = 0
        RETURN
    END IF
!
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
            IF (IER .NE. 0) RETURN
!
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
!
            I = NEAR(J)
4           IF (I .EQ. 0) GO TO 5
            NEXTI = NEXT(I)
!
! Test for the distance from I to K less than the distance
!   from I to J.
!
            D = -(X(I) * X(K) + Y(I) * Y(K) + Z(I) * Z(K))
            IF (D .LT. DIST(I)) THEN
!
! Replace J by K as the nearest triangulation node to I:
!   update NEAR(I) and DIST(I), and remove I from J's set
!   of unprocessed nodes and add it to K's set.
!
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
!
! Bottom of loop on I.
!
            I = NEXTI
            GO TO 4
!
! Bottom of loop on neighbors J.
!
5           IF (LP .NE. LPL) GO TO 3
6           CONTINUE
            RETURN
        END
