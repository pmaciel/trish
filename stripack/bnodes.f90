SUBROUTINE BNODES(N, LIST, LPTR, LEND, NODES, NB, NA, NT)
    INTEGER N, LIST(*), LPTR(*), LEND(N), NODES(*), NB, NA, NT
!
!***********************************************************
!   Given a triangulation of N nodes on the unit sphere
! created by Subroutine TRMESH, this subroutine returns an
! array containing the indexes (if any) of the counterclock-
! wise-ordered sequence of boundary nodes -- the nodes on
! the boundary of the convex hull of the set of nodes.  (The
! boundary is empty if the nodes do not lie in a single
! hemisphere.)  The numbers of boundary nodes, arcs, and
! triangles are also returned.
!
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to Subroutine
!                        TRMESH.
!
! The above parameters are not altered by this routine.
!
!       NODES = Integer array of length at least NB
!               (NB .LE. N).
!
! On output:
!
!       NODES = Ordered sequence of boundary node indexes
!               in the range 1 to N (in the first NB loca-
!               tions).
!
!       NB = Number of boundary nodes.
!
!       NA,NT = Number of arcs and triangles, respectively,
!               in the triangulation.
!
! Modules required by BNODES:  None
!
!***********************************************************
!
    INTEGER K, LP, N0, NN, NST
!
! Local parameters:
!
! K =   NODES index
! LP =  LIST pointer
! N0 =  Boundary node to be added to NODES
! NN =  Local copy of N
! NST = First element of nodes (arbitrarily chosen to be
!         the one with smallest index)
!
    NN = N
!
! Search for a boundary node.
!
    DO 1 NST = 1, NN
        LP = LEND(NST)
        IF (LIST(LP) .LT. 0) GO TO 2
1       CONTINUE
!
! The triangulation contains no boundary nodes.
!
        NB = 0
        NA = 3 * (NN - 2)
        NT = 2 * (NN - 2)
        RETURN
!
! NST is the first boundary node encountered.  Initialize
!   for traversal of the boundary.
!
2       NODES(1) = NST
        K = 1
        N0 = NST
!
! Traverse the boundary in counterclockwise order.
!
3       LP = LEND(N0)
        LP = LPTR(LP)
        N0 = LIST(LP)
        IF (N0 .EQ. NST) GO TO 4
        K = K + 1
        NODES(K) = N0
        GO TO 3
!
! Store the counts.
!
4       NB = K
        NT = 2 * N - NB - 2
        NA = NT + N - 1
        RETURN
    END
