SUBROUTINE CRLIST(N, NCOL, X, Y, Z, LIST, LEND, LPTR, LNEW, LTRI, LISTC, NB, XC, YC, ZC, RC, IER)
    INTEGER N, NCOL, LIST(*), LEND(N), LPTR(*), LNEW, LTRI(6, NCOL), LISTC(*), NB, IER
    REAL X(N), Y(N), Z(N), XC(*), YC(*), ZC(*), RC(*)
!
!***********************************************************
!   Given a Delaunay triangulation of nodes on the surface
! of the unit sphere, this subroutine returns the set of
! triangle circumcenters corresponding to Voronoi vertices,
! along with the circumradii and a list of triangle indexes
! LISTC stored in one-to-one correspondence with LIST/LPTR
! entries.
!
!   A triangle circumcenter is the point (unit vector) lying
! at the same angular distance from the three vertices and
! contained in the same hemisphere as the vertices.  (Note
! that the negative of a circumcenter is also equidistant
! from the vertices.)  If the triangulation covers the sur-
! face, the Voronoi vertices are the circumcenters of the
! triangles in the Delaunay triangulation.  LPTR, LEND, and
! LNEW are not altered in this case.
!
!   On the other hand, if the nodes are contained in a sin-
! gle hemisphere, the triangulation is implicitly extended
! to the entire surface by adding pseudo-arcs (of length
! greater than 180 degrees) between boundary nodes forming
! pseudo-triangles whose 'circumcenters' are included in the
! list.  This extension to the triangulation actually con-
! sists of a triangulation of the set of boundary nodes in
! which the swap test is reversed (a non-empty circumcircle
! test).  The negative circumcenters are stored as the
! pseudo-triangle 'circumcenters'.  LISTC, LPTR, LEND, and
! LNEW contain a data structure corresponding to the ex-
! tended triangulation (Voronoi diagram), but LIST is not
! altered in this case.  Thus, if it is necessary to retain
! the original (unextended) triangulation data structure,
! copies of LPTR and LNEW must be saved before calling this
! routine.
!
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 3.
!           Note that, if N = 3, there are only two Voronoi
!           vertices separated by 180 degrees, and the
!           Voronoi regions are not well defined.
!
!       NCOL = Number of columns reserved for LTRI.  This
!              must be at least NB-2, where NB is the number
!              of boundary nodes.
!
!       X,Y,Z = Arrays of length N containing the Cartesian
!               coordinates of the nodes (unit vectors).
!
!       LIST = Integer array containing the set of adjacency
!              lists.  Refer to Subroutine TRMESH.
!
!       LEND = Set of pointers to ends of adjacency lists.
!              Refer to Subroutine TRMESH.
!
! The above parameters are not altered by this routine.
!
!       LPTR = Array of pointers associated with LIST.  Re-
!              fer to Subroutine TRMESH.
!
!       LNEW = Pointer to the first empty location in LIST
!              and LPTR (list length plus one).
!
!       LTRI = Integer work space array dimensioned 6 by
!              NCOL, or unused dummy parameter if NB = 0.
!
!       LISTC = Integer array of length at least 3*NT, where
!               NT = 2*N-4 is the number of triangles in the
!               triangulation (after extending it to cover
!               the entire surface if necessary).
!
!       XC,YC,ZC,RC = Arrays of length NT = 2*N-4.
!
! On output:
!
!       LPTR = Array of pointers associated with LISTC:
!              updated for the addition of pseudo-triangles
!              if the original triangulation contains
!              boundary nodes (NB > 0).
!
!       LNEW = Pointer to the first empty location in LISTC
!              and LPTR (list length plus one).  LNEW is not
!              altered if NB = 0.
!
!       LTRI = Triangle list whose first NB-2 columns con-
!              tain the indexes of a clockwise-ordered
!              sequence of vertices (first three rows)
!              followed by the LTRI column indexes of the
!              triangles opposite the vertices (or 0
!              denoting the exterior region) in the last
!              three rows.  This array is not generally of
!              any use.
!
!       LISTC = Array containing triangle indexes (indexes
!               to XC, YC, ZC, and RC) stored in 1-1 corres-
!               pondence with LIST/LPTR entries (or entries
!               that would be stored in LIST for the
!               extended triangulation):  the index of tri-
!               angle (N1,N2,N3) is stored in LISTC(K),
!               LISTC(L), and LISTC(M), where LIST(K),
!               LIST(L), and LIST(M) are the indexes of N2
!               as a neighbor of N1, N3 as a neighbor of N2,
!               and N1 as a neighbor of N3.  The Voronoi
!               region associated with a node is defined by
!               the CCW-ordered sequence of circumcenters in
!               one-to-one correspondence with its adjacency
!               list (in the extended triangulation).
!
!       NB = Number of boundary nodes unless IER = 1.
!
!       XC,YC,ZC = Arrays containing the Cartesian coordi-
!                  nates of the triangle circumcenters
!                  (Voronoi vertices).  XC(I)**2 + YC(I)**2
!                  + ZC(I)**2 = 1.  The first NB-2 entries
!                  correspond to pseudo-triangles if NB > 0.
!
!       RC = Array containing circumradii (the arc lengths
!            or angles between the circumcenters and associ-
!            ated triangle vertices) in 1-1 correspondence
!            with circumcenters.
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N < 3.
!             IER = 2 if NCOL < NB-2.
!             IER = 3 if a triangle is degenerate (has ver-
!                     tices lying on a common geodesic).
!
!***********************************************************
!
    INTEGER LSTPTR
    INTEGER I1, I2, I3, I4, IERR, KT, KT1, KT2, KT11, KT12, KT21, KT22, LP, LPL, LPN, N0, N1, N2, N3, N4, NM2, NN, NT
    LOGICAL SWPTST
    LOGICAL SWP
    REAL C(3), T, V1(3), V2(3), V3(3)
!
! Local parameters:
!
! C =         Circumcenter returned by Subroutine CIRCUM
! I1,I2,I3 =  Permutation of (1,2,3):  LTRI row indexes
! I4 =        LTRI row index in the range 1 to 3
! IERR =      Error flag for calls to CIRCUM
! KT =        Triangle index
! KT1,KT2 =   Indexes of a pair of adjacent pseudo-triangles
! KT11,KT12 = Indexes of the pseudo-triangles opposite N1
!               and N2 as vertices of KT1
! KT21,KT22 = Indexes of the pseudo-triangles opposite N1
!               and N2 as vertices of KT2
! LP,LPN =    LIST pointers
! LPL =       LIST pointer of the last neighbor of N1
! N0 =        Index of the first boundary node (initial
!               value of N1) in the loop on boundary nodes
!               used to store the pseudo-triangle indexes
!               in LISTC
! N1,N2,N3 =  Nodal indexes defining a triangle (CCW order)
!               or pseudo-triangle (clockwise order)
! N4 =        Index of the node opposite N2 -> N1
! NM2 =       N-2
! NN =        Local copy of N
! NT =        Number of pseudo-triangles:  NB-2
! SWP =       Logical variable set to TRUE in each optimiza-
!               tion loop (loop on pseudo-arcs) iff a swap
!               is performed
! V1,V2,V3 =  Vertices of triangle KT = (N1,N2,N3) sent to
!               Subroutine CIRCUM
!
    NN = N
    NB = 0
    NT = 0
    IF (NN .LT. 3) GO TO 21
!
! Search for a boundary node N1.
!
    DO 1 N1 = 1, NN
        LP = LEND(N1)
        IF (LIST(LP) .LT. 0) GO TO 2
1       CONTINUE
!
! The triangulation already covers the sphere.
!
        GO TO 9
!
! There are NB .GE. 3 boundary nodes.  Add NB-2 pseudo-
!   triangles (N1,N2,N3) by connecting N3 to the NB-3
!   boundary nodes to which it is not already adjacent.
!
!   Set N3 and N2 to the first and last neighbors,
!     respectively, of N1.
!
2       N2 = -LIST(LP)
        LP = LPTR(LP)
        N3 = LIST(LP)
!
!   Loop on boundary arcs N1 -> N2 in clockwise order,
!     storing triangles (N1,N2,N3) in column NT of LTRI
!     along with the indexes of the triangles opposite
!     the vertices.
!
3       NT = NT + 1
        IF (NT .LE. NCOL) THEN
            LTRI(1, NT) = N1
            LTRI(2, NT) = N2
            LTRI(3, NT) = N3
            LTRI(4, NT) = NT + 1
            LTRI(5, NT) = NT - 1
            LTRI(6, NT) = 0
        END IF
        N1 = N2
        LP = LEND(N1)
        N2 = -LIST(LP)
        IF (N2 .NE. N3) GO TO 3
!
        NB = NT + 2
        IF (NCOL .LT. NT) GO TO 22
        LTRI(4, NT) = 0
        IF (NT .EQ. 1) GO TO 7
!
! Optimize the exterior triangulation (set of pseudo-
!   triangles) by applying swaps to the pseudo-arcs N1-N2
!   (pairs of adjacent pseudo-triangles KT1 and KT2 > KT1).
!   The loop on pseudo-arcs is repeated until no swaps are
!   performed.
!
4       SWP = .FALSE.
        DO 6 KT1 = 1, NT - 1
        DO 5 I3 = 1, 3
            KT2 = LTRI(I3 + 3, KT1)
            IF (KT2 .LE. KT1) GO TO 5
!
!   The LTRI row indexes (I1,I2,I3) of triangle KT1 =
!     (N1,N2,N3) are a cyclical permutation of (1,2,3).
!
            IF (I3 .EQ. 1) THEN
                I1 = 2
                I2 = 3
            ELSEIF (I3 .EQ. 2) THEN
                I1 = 3
                I2 = 1
            ELSE
                I1 = 1
                I2 = 2
            END IF
            N1 = LTRI(I1, KT1)
            N2 = LTRI(I2, KT1)
            N3 = LTRI(I3, KT1)
!
!   KT2 = (N2,N1,N4) for N4 = LTRI(I,KT2), where
!     LTRI(I+3,KT2) = KT1.
!
            IF (LTRI(4, KT2) .EQ. KT1) THEN
                I4 = 1
            ELSEIF (LTRI(5, KT2) .EQ. KT1) THEN
                I4 = 2
            ELSE
                I4 = 3
            END IF
            N4 = LTRI(I4, KT2)
!
!   The empty circumcircle test is reversed for the pseudo-
!     triangles.  The reversal is implicit in the clockwise
!     ordering of the vertices.
!
            IF (.NOT. SWPTST(N1, N2, N3, N4, X, Y, Z)) GO TO 5
!
!   Swap arc N1-N2 for N3-N4.  KTij is the triangle opposite
!     Nj as a vertex of KTi.
!
            SWP = .TRUE.
            KT11 = LTRI(I1 + 3, KT1)
            KT12 = LTRI(I2 + 3, KT1)
            IF (I4 .EQ. 1) THEN
                I2 = 2
                I1 = 3
            ELSEIF (I4 .EQ. 2) THEN
                I2 = 3
                I1 = 1
            ELSE
                I2 = 1
                I1 = 2
            END IF
            KT21 = LTRI(I1 + 3, KT2)
            KT22 = LTRI(I2 + 3, KT2)
            LTRI(1, KT1) = N4
            LTRI(2, KT1) = N3
            LTRI(3, KT1) = N1
            LTRI(4, KT1) = KT12
            LTRI(5, KT1) = KT22
            LTRI(6, KT1) = KT2
            LTRI(1, KT2) = N3
            LTRI(2, KT2) = N4
            LTRI(3, KT2) = N2
            LTRI(4, KT2) = KT21
            LTRI(5, KT2) = KT11
            LTRI(6, KT2) = KT1
!
!   Correct the KT11 and KT22 entries that changed.
!
            IF (KT11 .NE. 0) THEN
                I4 = 4
                IF (LTRI(4, KT11) .NE. KT1) THEN
                    I4 = 5
                    IF (LTRI(5, KT11) .NE. KT1) I4 = 6
                END IF
                LTRI(I4, KT11) = KT2
            END IF
            IF (KT22 .NE. 0) THEN
                I4 = 4
                IF (LTRI(4, KT22) .NE. KT2) THEN
                    I4 = 5
                    IF (LTRI(5, KT22) .NE. KT2) I4 = 6
                END IF
                LTRI(I4, KT22) = KT1
            END IF
5           CONTINUE
6           CONTINUE
            IF (SWP) GO TO 4
!
! Compute and store the negative circumcenters and radii of
!   the pseudo-triangles in the first NT positions.
!
7           DO 8 KT = 1, NT
                N1 = LTRI(1, KT)
                N2 = LTRI(2, KT)
                N3 = LTRI(3, KT)
                V1(1) = X(N1)
                V1(2) = Y(N1)
                V1(3) = Z(N1)
                V2(1) = X(N2)
                V2(2) = Y(N2)
                V2(3) = Z(N2)
                V3(1) = X(N3)
                V3(2) = Y(N3)
                V3(3) = Z(N3)
                CALL CIRCUM(V1, V2, V3, C, IERR)
                IF (IERR .NE. 0) GO TO 23
!
!   Store the negative circumcenter and radius (computed
!     from <V1,C>).
!
                XC(KT) = C(1)
                YC(KT) = C(2)
                ZC(KT) = C(3)
                T = V1(1) * C(1) + V1(2) * C(2) + V1(3) * C(3)
                IF (T .LT. -1.0) T = -1.0
                IF (T .GT. 1.0) T = 1.0
                RC(KT) = ACOS(T)
8               CONTINUE
!
! Compute and store the circumcenters and radii of the
!   actual triangles in positions KT = NT+1, NT+2, ...
!   Also, store the triangle indexes KT in the appropriate
!   LISTC positions.
!
9               KT = NT
!
!   Loop on nodes N1.
!
                NM2 = NN - 2
                DO 12 N1 = 1, NM2
                    LPL = LEND(N1)
                    LP = LPL
                    N3 = LIST(LP)
!
!   Loop on adjacent neighbors N2,N3 of N1 for which N2 > N1
!     and N3 > N1.
!
10                  LP = LPTR(LP)
                    N2 = N3
                    N3 = ABS(LIST(LP))
                    IF (N2 .LE. N1 .OR. N3 .LE. N1) GO TO 11
                    KT = KT + 1
!
!   Compute the circumcenter C of triangle KT = (N1,N2,N3).
!
                    V1(1) = X(N1)
                    V1(2) = Y(N1)
                    V1(3) = Z(N1)
                    V2(1) = X(N2)
                    V2(2) = Y(N2)
                    V2(3) = Z(N2)
                    V3(1) = X(N3)
                    V3(2) = Y(N3)
                    V3(3) = Z(N3)
                    CALL CIRCUM(V1, V2, V3, C, IERR)
                    IF (IERR .NE. 0) GO TO 23
!
!   Store the circumcenter, radius and triangle index.
!
                    XC(KT) = C(1)
                    YC(KT) = C(2)
                    ZC(KT) = C(3)
                    T = V1(1) * C(1) + V1(2) * C(2) + V1(3) * C(3)
                    IF (T .LT. -1.0) T = -1.0
                    IF (T .GT. 1.0) T = 1.0
                    RC(KT) = ACOS(T)
!
!   Store KT in LISTC(LPN), where Abs(LIST(LPN)) is the
!     index of N2 as a neighbor of N1, N3 as a neighbor
!     of N2, and N1 as a neighbor of N3.
!
                    LPN = LSTPTR(LPL, N2, LIST, LPTR)
                    LISTC(LPN) = KT
                    LPN = LSTPTR(LEND(N2), N3, LIST, LPTR)
                    LISTC(LPN) = KT
                    LPN = LSTPTR(LEND(N3), N1, LIST, LPTR)
                    LISTC(LPN) = KT
11                  IF (LP .NE. LPL) GO TO 10
12                  CONTINUE
                    IF (NT .EQ. 0) GO TO 20
!
! Store the first NT triangle indexes in LISTC.
!
!   Find a boundary triangle KT1 = (N1,N2,N3) with a
!     boundary arc opposite N3.
!
                    KT1 = 0
13                  KT1 = KT1 + 1
                    IF (LTRI(4, KT1) .EQ. 0) THEN
                        I1 = 2
                        I2 = 3
                        I3 = 1
                        GO TO 14
                    ELSEIF (LTRI(5, KT1) .EQ. 0) THEN
                        I1 = 3
                        I2 = 1
                        I3 = 2
                        GO TO 14
                    ELSEIF (LTRI(6, KT1) .EQ. 0) THEN
                        I1 = 1
                        I2 = 2
                        I3 = 3
                        GO TO 14
                    END IF
                    GO TO 13
14                  N1 = LTRI(I1, KT1)
                    N0 = N1
!
!   Loop on boundary nodes N1 in CCW order, storing the
!     indexes of the clockwise-ordered sequence of triangles
!     that contain N1.  The first triangle overwrites the
!     last neighbor position, and the remaining triangles,
!     if any, are appended to N1's adjacency list.
!
!   A pointer to the first neighbor of N1 is saved in LPN.
!
15                  LP = LEND(N1)
                    LPN = LPTR(LP)
                    LISTC(LP) = KT1
!
!   Loop on triangles KT2 containing N1.
!
16                  KT2 = LTRI(I2 + 3, KT1)
                    IF (KT2 .NE. 0) THEN
!
!   Append KT2 to N1's triangle list.
!
                        LPTR(LP) = LNEW
                        LP = LNEW
                        LISTC(LP) = KT2
                        LNEW = LNEW + 1
!
!   Set KT1 to KT2 and update (I1,I2,I3) such that
!     LTRI(I1,KT1) = N1.
!
                        KT1 = KT2
                        IF (LTRI(1, KT1) .EQ. N1) THEN
                            I1 = 1
                            I2 = 2
                            I3 = 3
                        ELSEIF (LTRI(2, KT1) .EQ. N1) THEN
                            I1 = 2
                            I2 = 3
                            I3 = 1
                        ELSE
                            I1 = 3
                            I2 = 1
                            I3 = 2
                        END IF
                        GO TO 16
                    END IF
!
!   Store the saved first-triangle pointer in LPTR(LP), set
!     N1 to the next boundary node, test for termination,
!     and permute the indexes:  the last triangle containing
!     a boundary node is the first triangle containing the
!     next boundary node.
!
                    LPTR(LP) = LPN
                    N1 = LTRI(I3, KT1)
                    IF (N1 .NE. N0) THEN
                        I4 = I3
                        I3 = I2
                        I2 = I1
                        I1 = I4
                        GO TO 15
                    END IF
!
! No errors encountered.
!
20                  IER = 0
                    RETURN
!
! N < 3.
!
21                  IER = 1
                    RETURN
!
! Insufficient space reserved for LTRI.
!
22                  IER = 2
                    RETURN
!
! Error flag returned by CIRCUM: KT indexes a null triangle.
!
23                  IER = 3
                    RETURN
                END
