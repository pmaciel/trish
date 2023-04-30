!--**--CH465--772--P:RW--21:7:1999
!--**--CH464--772--C:F--21:7:1999
!
!
!             STRITEST:  STRIPACK Test Program
!                         07/30/98
!
!
!   This driver tests software package STRIPACK for con-
! structing a Delaunay triangulation and Voronoi diagram
! of a set of nodes on the surface of the unit sphere.
! All STRIPACK subprograms are tested.
!
!   By default, a triangulation is created from a set of N
! nodes consisting of the north pole and N-1 points uniform-
! ly distributed around the 60-degree parallel (with
! constant longitudinal separation).  However, by enabling
! the READ statements below (C# in the first two columns),
! testing may be performed on an arbitrary set of up to NMAX
! nodes (refer to the PARAMETER statement below).  A data
! set consists of the following sequence of records:
!
!    N = Number of nodes (3 to NMAX):  format I5.
!    (RLAT(I),RLON(I), I = 1,N) = Nodal coordinates in
!                 degrees latitude (-90 to 90) and degrees
!                 longitude (-180 to 180):  format 2F15.9.
!
! This program must be linked to STRIPACK.
!
!
CHARACTER * 80 TITLET, TITLEV
INTEGER IER, IFLAG, K, KSUM, KT, LIN, LOUT, LP, LPL, LPLT, LPLV, LW, LWK, LNEW, N, N0, N1, N2, N3, NA, NB, NCOL, NMAX, NN
INTEGER NROW, NT, NT6, NTMX, NV
INTEGER NEARND
LOGICAL INSIDE, NUMBR
REAL A, AL, AREA, DEL, ELAT, ELON, P(3), PLTSIZ, SC, V1(3), V2(3), V3(3), VLAT, VLON, VNRM
REAL AREAS
!
PARAMETER(NMAX=100, NTMX=2 * NMAX, NT6=6 * NMAX, LWK=2 * NMAX, NCOL=NMAX, NROW=9)
PARAMETER(LIN=5, LOUT=6)
!
! Array storage for the triangulation, work space, and nodal
!   coordinates.
!
INTEGER LIST(NT6), LPTR(NT6), LEND(NMAX), IWK(LWK)
REAL DS(NMAX), RLAT(NMAX), RLON(NMAX), X(NMAX), Y(NMAX), Z(NMAX)
!
! Array storage for the Voronoi diagram:  adjacency array,
!   boundary triangle list, triangle circumcenters, and
!   circumradii.
!
INTEGER LISTC(NT6), LBTRI(6, NCOL)
REAL XC(NTMX), YC(NTMX), ZC(NTMX), RC(NTMX)
!
! Array storage for the triangle list.
!
INTEGER LTRI(NROW, NTMX)
!
! Plot size.
!
DATA PLTSIZ/7.5/
!
! Logical unit numbers for I/O:
!
DATA LPLT/3/, LPLV/4/
OPEN (LOUT, FILE='res')
OPEN (LPLT, FILE='res_3.eps')
OPEN (LPLV, FILE='res_4.eps')
!
! Store plot titles  They must be enclosed in parentheses.
!
TITLET = '(Triangulation created by STRITEST)'
TITLEV = '(Voronoi diagram created by STRITEST)'
!
! Generate the default set of nodes as latitudinal and lon-
!   gitudinal coordinates.  DEL is the separation in degrees
!   between the nodes on the 60-degree line of latitude.
!
N = 9
RLAT(1) = 90.
RLON(1) = 0.
DEL = 360./REAL(N - 1)
DO 1 K = 2, N
    RLAT(K) = 60.
    RLON(K) = REAL(K - 2) * DEL
1   CONTINUE
!
! *** Read a data set:  N, RLAT, and RLON.
!
!#    OPEN (LIN,FILE='stritest.dat',STATUS='OLD')
!#    READ (LIN,100,ERR=30) N
    IF (N .LT. 3 .OR. N .GT. NMAX) GO TO 31
!#    READ (LIN,110,ERR=30) (RLAT(K),RLON(K), K = 1,N)
!#100 FORMAT (I5)
!#110 FORMAT (2F15.9)
!
! Print a heading.
!
    WRITE (LOUT, 300)
300 FORMAT(///1X, 25X, 'STRIPACK test'///)
!
! Set X and Y to the values of RLON and RLAT, respectively,
!   in radians.  (RLON and RLAT are saved for printing by
!   Subroutine TRPRNT).
!
    SC = ATAN(1.) / 45.
    DO 2 K = 1, N
        X(K) = SC * RLON(K)
        Y(K) = SC * RLAT(K)
2       CONTINUE
!
! *** Transform spherical coordinates X and Y to Cartesian
!       coordinates (X,Y,Z) on the unit sphere (X**2 +
!       Y**2 + Z**2 = 1).
!
        CALL TRANS(N, Y, X, X, Y, Z)
!
! *** Create the triangulation and test the error flag.
!
        CALL TRMESH(N, X, Y, Z, LIST, LPTR, LEND, LNEW, IWK, IWK(N + 1), DS, IER)
        IF (IER .EQ. -2) THEN
            WRITE (LOUT, 510)
            STOP
        ELSEIF (IER .GT. 0) THEN
            WRITE (LOUT, 515)
            STOP
        END IF
!
! *** Print the spherical coordinates and adjacency informa-
!       tion on LOUT.  IFLAG > 0 indicates that RLON and
!       RLAT only are to be printed.
!
        IFLAG = 1
        CALL TRPRNT(N, RLON, RLAT, Z, IFLAG, LIST, LPTR, LEND, LOUT)
!
! *** Test TRLIST and TRLPRT by creating and printing a
!                 triangle list.
!
        CALL TRLIST(N, LIST, LPTR, LEND, NROW, NT, LTRI, IER)
        CALL TRLPRT(N, RLON, RLAT, Z, IFLAG, NROW, NT, LTRI, LOUT)
!
! *** Test TRPLOT by plotting the portion of the triangula-
!                 tion contained in the hemisphere centered
!                 at E = (ELAT,ELON), where ELAT and ELON
!                 are taken to be the center of the range of
!                 the nodal latitudes and longitudes.
!
        ELAT = RLAT(1)
        VLAT = ELAT
        ELON = RLON(1)
        VLON = ELON
        DO 3 K = 2, N
            IF (RLAT(K) .LT. ELAT) ELAT = RLAT(K)
            IF (RLAT(K) .GT. VLAT) VLAT = RLAT(K)
            IF (RLON(K) .LT. ELON) ELON = RLON(K)
            IF (RLON(K) .GT. VLON) VLON = RLON(K)
3           CONTINUE
            ELAT = (ELAT + VLAT) / 2.0
            ELON = (ELON + VLON) / 2.0
            A = 90.0
            NUMBR = N .LE. 200
            CALL TRPLOT(LPLT, PLTSIZ, ELAT, ELON, A, N, X, Y, Z, LIST, LPTR, LEND, TITLET, NUMBR, IER)
            IF (IER .EQ. 0) THEN
                WRITE (LOUT, 305)
305             FORMAT(/5X, 'A triangulation plot file was ', 'successfully created.'/)
            ELSE
                WRITE (LOUT, 590) IER
            END IF
!
! *** Test AREAS by computing and printing the area of the
!                convex hull of the nodes (sum of triangle
!                areas) relative to the total surface area
!                (4*Pi).
!
            AREA = 0.
            DO 4 KT = 1, NT
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
                AREA = AREA + AREAS(V1, V2, V3)
4               CONTINUE
                AREA = AREA / (16.0 * ATAN(1.0))
                WRITE (LOUT, 310) AREA
310             FORMAT(//5X, 'Area of convex hull relative to total ', 'surface area = ', F8.2)
!
! *** Test BNODES.  The ordered sequence of boundary nodes
!                   is stored in IWK.
!
                CALL BNODES(N, LIST, LPTR, LEND, IWK, NB, NA, NT)
                WRITE (LOUT, 320) NB, NA, NT
320             FORMAT(//5X, 'Output from BNODES:  ', I4, ' boundary nodes,  ', I4, ' arcs,  ', I4, ' triangles.'//)
!
! *** Test GETNP by ordering the nodes on distance from N0
!                and verifying the ordering.  The sequence
!                of nodal indexes is stored in IWK, and
!                the values of an increasing function (the
!                negative cosine) of angular distance is
!                stored in DS.
!
                N0 = N / 2
                IWK(1) = N0
                DS(1) = -1.0
                KSUM = N0
                DO 5 K = 2, N
                    CALL GETNP(X, Y, Z, LIST, LPTR, LEND, K, IWK, DS(K), IER)
                    IF (IER .NE. 0 .OR. DS(K) .LT. DS(K - 1)) THEN
                        WRITE (LOUT, 520)
                        STOP
                    END IF
                    KSUM = KSUM + IWK(K)
5                   CONTINUE
!
!   Test for all nodal indexes included in IWK.
!
                    IF (KSUM .NE. (N * (N + 1)) / 2) THEN
                        WRITE (LOUT, 520)
                        STOP
                    END IF
!
! *** Test NEARND by verifying that the nearest node to K is
!                 node K for K = 1 to N.
!
                    DO 6 K = 1, N
                        P(1) = X(K)
                        P(2) = Y(K)
                        P(3) = Z(K)
                        N0 = NEARND(P, 1, N, X, Y, Z, LIST, LPTR, LEND, AL)
                        IF (N0 .NE. K .OR. AL .GT. 1.E-3) THEN
                            WRITE (LOUT, 530)
                            STOP
                        END IF
6                       CONTINUE
!
! *** Test DELARC by removing a boundary arc if possible.
!                 The last two nodes define a boundary arc
!                 in the default data set.
!
                        N1 = N - 1
                        N2 = N
                        CALL DELARC(N, N1, N2, LIST, LPTR, LEND, LNEW, IER)
                        IF (IER .EQ. 1 .OR. IER .EQ. 4) THEN
                            WRITE (LOUT, 540) IER
                            STOP
                        END IF
                        IF (IER .NE. 0) THEN
                            WRITE (LOUT, 330) N1, N2
330                         FORMAT(5X, 'Subroutine DELARC not tested:'/5X, 'Nodes ', I4, ' and ', I4, &
                                   ' do not form a ', 'removable boundary arc.'//)
                        ELSE
                            CALL TRMESH(N, X, Y, Z, LIST, LPTR, LEND, LNEW, IWK, IWK(N + 1), DS, IER)
                        END IF
!
! *** Test CRLIST, VRPLOT, and SCOORD by constructing and
!                 plotting the Voronoi diagram, and printing
!                 the Voronoi region boundary (ordered
!                 sequence of Voronoi vertices) associated
!                 with N0.
!
!     Note that the triangulation data structure
!       is altered if NB > 0.
!
                        CALL CRLIST(N, NCOL, X, Y, Z, LIST, LEND, LPTR, LNEW, LBTRI, LISTC, NB, XC, YC, ZC, RC, IER)
                        IF (IER .NE. 0) THEN
                            WRITE (LOUT, 550) IER
                            STOP
                        END IF
!
!   Use the same parameter values that were used for the
!     triangulation plot (except the output unit and title).
!
                        NT = 2 * N - 4
                        CALL VRPLOT(LPLV, PLTSIZ, ELAT, ELON, A, N, X, Y, Z, NT, LISTC, LPTR, LEND, XC, YC, ZC, TITLEV, NUMBR, IER)
                        IF (IER .EQ. 0) THEN
                            WRITE (LOUT, 335)
335                         FORMAT(/5X, 'A Voronoi diagram plot file was ', 'successfully created.'/)
                        ELSE
                            WRITE (LOUT, 600) IER
                        END IF
                        N0 = 1
                        WRITE (LOUT, 340) N0
340                     FORMAT(17X, 'Voronoi region for node ', I4//5X, 'Triangle', 5X, 'Latitude', 5X, 'Longitude', 5X, &
                               'Circumradius'/)
!
!   Initialize for loop on Voronoi vertices (triangle cir-
!     cumcenters).  The number of vertices is accumulated
!     in NV, and the vertex indexes are stored in IWK.  The
!     vertices are converted to latitude and longitude in
!     degrees for printing.
!
                        NV = 0
                        LPL = LEND(N0)
                        LP = LPL
7                       LP = LPTR(LP)
                        KT = LISTC(LP)
                        NV = NV + 1
                        IWK(NV) = KT
                        CALL SCOORD(XC(KT), YC(KT), ZC(KT), VLAT, VLON, VNRM)
                        VLAT = VLAT / SC
                        VLON = VLON / SC
                        WRITE (LOUT, 345) KT, VLAT, VLON, RC(KT)
                        IF (LP .NE. LPL) GO TO 7
345                     FORMAT(9X, I4, 1X, F12.6, 2X, F12.6, 5X, F12.6)
!
! *** Test INSIDE by checking for node N0 inside its
!                 Voronoi region.
!
                        P(1) = X(N0)
                        P(2) = Y(N0)
                        P(3) = Z(N0)
                        IF (.NOT. INSIDE(P, NTMX, XC, YC, ZC, NV, IWK, IER)) WRITE (LOUT, 560) N0
                        IF (IER .NE. 0) WRITE (LOUT, 565) IER
!
! *** Recreate the triangulation and test the error flag.
!
                        CALL TRMESH(N, X, Y, Z, LIST, LPTR, LEND, LNEW, IWK, IWK(N + 1), DS, IER)
                        IF (IER .EQ. -2) THEN
                            WRITE (LOUT, 510)
                            STOP
                        ELSEIF (IER .GT. 0) THEN
                            WRITE (LOUT, 515)
                            STOP
                        END IF
!
! *** Test EDGE by forcing an edge between nodes N1=1 and
!               N2=N.  LW is the number of columns reserved
!               for a 2 by LW work space array (IWK).
!
                        N1 = 1
                        N2 = N
                        LW = LWK / 2
                        CALL EDGE(N1, N2, X, Y, Z, LW, IWK, LIST, LPTR, LEND, IER)
                        IF (IER .NE. 0 .AND. IER .NE. 5) THEN
                            WRITE (LOUT, 570) IER
                            STOP
                        END IF
!
! *** Test DELNOD by removing nodes 4 to N (in reverse
!                 order).  LW is the number of columns
!                 reserved for a 2 by LW work space array
!                 (IWK).
!
                        IF (N .LE. 3) THEN
                            WRITE (LOUT, 350)
350                         FORMAT(//5X, 'Subroutine DELNOD not tested:  ', 'N cannot be reduced below 3.'//)
                        ELSE
                            NN = N
8                           LW = LWK / 2
                            CALL DELNOD(NN, NN, X, Y, Z, LIST, LPTR, LEND, LNEW, LW, IWK, IER)
                            IF (IER .NE. 0 .AND. IER .NE. 5) THEN
                                WRITE (LOUT, 580) IER
                                STOP
                            END IF
                            IF (NN .GT. 3) GO TO 8
                        END IF
!
! Successful test.
!
                        WRITE (LOUT, 360)
360                     FORMAT(//5X, 'No errors encountered.'/)
                        STOP
!
! Error reading the data set.
!
!#   30 WRITE (*,500)
!#      STOP
!
! Invalid value of N.
!
31                      WRITE (*, 505) N
                        STOP
!
! Error message formats:
!
!#500 FORMAT (//5X,'*** Input data set invalid ***'/)
505                     FORMAT(//5X, '*** N is outside its valid ', 'range:  N =', I5, ' ***'/)
510                     FORMAT(//5X, '*** Error in TRMESH:  the first three ', 'nodes are collinear ***'/)
515                     FORMAT(//5X, '*** Error in TRMESH:  duplicate nodes ', 'encountered ***'/)
520                     FORMAT(//5X, '*** Error in GETNP ***'/)
530                     FORMAT(//5X, '*** Error in NEARND ***'/)
540                     FORMAT(//5X, '*** Error in DELARC:  IER = ', I1, ' ***'/)
550                     FORMAT(//5X, '*** Error in CRLIST:  IER = ', I1, ' ***'/)
560                     FORMAT(//5X, '*** Error in INSIDE:  node ', I4, ' is ', 'not contained in its Voronoi region ***'/)
565                     FORMAT(//5X, '*** Error in INSIDE:  IER = ', I1, ' ***'/)
570                     FORMAT(//5X, '*** Error in EDGE:  IER = ', I1, ' ***'/)
580                     FORMAT(//5X, '*** Error in DELNOD:  IER = ', I1, ' ***'/)
590                     FORMAT(//5X, '*** Error in TRPLOT:  IER = ', I1, ' ***'/)
600                     FORMAT(//5X, '*** Error in VRPLOT:  IER = ', I1, ' ***'/)
                    END
