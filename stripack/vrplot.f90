SUBROUTINE VRPLOT(LUN, PLTSIZ, ELAT, ELON, A, N, X, Y, Z, NT, LISTC, LPTR, LEND, XC, YC, ZC, TITLE, NUMBR, IER)
    CHARACTER * (*) TITLE
    INTEGER LUN, N, NT, LISTC(*), LPTR(*), LEND(N), IER
    LOGICAL NUMBR
    REAL PLTSIZ, ELAT, ELON, A, X(N), Y(N), Z(N), XC(NT), YC(NT), ZC(NT)
!
!***********************************************************
!   This subroutine creates a level-2 Encapsulated Post-
! script (EPS) file containing a graphical depiction of a
! Voronoi diagram of a set of nodes on the unit sphere.
! The visible vertices are projected onto the plane that
! contains the origin and has normal defined by a user-
! specified eye-position.  Projections of adjacent (visible)
! Voronoi vertices are connected by line segments.
!
!   The parameters defining the Voronoi diagram may be com-
! puted by Subroutine CRLIST.
!
!
! On input:
!
!       LUN = Logical unit number in the range 0 to 99.
!             The unit should be opened with an appropriate
!             file name before the call to this routine.
!
!       PLTSIZ = Plot size in inches.  A circular window in
!                the projection plane is mapped to a circu-
!                lar viewport with diameter equal to .88*
!                PLTSIZ (leaving room for labels outside the
!                viewport).  The viewport is centered on the
!                8.5 by 11 inch page, and its boundary is
!                drawn.  1.0 .LE. PLTSIZ .LE. 8.5.
!
!       ELAT,ELON = Latitude and longitude (in degrees) of
!                   the center of projection E (the center
!                   of the plot).  The projection plane is
!                   the plane that contains the origin and
!                   has E as unit normal.  In a rotated
!                   coordinate system for which E is the
!                   north pole, the projection plane con-
!                   tains the equator, and only northern
!                   hemisphere points are visible (from the
!                   point at infinity in the direction E).
!                   These are projected orthogonally onto
!                   the projection plane (by zeroing the z-
!                   component in the rotated coordinate
!                   system).  ELAT and ELON must be in the
!                   range -90 to 90 and -180 to 180, respec-
!                   tively.
!
!       A = Angular distance in degrees from E to the boun-
!           dary of a circular window against which the
!           Voronoi diagram is clipped.  The projected win-
!           dow is a disk of radius r = Sin(A) centered at
!           the origin, and only visible vertices whose
!           projections are within distance r of the origin
!           are included in the plot.  Thus, if A = 90, the
!           plot includes the entire hemisphere centered at
!           E.  0 .LT. A .LE. 90.
!
!       N = Number of nodes (Voronoi centers) and Voronoi
!           regions.  N .GE. 3.
!
!       X,Y,Z = Arrays of length N containing the Cartesian
!               coordinates of the nodes (unit vectors).
!
!       NT = Number of Voronoi region vertices (triangles,
!            including those in the extended triangulation
!            if the number of boundary nodes NB is nonzero):
!            NT = 2*N-4.
!
!       LISTC = Array of length 3*NT containing triangle
!               indexes (indexes to XC, YC, and ZC) stored
!               in 1-1 correspondence with LIST/LPTR entries
!               (or entries that would be stored in LIST for
!               the extended triangulation):  the index of
!               triangle (N1,N2,N3) is stored in LISTC(K),
!               LISTC(L), and LISTC(M), where LIST(K),
!               LIST(L), and LIST(M) are the indexes of N2
!               as a neighbor of N1, N3 as a neighbor of N2,
!               and N1 as a neighbor of N3.  The Voronoi
!               region associated with a node is defined by
!               the CCW-ordered sequence of circumcenters in
!               one-to-one correspondence with its adjacency
!               list (in the extended triangulation).
!
!       LPTR = Array of length 3*NT = 6*N-12 containing a
!              set of pointers (LISTC indexes) in one-to-one
!              correspondence with the elements of LISTC.
!              LISTC(LPTR(I)) indexes the triangle which
!              follows LISTC(I) in cyclical counterclockwise
!              order (the first neighbor follows the last
!              neighbor).
!
!       LEND = Array of length N containing a set of
!              pointers to triangle lists.  LP = LEND(K)
!              points to a triangle (indexed by LISTC(LP))
!              containing node K for K = 1 to N.
!
!       XC,YC,ZC = Arrays of length NT containing the
!                  Cartesian coordinates of the triangle
!                  circumcenters (Voronoi vertices).
!                  XC(I)**2 + YC(I)**2 + ZC(I)**2 = 1.
!
!       TITLE = Type CHARACTER variable or constant contain-
!               ing a string to be centered above the plot.
!               The string must be enclosed in parentheses;
!               i.e., the first and last characters must be
!               '(' and ')', respectively, but these are not
!               displayed.  TITLE may have at most 80 char-
!               acters including the parentheses.
!
!       NUMBR = Option indicator:  If NUMBR = TRUE, the
!               nodal indexes are plotted at the Voronoi
!               region centers.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if LUN, PLTSIZ, N, or NT is outside
!                     its valid range.
!             IER = 2 if ELAT, ELON, or A is outside its
!                     valid range.
!             IER = 3 if an error was encountered in writing
!                     to unit LUN.
!
!***********************************************************
!
    INTEGER IPX1, IPX2, IPY1, IPY2, IR, KV1, KV2, LP, LPL, N0
    LOGICAL ANNOT, IN1, IN2
    REAL CF, CT, EX, EY, EZ, FSIZN, FSIZT, R11, R12, R21, R22, R23, SF, T, TX, TY, WR, WRS, X0, X1, X2, Y0, Y1, Y2, Z1, Z2
!
    DATA ANNOT/.TRUE./, FSIZN/10.0/, FSIZT/16.0/
!
! Local parameters:
!
! ANNOT =     Logical variable with value TRUE iff the plot
!               is to be annotated with the values of ELAT,
!               ELON, and A
! CF =        Conversion factor for degrees to radians
! CT =        Cos(ELAT)
! EX,EY,EZ =  Cartesian coordinates of the eye-position E
! FSIZN =     Font size in points for labeling nodes with
!               their indexes if NUMBR = TRUE
! FSIZT =     Font size in points for the title (and
!               annotation if ANNOT = TRUE)
! IN1,IN2 =   Logical variables with value TRUE iff the
!               projections of vertices KV1 and KV2, respec-
!               tively, are inside the window
! IPX1,IPY1 = X and y coordinates (in points) of the lower
!               left corner of the bounding box or viewport
!               box
! IPX2,IPY2 = X and y coordinates (in points) of the upper
!               right corner of the bounding box or viewport
!               box
! IR =        Half the width (height) of the bounding box or
!               viewport box in points -- viewport radius
! KV1,KV2 =   Endpoint indexes of a Voronoi edge
! LP =        LIST index (pointer)
! LPL =       Pointer to the last neighbor of N0
! N0 =        Index of a node
! R11...R23 = Components of the first two rows of a rotation
!               that maps E to the north pole (0,0,1)
! SF =        Scale factor for mapping world coordinates
!               (window coordinates in [-WR,WR] X [-WR,WR])
!               to viewport coordinates in [IPX1,IPX2] X
!               [IPY1,IPY2]
! T =         Temporary variable
! TX,TY =     Translation vector for mapping world coordi-
!               nates to viewport coordinates
! WR =        Window radius r = Sin(A)
! WRS =       WR**2
! X0,Y0 =     Projection plane coordinates of node N0 or
!               label location
! X1,Y1,Z1 =  Coordinates of vertex KV1 in the rotated
!               coordinate system
! X2,Y2,Z2 =  Coordinates of vertex KV2 in the rotated
!               coordinate system or intersection of edge
!               KV1-KV2 with the equator (in the rotated
!               coordinate system)
!
!
! Test for invalid parameters.
!
    IF (LUN .LT. 0 .OR. LUN .GT. 99 .OR. PLTSIZ .LT. 1.0 .OR. PLTSIZ .GT. 8.5 .OR. N .LT. 3 .OR. NT .NE. 2 * N - 4) GO TO 11
    IF (ABS(ELAT) .GT. 90.0 .OR. ABS(ELON) .GT. 180.0 .OR. A .GT. 90.0) GO TO 12
!
! Compute a conversion factor CF for degrees to radians
!   and compute the window radius WR.
!
    CF = ATAN(1.0) / 45.0
    WR = SIN(CF * A)
    WRS = WR * WR
!
! Compute the lower left (IPX1,IPY1) and upper right
!   (IPX2,IPY2) corner coordinates of the bounding box.
!   The coordinates, specified in default user space units
!   (points, at 72 points/inch with origin at the lower
!   left corner of the page), are chosen to preserve the
!   square aspect ratio, and to center the plot on the 8.5
!   by 11 inch page.  The center of the page is (306,396),
!   and IR = PLTSIZ/2 in points.
!
    IR = NINT(36.0 * PLTSIZ)
    IPX1 = 306 - IR
    IPX2 = 306 + IR
    IPY1 = 396 - IR
    IPY2 = 396 + IR
!
! Output header comments.
!
    WRITE (LUN, 100, ERR=13) IPX1, IPY1, IPX2, IPY2
100 FORMAT('%!PS-Adobe-3.0 EPSF-3.0'/'%%BoundingBox:', 4I4/'%%Title:  Voronoi diagram'/'%%Creator:  STRIPACK'/'%%EndComments')
!
! Set (IPX1,IPY1) and (IPX2,IPY2) to the corner coordinates
!   of a viewport box obtained by shrinking the bounding box
!   by 12% in each dimension.
!
    IR = NINT(0.88 * REAL(IR))
    IPX1 = 306 - IR
    IPX2 = 306 + IR
    IPY1 = 396 - IR
    IPY2 = 396 + IR
!
! Set the line thickness to 2 points, and draw the
!   viewport boundary.
!
    T = 2.0
    WRITE (LUN, 110, ERR=13) T
    WRITE (LUN, 120, ERR=13) IR
    WRITE (LUN, 130, ERR=13)
110 FORMAT(F12.6, ' setlinewidth')
120 FORMAT('306 396 ', I3, ' 0 360 arc')
130 FORMAT('stroke')
!
! Set up an affine mapping from the window box [-WR,WR] X
!   [-WR,WR] to the viewport box.
!
    SF = REAL(IR) / WR
    TX = IPX1 + SF * WR
    TY = IPY1 + SF * WR
    WRITE (LUN, 140, ERR=13) TX, TY, SF, SF
140 FORMAT(2F12.6, ' translate'/2F12.6, ' scale')
!
! The line thickness must be changed to reflect the new
!   scaling which is applied to all subsequent output.
!   Set it to 1.0 point.
!
    T = 1.0 / SF
    WRITE (LUN, 110, ERR=13) T
!
! Save the current graphics state, and set the clip path to
!   the boundary of the window.
!
    WRITE (LUN, 150, ERR=13)
    WRITE (LUN, 160, ERR=13) WR
    WRITE (LUN, 170, ERR=13)
150 FORMAT('gsave')
160 FORMAT('0 0 ', F12.6, ' 0 360 arc')
170 FORMAT('clip newpath')
!
! Compute the Cartesian coordinates of E and the components
!   of a rotation R which maps E to the north pole (0,0,1).
!   R is taken to be a rotation about the z-axis (into the
!   yz-plane) followed by a rotation about the x-axis chosen
!   so that the view-up direction is (0,0,1), or (-1,0,0) if
!   E is the north or south pole.
!
!           ( R11  R12  0   )
!       R = ( R21  R22  R23 )
!           ( EX   EY   EZ  )
!
    T = CF * ELON
    CT = COS(CF * ELAT)
    EX = CT * COS(T)
    EY = CT * SIN(T)
    EZ = SIN(CF * ELAT)
    IF (CT .NE. 0.0) THEN
        R11 = -EY / CT
        R12 = EX / CT
    ELSE
        R11 = 0.0
        R12 = 1.0
    END IF
    R21 = -EZ * R12
    R22 = EZ * R11
    R23 = CT
!
! Loop on nodes (Voronoi centers) N0.
!   LPL indexes the last neighbor of N0.
!
    DO 3 N0 = 1, N
        LPL = LEND(N0)
!
! Set KV2 to the first (and last) vertex index and compute
!   its coordinates (X2,Y2,Z2) in the rotated coordinate
!   system.
!
        KV2 = LISTC(LPL)
        X2 = R11 * XC(KV2) + R12 * YC(KV2)
        Y2 = R21 * XC(KV2) + R22 * YC(KV2) + R23 * ZC(KV2)
        Z2 = EX * XC(KV2) + EY * YC(KV2) + EZ * ZC(KV2)
!
!   IN2 = TRUE iff KV2 is in the window.
!
        IN2 = Z2 .GE. 0. .AND. X2 * X2 + Y2 * Y2 .LE. WRS
!
! Loop on neighbors N1 of N0.  For each triangulation edge
!   N0-N1, KV1-KV2 is the corresponding Voronoi edge.
!
        LP = LPL
1       LP = LPTR(LP)
        KV1 = KV2
        X1 = X2
        Y1 = Y2
        Z1 = Z2
        IN1 = IN2
        KV2 = LISTC(LP)
!
!   Compute the new values of (X2,Y2,Z2) and IN2.
!
        X2 = R11 * XC(KV2) + R12 * YC(KV2)
        Y2 = R21 * XC(KV2) + R22 * YC(KV2) + R23 * ZC(KV2)
        Z2 = EX * XC(KV2) + EY * YC(KV2) + EZ * ZC(KV2)
        IN2 = Z2 .GE. 0. .AND. X2 * X2 + Y2 * Y2 .LE. WRS
!
! Add edge KV1-KV2 to the path iff both endpoints are inside
!   the window and KV2 > KV1, or KV1 is inside and KV2 is
!   outside (so that the edge is drawn only once).
!
        IF (.NOT. IN1 .OR. (IN2 .AND. KV2 .LE. KV1)) GO TO 2
        IF (Z2 .LT. 0.) THEN
!
!   KV2 is a 'southern hemisphere' point.  Move it to the
!     intersection of edge KV1-KV2 with the equator so that
!     the edge is clipped properly.  Z2 is implicitly set
!     to 0.
!
            X2 = Z1 * X2 - Z2 * X1
            Y2 = Z1 * Y2 - Z2 * Y1
            T = SQRT(X2 * X2 + Y2 * Y2)
            X2 = X2 / T
            Y2 = Y2 / T
        END IF
        WRITE (LUN, 180, ERR=13) X1, Y1, X2, Y2
180     FORMAT(2F12.6, ' moveto', 2F12.6, ' lineto')
!
! Bottom of loops.
!
2       IF (LP .NE. LPL) GO TO 1
3       CONTINUE
!
! Paint the path and restore the saved graphics state (with
!   no clip path).
!
        WRITE (LUN, 130, ERR=13)
        WRITE (LUN, 190, ERR=13)
190     FORMAT('grestore')
        IF (NUMBR) THEN
!
! Nodes in the window are to be labeled with their indexes.
!   Convert FSIZN from points to world coordinates, and
!   output the commands to select a font and scale it.
!
            T = FSIZN / SF
            WRITE (LUN, 200, ERR=13) T
200         FORMAT('/Helvetica findfont'/F12.6, ' scalefont setfont')
!
! Loop on visible nodes N0 that project to points (X0,Y0) in
!   the window.
!
            DO 4 N0 = 1, N
                IF (EX * X(N0) + EY * Y(N0) + EZ * Z(N0) .LT. 0.) GO TO 4
                X0 = R11 * X(N0) + R12 * Y(N0)
                Y0 = R21 * X(N0) + R22 * Y(N0) + R23 * Z(N0)
                IF (X0 * X0 + Y0 * Y0 .GT. WRS) GO TO 4
!
!   Move to (X0,Y0), and draw the label N0 with the origin
!     of the first character at (X0,Y0).
!
                WRITE (LUN, 210, ERR=13) X0, Y0
                WRITE (LUN, 220, ERR=13) N0
210             FORMAT(2F12.6, ' moveto')
220             FORMAT('(', I3, ') show')
4               CONTINUE
                END IF
!
! Convert FSIZT from points to world coordinates, and output
!   the commands to select a font and scale it.
!
                T = FSIZT / SF
                WRITE (LUN, 200, ERR=13) T
!
! Display TITLE centered above the plot:
!
                Y0 = WR + 3.0 * T
                WRITE (LUN, 230, ERR=13) TITLE, Y0
230             FORMAT(A80/'  stringwidth pop 2 div neg ', F12.6, ' moveto')
                WRITE (LUN, 240, ERR=13) TITLE
240             FORMAT(A80/'  show')
                IF (ANNOT) THEN
!
! Display the window center and radius below the plot.
!
                    X0 = -WR
                    Y0 = -WR - 50.0 / SF
                    WRITE (LUN, 210, ERR=13) X0, Y0
                    WRITE (LUN, 250, ERR=13) ELAT, ELON
                    Y0 = Y0 - 2.0 * T
                    WRITE (LUN, 210, ERR=13) X0, Y0
                    WRITE (LUN, 260, ERR=13) A
250                 FORMAT('(Window center:  ELAT = ', F7.2, ',  ELON = ', F8.2, ') show')
260                 FORMAT('(Angular extent:  A = ', F5.2, ') show')
                END IF
!
! Paint the path and output the showpage command and
!   end-of-file indicator.
!
                WRITE (LUN, 270, ERR=13)
270             FORMAT('stroke'/'showpage'/'%%EOF')
!
! HP's interpreters require a one-byte End-of-PostScript-Job
!   indicator (to eliminate a timeout error message):
!   ASCII 4.
!
                WRITE (LUN, 280, ERR=13) CHAR(4)
280             FORMAT(A1)
!
! No error encountered.
!
                IER = 0
                RETURN
!
! Invalid input parameter LUN, PLTSIZ, N, or NT.
!
11              IER = 1
                RETURN
!
! Invalid input parameter ELAT, ELON, or A.
!
12              IER = 2
                RETURN
!
! Error writing to unit LUN.
!
13              IER = 3
                RETURN
            END
