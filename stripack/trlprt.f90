SUBROUTINE TRLPRT(N, X, Y, Z, IFLAG, NROW, NT, LTRI, LOUT)
    INTEGER N, IFLAG, NROW, NT, LTRI(NROW, NT), LOUT
    REAL X(N), Y(N), Z(N)
!
!***********************************************************
!   This subroutine prints the triangle list created by Sub-
! routine TRLIST and, optionally, the nodal coordinates
! (either latitude and longitude or Cartesian coordinates)
! on logical unit LOUT.  The numbers of boundary nodes,
! triangles, and arcs are also printed.
!
!
! On input:
!
!       N = Number of nodes in the triangulation.
!           3 .LE. N .LE. 9999.
!
!       X,Y,Z = Arrays of length N containing the Cartesian
!               coordinates of the nodes if IFLAG = 0, or
!               (X and Y only) arrays of length N containing
!               longitude and latitude, respectively, if
!               IFLAG > 0, or unused dummy parameters if
!               IFLAG < 0.
!
!       IFLAG = Nodal coordinate option indicator:
!               IFLAG = 0 if X, Y, and Z (assumed to contain
!                         Cartesian coordinates) are to be
!                         printed (to 6 decimal places).
!               IFLAG > 0 if only X and Y (assumed to con-
!                         tain longitude and latitude) are
!                         to be printed (to 6 decimal
!                         places).
!               IFLAG < 0 if only the adjacency lists are to
!                         be printed.
!
!       NROW = Number of rows (entries per triangle) re-
!              served for the triangle list LTRI.  The value
!              must be 6 if only the vertex indexes and
!              neighboring triangle indexes are stored, or 9
!              if arc indexes are also stored.
!
!       NT = Number of triangles in the triangulation.
!            1 .LE. NT .LE. 9999.
!
!       LTRI = NROW by NT array whose J-th column contains
!              the vertex nodal indexes (first three rows),
!              neighboring triangle indexes (second three
!              rows), and, if NROW = 9, arc indexes (last
!              three rows) associated with triangle J for
!              J = 1,...,NT.
!
!       LOUT = Logical unit number for output.  If LOUT is
!              not in the range 0 to 99, output is written
!              to unit 6.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!   The triangle list and nodal coordinates (as specified by
! IFLAG) are written to unit LOUT.
!
! Modules required by TRLPRT:  None
!
!***********************************************************
!
    INTEGER I, K, LUN, NA, NB, NL, NLMAX, NMAX
    DATA NMAX/9999/, NLMAX/58/
!
! Local parameters:
!
! I =     DO-loop, nodal index, and row index for LTRI
! K =     DO-loop and triangle index
! LUN =   Logical unit number for output
! NA =    Number of triangulation arcs
! NB =    Number of boundary nodes
! NL =    Number of lines printed on the current page
! NLMAX = Maximum number of print lines per page (except
!           for the last page which may have two addi-
!           tional lines)
! NMAX =  Maximum value of N and NT (4-digit format)
!
    LUN = LOUT
    IF (LUN .LT. 0 .OR. LUN .GT. 99) LUN = 6
!
! Print a heading and test for invalid input.
!
    WRITE (LUN, 100) N
    NL = 3
    IF (N .LT. 3 .OR. N .GT. NMAX .OR. (NROW .NE. 6 .AND. NROW .NE. 9) .OR. NT .LT. 1 .OR. NT .GT. NMAX) THEN
!
! Print an error message and exit.
!
        WRITE (LUN, 110) N, NROW, NT
        RETURN
    END IF
    IF (IFLAG .EQ. 0) THEN
!
! Print X, Y, and Z.
!
        WRITE (LUN, 101)
        NL = 6
        DO 1 I = 1, N
        IF (NL .GE. NLMAX) THEN
            WRITE (LUN, 108)
            NL = 0
        END IF
        WRITE (LUN, 103) I, X(I), Y(I), Z(I)
        NL = NL + 1
1       CONTINUE
        ELSEIF (IFLAG .GT. 0) THEN
!
! Print X (longitude) and Y (latitude).
!
        WRITE (LUN, 102)
        NL = 6
        DO 2 I = 1, N
        IF (NL .GE. NLMAX) THEN
            WRITE (LUN, 108)
            NL = 0
        END IF
        WRITE (LUN, 104) I, X(I), Y(I)
        NL = NL + 1
2       CONTINUE
        END IF
!
! Print the triangulation LTRI.
!
        IF (NL .GT. NLMAX / 2) THEN
            WRITE (LUN, 108)
            NL = 0
        END IF
        IF (NROW .EQ. 6) THEN
            WRITE (LUN, 105)
        ELSE
            WRITE (LUN, 106)
        END IF
        NL = NL + 5
        DO 3 K = 1, NT
        IF (NL .GE. NLMAX) THEN
            WRITE (LUN, 108)
            NL = 0
        END IF
        WRITE (LUN, 107) K, (LTRI(I, K), I=1, NROW)
        NL = NL + 1
3       CONTINUE
!
! Print NB, NA, and NT (boundary nodes, arcs, and
!   triangles).
!
        NB = 2 * N - NT - 2
        IF (NB .LT. 3) THEN
            NB = 0
            NA = 3 * N - 6
        ELSE
            NA = NT + N - 1
        END IF
        WRITE (LUN, 109) NB, NA, NT
        RETURN
!
! Print formats:
!
100     FORMAT(///18X, 'STRIPACK (TRLIST) Output,  N = ', I4)
101     FORMAT(//8X, 'Node', 10X, 'X(Node)', 10X, 'Y(Node)', 10X, 'Z(Node)'//)
102     FORMAT(//16X, 'Node', 8X, 'Longitude', 9X, 'Latitude'//)
103     FORMAT(8X, I4, 3E17.6)
104     FORMAT(16X, I4, 2E17.6)
105     FORMAT(//1X, 'Triangle', 8X, 'Vertices', 12X, 'Neighbors'/4X, 'KT', 7X, 'N1', &
                5X, 'N2', 5X, 'N3', 4X, 'KT1', 4X, 'KT2', 4X, 'KT3'/)
106     FORMAT(//1X, 'Triangle', 8X, 'Vertices', 12X, 'Neighbors', 14X, 'Arcs'/4X, &
                'KT', 7X, 'N1', 5X, 'N2', 5X, 'N3', 4X, 'KT1', 4X, 'KT2', 4X, 'KT3', 4X, 'KA1', 4X, 'KA2', 4X, 'KA3'/)
107     FORMAT(2X, I4, 2X, 6(3X, I4), 3(2X, I5))
108     FORMAT(///)
109     FORMAT(/1X, 'NB = ', I4, ' Boundary Nodes', 5X, 'NA = ', I5, ' Arcs', 5X, 'NT = ', I5, ' Triangles')
110     FORMAT(//1X, 10X, '*** Invalid Parameter:  N =', I5, ', NROW =', I5, ', NT =', I5, ' ***')
        END
