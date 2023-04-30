SUBROUTINE TRPRNT(N, X, Y, Z, IFLAG, LIST, LPTR, LEND, LOUT)
    INTEGER N, IFLAG, LIST(*), LPTR(*), LEND(N), LOUT
    REAL X(N), Y(N), Z(N)
!
!***********************************************************
!   This subroutine prints the triangulation adjacency lists
! created by Subroutine TRMESH and, optionally, the nodal
! coordinates (either latitude and longitude or Cartesian
! coordinates) on logical unit LOUT.  The list of neighbors
! of a boundary node is followed by index 0.  The numbers of
! boundary nodes, triangles, and arcs are also printed.
!
!
! On input:
!
!       N = Number of nodes in the triangulation.  N .GE. 3
!           and N .LE. 9999.
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
!       LIST,LPTR,LEND = Data structure defining the trian-
!                        gulation.  Refer to Subroutine
!                        TRMESH.
!
!       LOUT = Logical unit for output.  If LOUT is not in
!              the range 0 to 99, output is written to
!              logical unit 6.
!
! Input parameters are not altered by this routine.
!
! On output:
!
!   The adjacency lists and nodal coordinates (as specified
! by IFLAG) are written to unit LOUT.
!
! Modules required by TRPRNT:  None
!
!***********************************************************
!
    INTEGER I, INC, K, LP, LPL, LUN, NA, NABOR(400), NB, ND, NL, NLMAX, NMAX, NODE, NN, NT
    DATA NMAX/9999/, NLMAX/58/
!
! Local parameters:
!
! I =     NABOR index (1 to K)
! INC =   Increment for NL associated with an adjacency list
! K =     Counter and number of neighbors of NODE
! LP =    LIST pointer of a neighbor of NODE
! LPL =   Pointer to the last neighbor of NODE
! LUN =   Logical unit for output (copy of LOUT)
! NA =    Number of arcs in the triangulation
! NABOR = Array containing the adjacency list associated
!           with NODE, with zero appended if NODE is a
!           boundary node
! NB =    Number of boundary nodes encountered
! ND =    Index of a neighbor of NODE (or negative index)
! NL =    Number of lines that have been printed on the
!           current page
! NLMAX = Maximum number of print lines per page (except
!           for the last page which may have two addi-
!           tional lines)
! NMAX =  Upper bound on N (allows 4-digit indexes)
! NODE =  Index of a node and DO-loop index (1 to N)
! NN =    Local copy of N
! NT =    Number of triangles in the triangulation
!
    NN = N
    LUN = LOUT
    IF (LUN .LT. 0 .OR. LUN .GT. 99) LUN = 6
!
! Print a heading and test the range of N.
!
    WRITE (LUN, 100) NN
    IF (NN .LT. 3 .OR. NN .GT. NMAX) THEN
!
! N is outside its valid range.
!
        WRITE (LUN, 110)
        RETURN
    END IF
!
! Initialize NL (the number of lines printed on the current
!   page) and NB (the number of boundary nodes encountered).
!
    NL = 6
    NB = 0
    IF (IFLAG .LT. 0) THEN
!
! Print LIST only.  K is the number of neighbors of NODE
!   that have been stored in NABOR.
!
        WRITE (LUN, 101)
        DO 2 NODE = 1, NN
            LPL = LEND(NODE)
            LP = LPL
            K = 0
!
1           K = K + 1
            LP = LPTR(LP)
            ND = LIST(LP)
            NABOR(K) = ND
            IF (LP .NE. LPL) GO TO 1
            IF (ND .LE. 0) THEN
!
!   NODE is a boundary node.  Correct the sign of the last
!     neighbor, add 0 to the end of the list, and increment
!     NB.
!
                NABOR(K) = -ND
                K = K + 1
                NABOR(K) = 0
                NB = NB + 1
            END IF
!
!   Increment NL and print the list of neighbors.
!
            INC = (K - 1) / 14 + 2
            NL = NL + INC
            IF (NL .GT. NLMAX) THEN
                WRITE (LUN, 108)
                NL = INC
            END IF
            WRITE (LUN, 104) NODE, (NABOR(I), I=1, K)
            IF (K .NE. 14) WRITE (LUN, 107)
2           CONTINUE
            ELSEIF (IFLAG .GT. 0) THEN
!
! Print X (longitude), Y (latitude), and LIST.
!
            WRITE (LUN, 102)
            DO 4 NODE = 1, NN
                LPL = LEND(NODE)
                LP = LPL
                K = 0
!
3               K = K + 1
                LP = LPTR(LP)
                ND = LIST(LP)
                NABOR(K) = ND
                IF (LP .NE. LPL) GO TO 3
                IF (ND .LE. 0) THEN
!
!   NODE is a boundary node.
!
                    NABOR(K) = -ND
                    K = K + 1
                    NABOR(K) = 0
                    NB = NB + 1
                END IF
!
!   Increment NL and print X, Y, and NABOR.
!
                INC = (K - 1) / 8 + 2
                NL = NL + INC
                IF (NL .GT. NLMAX) THEN
                    WRITE (LUN, 108)
                    NL = INC
                END IF
                WRITE (LUN, 105) NODE, X(NODE), Y(NODE), (NABOR(I), I=1, K)
                IF (K .NE. 8) WRITE (LUN, 107)
4               CONTINUE
                ELSE
!
! Print X, Y, Z, and LIST.
!
                WRITE (LUN, 103)
                DO 6 NODE = 1, NN
                    LPL = LEND(NODE)
                    LP = LPL
                    K = 0
!
5                   K = K + 1
                    LP = LPTR(LP)
                    ND = LIST(LP)
                    NABOR(K) = ND
                    IF (LP .NE. LPL) GO TO 5
                    IF (ND .LE. 0) THEN
!
!   NODE is a boundary node.
!
                        NABOR(K) = -ND
                        K = K + 1
                        NABOR(K) = 0
                        NB = NB + 1
                    END IF
!
!   Increment NL and print X, Y, Z, and NABOR.
!
                    INC = (K - 1) / 5 + 2
                    NL = NL + INC
                    IF (NL .GT. NLMAX) THEN
                        WRITE (LUN, 108)
                        NL = INC
                    END IF
                    WRITE (LUN, 106) NODE, X(NODE), Y(NODE), Z(NODE), (NABOR(I), I=1, K)
                    IF (K .NE. 5) WRITE (LUN, 107)
6                   CONTINUE
                    END IF
!
! Print NB, NA, and NT (boundary nodes, arcs, and
!   triangles).
!
                    IF (NB .NE. 0) THEN
                        NA = 3 * NN - NB - 3
                        NT = 2 * NN - NB - 2
                    ELSE
                        NA = 3 * NN - 6
                        NT = 2 * NN - 4
                    END IF
                    WRITE (LUN, 109) NB, NA, NT
                    RETURN
!
! Print formats:
!
100                 FORMAT(///15X, 'STRIPACK Triangulation Data ', 'Structure,  N = ', I5//)
101                 FORMAT(1X, 'Node', 31X, 'Neighbors of Node'//)
102                 FORMAT(1X, 'Node', 5X, 'Longitude', 6X, 'Latitude', 18X, 'Neighbors of Node'//)
103                 FORMAT(1X, 'Node', 5X, 'X(Node)', 8X, 'Y(Node)', 8X, 'Z(Node)', 11X, 'Neighbors of Node'//)
104                 FORMAT(1X, I4, 4X, 14I5 / (1X, 8X, 14I5))
105                 FORMAT(1X, I4, 2E15.6, 4X, 8I5 / (1X, 38X, 8I5))
106                 FORMAT(1X, I4, 3E15.6, 4X, 5I5 / (1X, 53X, 5I5))
107                 FORMAT(1X)
108                 FORMAT(///)
109                 FORMAT(/1X, 'NB = ', I4, ' Boundary Nodes', 5X, 'NA = ', I5, ' Arcs', 5X, 'NT = ', I5, ' Triangles')
110                 FORMAT(1X, 10X, '*** N is outside its valid', ' range ***')
                END
