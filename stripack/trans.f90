SUBROUTINE TRANS(N, RLAT, RLON, X, Y, Z)
    INTEGER N
    REAL RLAT(N), RLON(N), X(N), Y(N), Z(N)
!
!***********************************************************
!   This subroutine transforms spherical coordinates into
! Cartesian coordinates on the unit sphere for input to
! Subroutine TRMESH.  Storage for X and Y may coincide with
! storage for RLAT and RLON if the latter need not be saved.
!
!
! On input:
!
!       N = Number of nodes (points on the unit sphere)
!           whose coordinates are to be transformed.
!
!       RLAT = Array of length N containing latitudinal
!              coordinates of the nodes in radians.
!
!       RLON = Array of length N containing longitudinal
!              coordinates of the nodes in radians.
!
! The above parameters are not altered by this routine.
!
!       X,Y,Z = Arrays of length at least N.
!
! On output:
!
!       X,Y,Z = Cartesian coordinates in the range -1 to 1.
!               X(I)**2 + Y(I)**2 + Z(I)**2 = 1 for I = 1
!               to N.
!
!***********************************************************
!
    INTEGER I, NN
    REAL COSPHI, PHI, THETA
!
! Local parameters:
!
! COSPHI = cos(PHI)
! I =      DO-loop index
! NN =     Local copy of N
! PHI =    Latitude
! THETA =  Longitude
!
    NN = N
    DO 1 I = 1, NN
        PHI = RLAT(I)
        THETA = RLON(I)
        COSPHI = COS(PHI)
        X(I) = COSPHI * COS(THETA)
        Y(I) = COSPHI * SIN(THETA)
        Z(I) = SIN(PHI)
1       CONTINUE
        RETURN
    END
