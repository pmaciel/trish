INTEGER FUNCTION JRAND(N, IX, IY, IZ)
    INTEGER N, IX, IY, IZ
!
!***********************************************************
!   This function returns a uniformly distributed pseudo-
! random integer in the range 1 to N.
!
!
! On input:
!
!       N = Maximum value to be returned.
!
! N is not altered by this function.
!
!       IX,IY,IZ = Integer seeds initialized to values in
!                  the range 1 to 30,000 before the first
!                  call to JRAND, and not altered between
!                  subsequent calls (unless a sequence of
!                  random numbers is to be repeated by
!                  reinitializing the seeds).
!
! On output:
!
!       IX,IY,IZ = Updated integer seeds.
!
!       JRAND = Random integer in the range 1 to N.
!
! Reference:  B. A. Wichmann and I. D. Hill, "An Efficient
!             and Portable Pseudo-random Number Generator",
!             Applied Statistics, Vol. 31, No. 2, 1982,
!             pp. 188-190.
!
!***********************************************************
!
    REAL U, X
!
! Local parameters:
!
! U = Pseudo-random number uniformly distributed in the
!     interval (0,1).
! X = Pseudo-random number in the range 0 to 3 whose frac-
!       tional part is U.
!
    IX = MOD(171 * IX, 30269)
    IY = MOD(172 * IY, 30307)
    IZ = MOD(170 * IZ, 30323)
    X = (REAL(IX) / 30269.) + (REAL(IY) / 30307.) + (REAL(IZ) / 30323.)
    U = X - INT(X)
    JRAND = REAL(N) * U + 1.
    RETURN
END
