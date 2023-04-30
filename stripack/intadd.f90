SUBROUTINE INTADD(KK, I1, I2, I3, LIST, LPTR, LEND, LNEW)
    INTEGER KK, I1, I2, I3, LIST(*), LPTR(*), LEND(*), LNEW
!
!***********************************************************
!   This subroutine adds an interior node to a triangulation
! of a set of points on the unit sphere.  The data structure
! is updated with the insertion of node KK into the triangle
! whose vertices are I1, I2, and I3.  No optimization of the
! triangulation is performed.
!
!   This routine is identical to the similarly named routine
! in TRIPACK.
!
!
! On input:
!
!       KK = Index of the node to be inserted.  KK .GE. 1
!            and KK must not be equal to I1, I2, or I3.
!
!       I1,I2,I3 = Indexes of the counterclockwise-ordered
!                  sequence of vertices of a triangle which
!                  contains node KK.
!
! The above parameters are not altered by this routine.
!
!       LIST,LPTR,LEND,LNEW = Data structure defining the
!                             triangulation.  Refer to Sub-
!                             routine TRMESH.  Triangle
!                             (I1,I2,I3) must be included
!                             in the triangulation.
!
! On output:
!
!       LIST,LPTR,LEND,LNEW = Data structure updated with
!                             the addition of node KK.  KK
!                             will be connected to nodes I1,
!                             I2, and I3.
!
! Modules required by INTADD:  INSERT, LSTPTR
!
!***********************************************************
!
    INTEGER LSTPTR
    INTEGER K, LP, N1, N2, N3
!
! Local parameters:
!
! K =        Local copy of KK
! LP =       LIST pointer
! N1,N2,N3 = Local copies of I1, I2, and I3
!
    K = KK
!
! Initialization.
!
    N1 = I1
    N2 = I2
    N3 = I3
!
! Add K as a neighbor of I1, I2, and I3.
!
    LP = LSTPTR(LEND(N1), N2, LIST, LPTR)
    CALL INSERT(K, LP, LIST, LPTR, LNEW)
    LP = LSTPTR(LEND(N2), N3, LIST, LPTR)
    CALL INSERT(K, LP, LIST, LPTR, LNEW)
    LP = LSTPTR(LEND(N3), N1, LIST, LPTR)
    CALL INSERT(K, LP, LIST, LPTR, LNEW)
!
! Add I1, I2, and I3 as neighbors of K.
!
    LIST(LNEW) = N1
    LIST(LNEW + 1) = N2
    LIST(LNEW + 2) = N3
    LPTR(LNEW) = LNEW + 1
    LPTR(LNEW + 1) = LNEW + 2
    LPTR(LNEW + 2) = LNEW
    LEND(K) = LNEW + 2
    LNEW = LNEW + 3
    RETURN
END
