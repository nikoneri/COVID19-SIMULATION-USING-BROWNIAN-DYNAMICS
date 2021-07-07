SUBROUTINE CONFIGINIALE(X,Y,N,BoxL)
!*********************************************************************************************************************
!   1. AUTHORS AND E-MAIL CONTACT:
!   -Isaac Neri Gomez Sarmiento (isaacneri.gs@gmail.com)
!   -César Omar Ramírez Álvarez (cesaromarramirezalvarez@gmail.com)
!   -Jonas Valenzuela Terán     (24jonass@gmail.com)
!   -We thank the support and advice of our molecular simulations teacher Dr. Laura Lorenia Yeomans Reyna
!*********************************************************************************************************************
!   2. DESCRIPTION:
!       THIS SUBROUTINE ALLOW US TO START OUR SYSTEM OF PARTICLES (PEOPLE) WITH A RANDOM INITIAL CONFIGURATION
!       THIS MEANS THAT THEY WILL BE PLACED INSIDE THE SIMULATION BOX IN A RANDOM WAY
!*********************************************************************************************************************

!***********************
!VARIABLE DECLARATION
!***********************
IMPLICIT NONE
REAL, DIMENSION(N):: X, Y !Arrays for position coordinates
INTEGER, INTENT(IN):: N !Number of particles
REAL, INTENT(IN):: BOXL !Box size
INTEGER :: i, j !Iterators
REAL :: RANX, RANY
REAL:: R, S
REAL :: XIJ, YIJ !Difference between a pair of particles' axis values
REAL :: R0 !Distance between 2 particles


!Creating a new file where we'll store the initial configuration
OPEN (10, FILE="CIniA.dat", STATUS='unknown')

!Iterating over each particle
DO i = 1, N

    !Here we're reinitializing the seed of the random number generator
    CALL initRandomSeed()

    !Generating uniform random numbers
69   CALL RANDOM_NUMBER(RANX)
     CALL RANDOM_NUMBER(RANY)

    !The following is for placing the particles inside the box, leaving a blank frame of width 0.5, also called mat.
     R = RANX - 0.5D0
     S = RANY - 0.5D0
     X(i) = R*(BOXL-1)
     Y(i) = S*(BOXL-1)


    DO j = 1, i-1
        XIJ = X(i) - X(j)
        YIJ = Y(i) - Y(j)
        R0 = SQRT(XIJ**2 + YIJ**2)

        !If distance between pair of particles is less than the diameter=1 of a particle, then they overlap and then new random numbers should be computed
        IF (R0 <= 1.0) THEN
            WRITE(*,*) "Overlap", i, j
            GO TO 69
        END IF
    END DO

    !Here we are writing in a file the initial configuration
    WRITE(10,*) SNGL(X(I)), SNGL(Y(I))

END DO

END SUBROUTINE
