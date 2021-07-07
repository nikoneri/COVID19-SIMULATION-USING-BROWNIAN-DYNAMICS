SUBROUTINE CONFIGINIREG(X,Y,N,BoxL)
!*********************************************************************************************************************
!   1. AUTHORS AND E-MAIL CONTACT:
!   -Isaac Neri Gomez Sarmiento (isaacneri.gs@gmail.com)
!   -César Omar Ramírez Álvarez (cesaromarramirezalvarez@gmail.com)
!   -Jonas Valenzuela Terán     (24jonass@gmail.com)
!   -We thank the support and advice of our molecular simulations teacher Dr. Laura Lorenia Yeomans Reyna
!*********************************************************************************************************************
!   2. DESCRIPTION:
!       THIS SUBROUTINE ALLOW US TO START OUR SYSTEM OF PARTICLES (PEOPLE) WITH A REGULAR INITIAL CONFIGURATION
!       THIS MEANS THAT THEY WILL BE PLACED INSIDE THE SIMULATION BOX IN AN ORDERED WAY
!*********************************************************************************************************************

!***********************
!VARIABLE DECLARATION
!***********************
IMPLICIT NONE
INTEGER, INTENT(IN) :: N !Number of particles
REAL, INTENT(IN) :: BoxL !Box size
INTEGER :: i, j, l !Counters
INTEGER :: m !m=int(sqrt(N)). This is the length for a side of the square simulation box
REAL:: Sigma !Diameter of a particle. Default is 1
REAL :: XL, YL !Position coordinates
REAL :: Dsep !Distance between particles
REAL, DIMENSION(N) :: X, Y !Arrays that holds position coordinates for every particle

!Creating a new file to save initial position coordinates for every particle in the system.
OPEN (40, FILE='CIniR.dat', STATUS='unknown')

!Initializing position coordinates arrays to zero.
X = 0
Y = 0

Sigma = 1
l = 0
m = INT(N**(1./2.))
Dsep = (BoxL - Sigma)/(m - 1)

!Placing particles in the simulation box
DO i = 1, m
    XL = (-BoxL/2 + (Sigma/2.)) + (i-1)*Dsep
    DO j = 1, m
         YL = (-BoxL/2 + (Sigma/2.)) + (j-1)*Dsep
             l=l+1
             X(l)=XL
             Y(l)=YL
            !Writing initial configuration in a file
            WRITE(40,*) SNGL(X(l)), SNGL(Y(l))
    END DO
END DO

END SUBROUTINE
