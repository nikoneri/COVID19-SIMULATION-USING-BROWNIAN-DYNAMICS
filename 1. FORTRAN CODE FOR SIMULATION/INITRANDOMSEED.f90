SUBROUTINE initRandomSeed()
!*********************************************************************************************************************
!   1. AUTHORS AND E-MAIL CONTACT:
!   -Isaac Neri Gomez Sarmiento (isaacneri.gs@gmail.com)
!   -César Omar Ramírez Álvarez (cesaromarramirezalvarez@gmail.com)
!   -Jonas Valenzuela Terán     (24jonass@gmail.com)
!   -We thank the support and advice of our molecular simulations teacher Dr. Laura Lorenia Yeomans Reyna
!*********************************************************************************************************************
!   2. DESCRIPTION:
!       THIS SUBROUTINE ALLOWS US TO REINITIALIZE A RANDOM NUMBER GENERATOR SEED WITH THE HELP OF THE SYSTEM CLOCK
!       SO WE AVOID HAVING REPEATING PATTERNS OF NUMBERS WHEN MAKING NEW RANDOM NUMBERS IN A LOOP
!
!   3. SOME REMARKS:
!       This particular subroutine was obtained from:
!       https://eva.fcien.udelar.edu.uy/pluginfile.php/19200/mod_resource/content/0/randomNumbers.html
!*********************************************************************************************************************

!********************
!VARIABLE DECLARATION
!********************
IMPLICIT NONE
INTEGER :: i, n, Clock
INTEGER, DIMENSION(:), ALLOCATABLE :: Seed

CALL Random_Seed(Size = n)

ALLOCATE(Seed(n))

CALL System_Clock(Count = Clock)

Seed = Clock + 37*(/(i-1,i=1,n)/)

CALL Random_Seed(put = Seed)

DEALLOCATE(Seed)

END SUBROUTINE
