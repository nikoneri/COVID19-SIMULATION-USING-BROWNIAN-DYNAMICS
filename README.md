# COVID19-SIMULATION-USING-BROWNIAN-DYNAMICS

   ## 1. AUTHORS AND E-MAIL CONTACT:
   <ul>
   <li>Isaac Neri Gomez Sarmiento (isaacneri.gs@gmail.com) </li>
   <li> César Omar Ramírez Álvarez (cesaromarramirezalvarez@gmail.com) </li>
   <li> Jonas Valenzuela Terán     (24jonass@gmail.com) </li>
   <li>We thank the support and advice of our molecular simulations teacher Dr. Laura Lorenia Yeomans Reyna </li>
   </ul>
   ## 2. DESCRIPTION:
   *In this program we obtained a simple contagion model, when making certain analogies between a human population and
   a two-dimensional coloidal suspension system under brownian motion.
   We consider Hermosillo's population density (a city located in Sonora, Mexico).
   In this model we include sanitary measurements such as "Quedate en Casa" or Stay Home
   and "Sana Distancia" or Keep your Distance.

   *We can compare different scenarios such as :
   1. Not implementing any sanitary measurements (the parameters that allow this scenario are set by default in this code).
   2. Implementing just Stay Home (with a gradual intensity).
   3. Implementing just Keep your Distance (with a gradual intensity).
   4. Implementing both Stay Home and Keep your Distance (with a gradual intensity).


   ## 3. SOME REMARKS:
   *We use the Ermak algorithm to move particles (people) and solve overdamped Langevin equations or in the diffusive regime.
   *This work is originally in Spanish, so some of the variables and subroutines are in that language.
   *This simple contagion model doesn't include recovered people or inmmune people. Including this in the future will improve the model.
   *We consider people to be 2D colloids with a diameter of 1 m.
   *We work with reduced units, also called dimensionless units.
   *We work with a repulsive potential of the type U=(r0/r)^v.
   *Besides computing the number of infected people, this program can also compute structural and dynamical properties:
       1. Radial distribution function g(r), which gives us an insight on the most likely average distance between particles.
       2. Mean Square Displacement W(t), which gives us an insight on how much area the particles spread out.
       3. Diffussion Coefficient D(t), which gives us an insight on the mobility (speed) of particles.
   *This work has been presented at two congress:
       1. LXIII National Congress of Physics in Mexico.
       2. XXXI National Week of Research and Teaching in Mathematics in Hermosillo, Sonora, Mexico.
