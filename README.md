# COVID19-SIMULATION-USING-BROWNIAN-DYNAMICS

## 1. AUTHORS AND E-MAIL CONTACT:
<ul>
<li>Isaac Neri Gomez Sarmiento (isaacneri.gs@gmail.com) </li>
<li> César Omar Ramírez Álvarez (cesaromarramirezalvarez@gmail.com) </li>
<li> Jonas Valenzuela Terán     (24jonass@gmail.com) </li>
<li>We thank the support and advice of our molecular simulations teacher Dr. Laura Lorenia Yeomans Reyna </li>
</ul>


## 2. DESCRIPTION:
In this project we obtained a simple contagion model, when making certain analogies between a human population and
a two-dimensional coloidal suspension system under brownian motion.
We consider Hermosillo's population density (a city located in Sonora, Mexico).
In this model we include sanitary measurements such as "Quedate en Casa" or Stay Home
and "Sana Distancia" or Keep your Distance.

<b>We can compare different scenarios such as :</b>
<ol>
<li> Not implementing any sanitary measurements (the parameters that allow this scenario are set by default in the code uploaded in this repo).</li>
<li> Implementing just Stay Home (with a gradual intensity).</li>
<li> Implementing just Keep your Distance (with a gradual intensity).</li>
<li> Implementing both Stay Home and Keep your Distance (with a gradual intensity).</li>
</ol>

## 3. SOME REMARKS:
<ul>
<li>We use the Ermak algorithm to move particles (people) and solve overdamped Langevin equations or in the diffusive regime. </li>
<li>This work is originally in Spanish, so some of the variables and subroutines are in that language.</li>
<li>This simple contagion model doesn't include recovered people or inmmune people. Including this in the future will improve the model.</li>
<li>We consider people to be 2D colloids with a diameter of 1 m.</li>
<li>We work with reduced units, also called dimensionless units.</li>
<li>We work with a repulsive potential of the type U=(r0/r)^v.</li>
<li>Besides computing the number of infected people, this program can also compute structural and dynamical properties:
    <ol>
    <li> Radial distribution function g(r), which gives us an insight on the most likely average distance between particles.</li>
    <li> Mean Square Displacement W(t), which gives us an insight on how much area the particles spread out.</li>
    <li> Diffussion Coefficient D(t), which gives us an insight on the mobility (speed) of particles.</li>
    </ol>
    </li>
<li> This work has been presented at two congress:
    <ol>
    <li>LXIII National Congress of Physics in Mexico (you can watch our presentation (in spanish) at: https://youtu.be/47bXLXEzZDo?t=8183). </li>
    <li> XXXI National Week of Research and Teaching in Mathematics in Hermosillo, Sonora, Mexico. </li>
    </ol>
    </li>
</ul>
