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

#### <b> GIF for the first scenario</b>

<p align="center">
<img src="https://github.com/nikoneri/COVID19-SIMULATION-USING-BROWNIAN-DYNAMICS/blob/main/4.%20GRAPHS%20AND%20GIF%20ANIMATION/CASO1_aleatorio10%25.gif" width="450px" >
</p>

### 2.1 HIGHLIGHTS
<ul>
    <li>It was possible to find analogies between the physical model and the model of infection, enabling us to explore the effectivity of sanitary measurements using brownian dynamics.</li>
    
<li>Based in our results, following both sanitary measures (Keep your Distance and Stay Home) reduce considerably the number of new infected people.</li>
<li>By itself, the sanitary measure Keep your Distance was more effective than Stay Home, according to the chosen parameters. In the ideal case where Keep your Distance is being followed by all people at any moment and place, it might considerably reduce the velocity of infection, because people would be out of the distance of infection. </li>
    <li>This work serves as a teaching tool for learning Brownian Dynamics and make people aware of the importance of following the sanitary measurements to face the COVID-19 pandemic.</li>
</ul>
## 3. SOME REMARKS:
<ul>
<li>This work is originally in Spanish, so some of the variables and subroutines are in that language.</li>
<li>We use the Ermak algorithm to move particles (people) and solve overdamped Langevin equations or in the diffusive regime. </li>
<li>This simple contagion model doesn't include recovered people or inmmune people. Including this in the future will improve the model.</li>
<li>We consider people to be 2D colloids with a diameter of 1 m.</li>
 <li>The infection distance is by deafult DINF=2.5.</li>
<li>We work with reduced units, also called dimensionless units.</li>
<li>We work with a repulsive potential of the type U=(r0/r)^v.</li>
<li>Besides computing the number of infected people, this program can also compute structural and dynamical properties:
    <ol>
    <li> Radial distribution function g(r), which gives us an insight on the most likely average distance between particles.</li>
    <li> Mean Square Displacement W(t), which gives us an insight on how much area the particles spread out.</li>
    <li> Diffussion Coefficient D(t), which gives us an insight on the mobility (speed) of particles.</li>
    </ol>
    </li>
<li> This work has been presented at two congress. In this repo we attached 2 PDF for our presentation and poster, where  you'll find a detailed description of our modelling and more results for each of the scenarios (in Spanish).
    <ol>
    <li>Presentation: LXIII National Congress of Physics in Mexico (you can watch our presentation at: https://youtu.be/47bXLXEzZDo?t=8183). </li>
    <li>Poster: XXXI National Week of Research and Teaching in Mathematics in Hermosillo, Sonora, Mexico. </li>
    </ol>
    </li>
</ul>
