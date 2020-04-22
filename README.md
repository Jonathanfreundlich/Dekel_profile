<h2> Analytic expressions for the Dekel+ profile </h2>

<p align="justify">
This python program enables to compute the analytical expressions of the Dekel+ dark matter density profile. It is based on the following article: 
</p>

<p align="justify">
<a href="https://ui.adsabs.harvard.edu/abs/2020arXiv200408395F"  style="text-decoration:none" class="type1"><b>The Dekel+ profile: a mass-dependent dark-matter density profile with flexible inner slope and analytic potential, velocity dispersion, and lensing properties</b></a> 
<a href="https://ui.adsabs.harvard.edu/link_gateway/2020arXiv200408395F/EPRINT_PDF" style="text-decoration:none" class="type1"> [PDF] </a>
</p>

<p align="justify">
Quantities expressed analytically in this program: 
<ul>
       <li>The density </li>
       <li>The average density</li>
       <li>The enclosed mass</li>
       <li>The orbital velocity</li>
       <li>The logarithmic density slope</li>
       <li>The gravitational potential</li>
       <li>The velocity dispersion</li>
</ul>       
</p>

<p align="justify">
For the velocity dispersion, it also provides two expressions in terms of finite sums (<a href="https://ui.adsabs.harvard.edu/abs/1996MNRAS.278..488Z/abstract"  style="text-decoration:none" class="type1">Zhao 1996,</a> Eq. 19, and <a href="https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.4523F/abstract"  style="text-decoration:none" class="type1">Freundlich et al. 2020a</a>, Eqs. B8 and B10). Some of the quantities are also defined for other profiles (NFW, Einasto, Di Cintio+, generalised NFW with flexible inner slope). 
</p>

<p align="justify">
Quantities integrated numerically in this program: 
<ul>
       <li>The distribution function</li>
       <li>The projected surface density</li>
       <li>The projected mass</li>
       <li>The convergence</li>
       <li>The shear</li>
</ul>       
</p>

<p align="justify">
       The program also implements the <b>mass-dependent prescriptions</b> for the Dekel+ and the Di Cintio+ profiles. 
</p>

<h3 align="justify">
       Example:  
</h3>

<pre><code>
import numpy as np
import matplotlib.pyplot as plt
import Dekel_profile as prf
</code></pre>

