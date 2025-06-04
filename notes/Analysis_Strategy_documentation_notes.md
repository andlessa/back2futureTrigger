### Notes on Analysis Strategy


1. Maybe the readout time for particles produced within the tracker is not always right? 
   If the particle is produced from a LLP decay within the tracker, the time for the LLP to travel
   from the interaction point to its decay is not $`\Delta t = t_0 - |\vec{r}_0|/c`$. since the LLP might 
   not travel at the speed of light.

2. How do you define the momenta for each particle in a given individual TT cell in Eq.5? Is it something like: $`(p_x,p_y) = \frac{E}{\cosh\eta_{J}}(\cos\phi_{J},\sin\phi_{J})`$, where $`J`$ enumerates the TT cells? If so, then the TT cell momentum is simply   $`p^{J}_{x,y} = \frac{E_J}{\cosh\eta_{J}}(\cos\phi_{J},\sin\phi_{J})`$, where $`E_J`$ is the total cell deposit. Is this correct?

3. The EMF and $`\Delta R(jet,tracks)`$ requirements for the HLT seems to require that the LLP decay in the $`N`$ event record takes places within the HCAL. Otherwise the LLP decay would generate visible particles assigned to the ECAL TTs and we would have $`E_{EM} > E_{HAD}`$. Is this correct?

4. Why was the $`\Delta \phi`$ cut modified?

