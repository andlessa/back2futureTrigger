### Notes on Analysis Strategy


1. Maybe the readout time for particles produced within the tracker is not always right? 
   If the particle is produced from a LLP decay within the tracker, the time for the LLP to travel
   from the interaction point to its decay is not $`\Delta t = t_0 - |\vec{r}_0|/c`$. since the LLP might 
   not travel at the speed of light.
> Fair point, something which I have not thought about in detail -- Also, when a particle is produced (by the LLP decay) in the ECAL or HCAL, we're currently assuming that the LLP has travelled with speed of light from the interaction point to its decay vertex. I think we could do better, and check in each extrapolation if a $`\chi_{1}`$ was involved in the previous production chain of this particle, and if so, we could account for the LLP velocity.

3. How do you define the momenta for each particle in a given individual TT cell in Eq.5? Is it something like: $`(p_x,p_y) = \frac{E}{\cosh\eta_{J}}(\cos\phi_{J},\sin\phi_{J})`$, where $`J`$ enumerates the TT cells? If so, then the TT cell momentum is simply   $`p^{J}_{x,y} = \frac{E_J}{\cosh\eta_{J}}(\cos\phi_{J},\sin\phi_{J})`$, where $`E_J`$ is the total cell deposit. Is this correct?
> Actually, we currently just take the 4 momentum of the individual particles to get $`(p_x,p_y)`$. But some time ago, I compared this approach to the method you mentioned above, and if I remember correctly there was not much difference. I'll check this again, and perhaps we should tend to go for the "energy" approach rather than the "4 momentum", because the first one is closer to what is done on trigger level.  

4. The EMF and $`\Delta R(jet,tracks)`$ requirements for the HLT seems to require that the LLP decay in the $`N`$ event record takes places within the HCAL. Otherwise the LLP decay would generate visible particles assigned to the ECAL TTs and we would have $`E_{EM} > E_{HAD}`$. Is this correct?
> When the LLP decay is taking place in the HCAL, the EMF requirement is trivially fulfilled;
> Also, for LLP decays in the ECAL, the EMF condition might be fulfilled, in case the decay position is close to the outer boundary of the ECAL -- we currently require (guesstimate) that the remainig ECAL distance must be at most one hadronic interaction length.
> LLP decays in the trackwer would not fulfill the EMF condition.
> The $`\Delta R(jet,tracks)`$ condition additionally requires that there are no tracks in the tracker pointing to the jet.

5. Why was the $`\Delta \phi`$ cut modified?
> You're probably referring to $`\Delta \phi (\text{MET,jet}) < 1.0`$ v.s. $`\Delta \phi (\text{jet,jet}) > \pi - 1.0`$. The trigger is actually using the former version of this cut, i.e. it is searching for MET in N-1 followed by a jet in N relatively close in phi; in contrast, we previously approximated this cut by $`\mathrm{MET} = -E_T`$, and therefore selecting events with $`\Delta \phi (\text{jet}_N , \text{jet}_{N-1}) > \pi - 1`$. In the updated analysis strategy, we're treating the MET in N-1 properly, and can therefore directly apply the $`\Delta \phi`$ cut as it is implemented in the Trigger.
