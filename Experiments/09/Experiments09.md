# Experiments 09 

- This experiment was initialized with the sea ice state at 1993-12-31:00:00:00. This sea ice state was obtained by running the 
Simulation at 40km resolution with the controlled parameters for the VP rheology using the ideal configuration. This experiment number is #11. 

- It was runned for a month at the 40 km resolution. 

- I set $\zeta$ = 0 

- I used the following parameters

    | Parameter    | Value |
    | :--------: | :-------: |
    | $\mu_0$  | 1.3d-01     |
    | $\mu_\infty$ | 4.0d-01     |
    | $\mu_b$ | 7.3d-01 |
    | $I_0$    | 6.8d-03    |
    | $\overline d$    | 40d03*7    |
    | $c_\Phi$    | 1    |
    | $\Phi_0$    | 1   |


## Questions
- I'm capping zeta at 0 if the shear is 0, should I do this? 

- Also, as specified like 
$$
    p_{eq} = \rho \left(\frac{\overline{d} \dot{\epsilon_{II}}}{\Phi(I) - \Phi_0}\right)^2
$$
The pressure has units of $[p_{eq}] = \text{Nm}^{-2}$ which makes it non comparable to the pressure specified by Hibler which $[P_{max}] = \text{Nm}^{-1}$. So, I multiplied $I$ and $p_{eq}$ as followed to obtain the right units. 


- There are in fact quite a lot of singularities possible in this model.
  
1. In $p_{eq}$: 
$$
    p_{eq} = \rho h \left(\frac{\overline{d} \dot{\epsilon_{II}}}{\Phi(I) - \Phi_0}\right)^2
$$
To what value should I set the pressure if $\Phi(I) - \Phi_0 \rightarrow 0$?

2. In $I$: 
$$
    I = \overline{d} \dot{\epsilon}_{II} \sqrt{\frac{h \rho}{p_{eq}}}
$$
To what value I set I if $p_{eq} \rightarrow 0$?

3. In $\zeta$ and $\eta$: 
$$
    \zeta = \frac{\mu(I)}{\dot{\epsilon_{II}}}p_{eq} \\
    \eta = \frac{\mu_b}{\dot{\epsilon_{II}}}p_{eq}
$$
To what value I set $\zeta$ and $\eta$ as $\dot{\epsilon_{II}} \rightarrow 0$? 

