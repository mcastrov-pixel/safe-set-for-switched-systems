# safe-set-for-switched-systems
Generate Maximal Output Admissible Set (MOAS) for switched systems (with controlled switching signal)

This code (and examples) are taken from our work [1]. Please make sure to properly cite us if you found this useful :)



Two examples are available:
  main files are  **dynamics1_main.m** and **dynamics2_main.m**

**(Figures 1 to 4)** The figures below illustrate how switches between different modes can be beneficial to increase the set of safe states for a switched linear system with two modes (takend from dynamics2_main)
**(Figure 5)** With the proposed approach we can also construct safe trajectories that prioritize different objectives: maybe we wish to always stay in a specfic mode or maybe we want to get to the origin faster


<figure align="center">
  <img src="https://github.com/user-attachments/files/29142780/MOAS_setS_dyn2_1.pdf" width="75%">
  <figcaption><b>Figure 1:</b> MOAS obtained when considering a switching signal that has no changes (either only mode 1 or only mode 2) .</figcaption>
</figure>


<figure align="center">
  <img src="https://github.com/user-attachments/files/29142781/MOAS_setS_dyn2_2.pdf" width="75%">
  <figcaption><b>Figure 2:</b> MOAS obtained  when considering a switching signal that has two changes before becoming constant (1,2,1,1,1... or 2,1,2,2,2...).</figcaption>
</figure>



  
<figure align="center">
  <img src="https://github.com/user-attachments/files/29142782/MOAS_setS_dyn2_3.pdf" width="75%">
  <figcaption><b>Figure 3:</b> MOAS obtained  when considering a switching signal that has four changes before becoming constant (similar to previous case).</figcaption>
</figure>



<figure align="center">
  <img src="https://github.com/user-attachments/files/29142373/MOAS_setS_dyn2_4.pdf" width="75%">
  <figcaption><b>Figure 4:</b> MOAS obtained  when considering a switching signal that has six changes before becoming constant(similar to previous case).</figcaption>
</figure>




<figure align="center">
  <img src="https://github.com/user-attachments/files/29142600/trajectories.pdf" width="75%">
  <figcaption><b>Figure 5:</b> Both left and right show trajectories for the same dynamics but with different objectives guiding the safe switching signal generation. On the left we try to stay in mode 1 on the right we try to go to the origin ASAP).</figcaption>
</figure>



**uses the following toolboxes**

Optimization problem are solved with quadprog using Yalmip [https://yalmip.github.io/download/]
uses MPT3 for polytope ploting: https://www.mpt3.org/pmwiki.php/Main/Installation
Cora toolbox for some light operations on polytopes: https://tumcps.github.io/CORA/
