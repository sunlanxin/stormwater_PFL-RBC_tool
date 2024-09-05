# stormwater_PFL-RBC_tool

The `stormwater_PFL-RBC_tool` predictive fuzzy logic and rule-based control (PFL-RBC) approach for practitioners to address the real-time continuous operation of stormwater storage systems in urban areas. This approach enables more practical incorporation with accessible predictive information (i.e., the total rainfall depth forecast) and more time-efficient training process for adapting a new environment. Our codes also provide a framework for compare the performances of PFL-RBC with other control models (e.g., the static RBC and the optimization based MPC)


## Code

The main modules were coded in MATLAB, as listed below:

1. `Qtarget_FLC.m` to train a fuzzy logic controller in the PFL-RBC.

2. `Main.m` to implement an RTC model (e.g., the PFL-RBC approach) in a single storage tank case.

3. `control_model_wet.m` to generate a control strategy during wet period.

4. `control_model_dry.m` to generate a control strategy during dry period.

5. `TVGM_URBAN.p` to simulate the total inflow from upstream.

6. `SWMM.p` to simulate the dynamics of controlled storage systems. 

Noted that:

- the fuzzy logic controller was simulated using the `fuzzy logic designer` in MATLAB.

- the stepwise process of control decision-making and execution were coded in the framework of `MatSWMM` (developed by Riano-Briceno et al. in 2016), which alows the SWMM model to be paused at each control time step to extract the sataes and adjustthe set points of the orifice and pumps.

- the total inflow can also be simulated using the `SWMM` if the required modeling data is available.


## Requried Data

-  the SWMM.inp file of the demonstration case
-  the forecasted and monitored rainfall data


## Running steps

1. download the `MATLAB` software (version not lower than R2021a is suggested). 

2. download the `stormwater_PFL-RBC_tool` and open it at the MATLAB. 

3. prepare a demonstration case for control simulation (refer to the `Requried Data`)

4. run the `Qtarget_FLC.m` to get a trained FLC controller (i.e., `target flow.fis`) .

5. change the parameters in the `Main.m` according to your case.

6. run the `Main.m` and output the simulation results.
