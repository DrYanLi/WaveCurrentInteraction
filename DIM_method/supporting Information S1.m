% readme for explanations of different files; Any questions can be directed
% to Dr Li Yan, lyhust306@gmail.com. 

%% EXPLANATIONS %%
% here follow explanationas of the functions in the folder named
% 'the_DIM_text_examples'. Examples for a sloping seabed can be found in
% the folder named 'slopingbed'

%% get_w.m %%
% get_w.m returns the normalized vertical velocity
% w(k_x,k_y,z)/w(k_x,k_y,0) with a given guess of the phase velocity c
% defined relative to the surface velocity of a shear current, if any.

%% information w.r.t a shear current
% shear_data.m and shearprofile.m return the information of an example shear current
% where U(z) is a 6th-order polynomial. More general profiles are
% possible. As for other general profiles, one needs to change (or combine) the two
% functions accordingly and the obtain_data.m function; 

%% KC and EL first-order approximation: obtain_kc_EL_simpsons.m %%
% obtain_kc_EL_simpsons.m returns the 1st-order approximate phase velocity from Kirby
% and Chen (1989) and Ellingsen and Li (2017) by using Simpson's method
% for the evaluation of the integral involved

%% the DIM method %%
% 'solve_omega_iter.m' and 'solve_omega_iterE.m' return numerical results based
% on the DIM method, the former of which limites the total number of
% iterations for convergence whereas the latter presets a tolerance for
% convergence. The initial guess of phase velocity required for the DIM
% method can be determined by selecting different choices made possible in
% the function obtain_data.m. The choices includes c_0 =
% \sqrt(g\tanh(kh)/k), the 1st-order approxiamtions from KC or EL. 

%% test cases %%
% the function named 'script_plot_final_...' are test cases wherein 4 or
% at least 3 methods are employed, including KC, EL, DIM, and the piecewise linear
% approxiamtion (PLA). Data from the PLA implementation of Smeltzer & Ellingsen (2017) 
% is included in the .dat files.  



