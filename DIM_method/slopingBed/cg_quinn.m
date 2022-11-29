function [omega,cg_Q] = cg_quinn(N,kv,h,theta,Fr,alpha)
sgn_initial_TC = 2;      
g              = 9.8;
% h              = 1;
eps            = 10^-6;

%% define short and long waves
nk            = length(kv);
Ik.k          = kv;
Ik.kx         = kv.*cos(theta);  
Ik.ky         = kv.*sin(theta);
Ik.tht        = theta;

%%
v_z           = 0:-1/N:-1;
kz_cutv(1:nk) = 3.5 + 2*log(N/7);
kh_z          = min(kv.*h,kz_cutv);  %check here if it gives a vector of size nk
[k_i,dm_nz]   = meshgrid(kv,v_z); 
[nkz_cut,~]   = meshgrid(kh_z,v_z); 
dm_nz         = dm_nz.*nkz_cut./k_i;
[Ik.k_i,~]    = meshgrid(kv,dm_nz(:,1));
[Ik.kx_i,~]   = meshgrid(Ik.kx,dm_nz(:,1)); 
Ik.ky_i       = sign(Ik.ky).*sqrt(Ik.k_i.^2-Ik.kx_i.^2); 
%%
cst_str.g     = g;
cst_str.h     = h;
cst_str.Fr    = Fr;
cst_str.alpha = alpha; 

[I,z]         = shearprofile(cst_str,dm_nz);
%%
I.h           = h;
I.N           = N;
z.dz          = nkz_cut/N./k_i;
%%
[tc_kc,tc_EL]  = obtain_kc_EL_simpsons(I,Ik,z);
%%
switch sgn_initial_TC
    case 0
        tc            = sqrt(g.*tanh(kv*h)./kv);
        tc0.tc        = tc;
        tc0.tc_i      = sqrt(g.*tanh(k_i*h)./k_i);
    case 1
        tc0.tc        = tc_kc;
        [tc0.tc_i,~]  = meshgrid(tc0.tc,dm_nz(:,1));
    case 2
        tc0.tc        = tc_EL;
        [tc0.tc_i,~]  = meshgrid(tc0.tc,dm_nz(:,1));
end

%% use c0 as initial guess
[tc_dim,~,cg_Q]  = solve_omega_cgQ(I,Ik,tc0,z,eps,50);  
omega            = kv.*tc_dim;

%% outputs
end
