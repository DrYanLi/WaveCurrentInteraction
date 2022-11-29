function [tc_str, R_str, U_str] = obtain_data(N,sgn,theta)

M              = 1; % total number of iterations for the DIM
sgn_initial_TC = 0; % preset the initial guess requried for the DIM, 0:c_0; 1: KC, 2: EL

%% different wind-induced profiles
if sgn == 1
    c_b = load('c_b1.dat');
    kv  = load('kv1.dat');
end
if sgn==2
    c_b = load('c_b2.dat');
    kv  = load('kv2.dat');
end
if sgn==3
    c_b = load('c_b3.dat');
    kv  = load('kv3.dat');
end
        
[px,py]   = shear_data(sgn);
g         = 9.8;
h         = 1;
eps       = 10^-6;

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
[I,z]         = shearprofile(px,py,dm_nz);
%%
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
[tc_dimM0_1,~,R01]  = solve_omega_iter(I,Ik,tc0,z,eps,M);  
M                   = 2;
[tc_dimM0_2,~,R02]  = solve_omega_iter(I,Ik,tc0,z,eps,M);  

%% 
sgn_initial_TC    = 1;
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
M                     = 1;
[tc_dimM1_1,~,R11]    = solve_omega_iter(I,Ik,tc0,z,eps,M);
[~,~,R_kc]            = solve_omega_iter(I,Ik,tc0,z,eps,0);

%%
sgn_initial_TC        = 2;
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
[~,~,R_EL]            = solve_omega_iter(I,Ik,tc0,z,eps,0);
%% 
tc_PLA      = c_b -(I.U0(1)*kv*cos(theta)+I.U0(2)*kv*sin(theta))./kv;
R01_e       = (tc_dimM0_1-tc_PLA)./tc_PLA;
R02_e       = (tc_dimM0_2-tc_PLA)./tc_PLA;
R11_e       = (tc_dimM1_1-tc_PLA)./tc_PLA;

%% outputs
tc_str.PLA   = tc_PLA;
tc_str.kc    = tc_kc;
tc_str.EL    = tc_EL;
tc_str.kv    = kv;
tc_str.dim01 = tc_dimM0_1;
tc_str.dim02 = tc_dimM0_2;
tc_str.dim11 = tc_dimM1_1;

R_str.R01    = R01;
R_str.R02    = R02;
R_str.R11    = R11;
R_str.R01_e  = R01_e;
R_str.R02_e  = R02_e;
R_str.R11_e  = R11_e;
R_str.R_kc   = R_kc;
R_str.R_EL   = R_EL;

U_str.Uz     = I.Ukx_i(:,1)/sqrt(g*h);
U_str.z      = dm_nz(:,1);
end
