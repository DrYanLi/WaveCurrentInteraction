function [tc,wz,R_err0] = solve_omega_iter(I,Ik,tc0,z,eps,M)
%% size(k)*N
g       = I.g;
h       = I.h;
k_z     = Ik.k_i;
k       = Ik.k;
kx      = Ik.kx;
ky      = Ik.ky;
kx_z    = Ik.kx_i;
ky_z    = Ik.ky_i;
z_k     = z.z_k;
dz      = z.dz;      % now a vector of the same size of k
%%
% - % -% - % - % - % - %
%     inputs above     %
% - % -% - % - % - % - % 
%%
thdk    = tanh(k*h)./k;
kdu0    = kx*I.dU0(1)+ky*I.dU0(2);
ku0_z   = (kx_z*I.U0(1)+ky_z*I.U0(2));
kuz_z   = (kx_z.*I.Ukx_i+ky_z.*I.Uky_i);
kDtu0_z =  kuz_z-ku0_z;
kdduz_z = (kx_z.*I.ddUkx_i+ky_z.*I.ddUky_i); 
% kduz    = (kx_z.*I.dUkx_i+ky_z.*I.dUky_i);
kz_i    = k_z.*z_k;
shdch   = (exp(kz_i)-exp(-kz_i-2*k_z*h))./(1+exp(-2*k_z*h));

%% obtian the initial guess
% [tc_app]     = c_app(g,h,k,kduz,dm_nz,dz); %% To approximate an initial c_0
% [tc_app_z,~] = meshgrid(tc_app,dm_nz);

%% output the integral as a function of k and c
tc           = tc0.tc;
tc_z         = tc0.tc_i;
ku_minus_kcz = kDtu0_z-k_z.*tc_z;
kddu_dv_kc   = kdduz_z./ku_minus_kcz;
kddu_dv_kc(isnan(kddu_dv_kc))=0;  % avoid a critical layer

%% estimate the initial error
kbcz      = kx_z.^2+ky_z.^2+ kddu_dv_kc;
[wz]      = get_w(I.N,dz,kbcz);
ig_itg    = kdduz_z.*wz.*shdch./k_z./ku_minus_kcz;
dig_itg   = kdduz_z.*wz.*shdch./ku_minus_kcz.^2;
ig_itg(isnan(ig_itg))  = 0;
dig_itg(isnan(dig_itg))= 0;

d_ig      = simpsons(dig_itg,dz);
ig        = simpsons(ig_itg,dz);
delta_c   = (1+ig).*tc.^2+tc.*thdk.*kdu0./k-g.*thdk;
d_delta_c = 2*(1+ig).*tc+thdk.*kdu0./k+d_ig.*tc.^2;
converg   = abs(delta_c./tc./d_delta_c);
R_err0    = converg;
max_convg = max(converg);
j         = M;  %iterative number

%% start a loop to calculate an exact solution
while max_convg>eps 
    if j<1
        break
    end
    j            = j-1;
    epw          = 1; 
    tc           = tc-epw.*delta_c./d_delta_c;  
%     [tc_z,~]     = meshgrid(tc,dm_nz);
%     ku_minus_kcz = kDtu0_z-k_z.*tc_z;

    ku_minus_kcz = kDtu0_z-k_z.*tc;
    kddu_dv_kc   = kdduz_z./ku_minus_kcz;
    
    kddu_dv_kc(isnan(kddu_dv_kc))=0;  % avoid a critical layer   
    
    kbcz         = kx_z.^2+ky_z.^2+ kddu_dv_kc;
    [wz]         = get_w(I.N,dz,kbcz);
    ig_itg       = kdduz_z.*wz.*shdch./k_z./ku_minus_kcz;
    dig_itg      = kdduz_z.*wz.*shdch./ku_minus_kcz.^2;
    
    ig_itg( isnan(ig_itg) )  = 0;
    dig_itg( isnan(dig_itg) )= 0;
    
    d_ig         = simpsons( dig_itg, dz );
    ig           = simpsons( ig_itg, dz );
    delta_c      = (1+ig).*tc.^2+tc.*thdk.*kdu0./k-g.*thdk;
    d_delta_c    = 2*(1+ig).*tc+thdk.*kdu0./k+d_ig.*tc.^2;
    converg      = abs(delta_c./tc./d_delta_c);   
    max_convg    = max(abs(converg));
    R_err0       = converg;
end
end






