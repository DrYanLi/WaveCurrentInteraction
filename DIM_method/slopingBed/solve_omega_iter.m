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
th      = tanh(k*h);
thdk    = th./k;
kdu0    = kx*I.dU0(1)+ky*I.dU0(2);
ku0_z   = (kx_z*I.U0(1)+ky_z*I.U0(2));
kuz_z   = (kx_z.*I.Ukx_i+ky_z.*I.Uky_i);
kDtu0_z =  kuz_z-ku0_z;
kdduz_z = (kx_z.*I.ddUkx_i+ky_z.*I.ddUky_i); 
% kduz    = (kx_z.*I.dUkx_i+ky_z.*I.dUky_i);
kz_i    = k_z.*z_k;
shdch   = (exp(kz_i)-exp(-kz_i-2*k_z*h))./(1+exp(-2*k_z*h));

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

%% start a loop to calculate the exact solution
while max_convg>eps 
    j            = j-1;
    epw          = 0.5; 
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
    converg      = (delta_c./tc./d_delta_c);   
    max_convg    = max(abs(converg));
    R_err0       = converg;
    if j<1
        fprintf('M<0')
        break
    end
end

%% NOW calculate c_g FROM the DIM paper
% tsgm     = tc.*k; 
% tsgmsq   = tsgm.^2; 
% sigm0sq  = g.*thdk.*k.^2;
% sigm_z   = -ku_minus_kcz; 
% wzsq     = wz.^2;
% 
% S_k      = kdu0./k./tsgm;
% S_x      = I.dU0(1)./2./tsgm;
% S_y      = I.dU0(2)./2./tsgm;
% 
% Qz_x     = I.dUkx_i./sigm_z;
% Qz_y     = I.dUky_i./sigm_z;
% Qz_k     = kdduz_z./sigm_z./k_z/2; 
% 
% sigdsig  = sigm0sq./tsgmsq./th; 
% 
% itg_ivs  = k_z.*tc./sigm_z.*Qz_k.*wzsq; 
% itg_Nx   = (Qz_k.*k_z.*I.Ukx_i./sigm_z+Qz_x-kx_z).*wzsq;
% itg_Ny   = (Qz_k.*k_z.*I.Uky_i./sigm_z+Qz_y-ky_z).*wzsq;
% 
% itg_ivs(isnan(itg_ivs)) = 0; 
% itg_Nx(isnan(itg_Nx))   = 0;
% itg_Ny(isnan(itg_Ny))   = 0; 
% 
% Igr_ivs  = simpsons( itg_ivs, dz );
% Igr_Nx   = simpsons( itg_Nx,  dz );
% Igr_Ny   = simpsons( itg_Ny,  dz );
% 
% bI_vs    = Igr_ivs + sigdsig-S_k; 
% N_x      = Igr_Nx+ (sigdsig-S_k).*I.U0(1)./tc-S_x+sigdsig.*kx./k; 
% N_y      = Igr_Ny+ (sigdsig-S_k).*I.U0(2)./tc-S_y+sigdsig.*ky./k; 
% 
% cg_Q     = sqrt(N_x.^2+N_y.^2)./bI_vs.*tc; 
end






