function [ output_args ] = test()
omegav         = 1:0.05:7;
h_1 = 3;
h_2 = 1;
g = 9.8;
theta = 0;
Fr_0 = 0;
for i=1:length(omegav)
omega_fnvs_1   = @(k) sqrt(tanh(k.*h_1)*g.*k)+k.*Fr_0.*sqrt(g)*cos(theta)-omegav(i);
omega_fnvs_2   = @(k) sqrt(tanh(k.*h_2)*g.*k)+k.*Fr_0.*sqrt(g)*cos(theta)-omegav(i);

k_1           =omegav(i)^2/g; 

k_1nvs         = fzero(omega_fnvs_1,k_1);
k_2nvs         = fzero(omega_fnvs_2,k_1);

k_2 = k_1;
cg_nvs_1       = sqrt(tanh(k_1.*h_1).*g./(k_1))/2.*(1+2*(k_1.*h_1)./sinh(2*(k_1.*h_1)))+Fr_0.*sqrt(g)*cos(theta);
cg_nvs_2       = sqrt(tanh(k_2.*h_2).*g./(k_2))/2.*(1+2*(k_2.*h_2)./sinh(2*(k_2.*h_2)))+Fr_0.*sqrt(g)*cos(theta);

AR(i) = cg_nvs_1/cg_nvs_2;
end
plot(omegav,AR)
end

