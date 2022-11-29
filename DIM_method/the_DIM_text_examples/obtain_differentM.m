function [R_M1,R_kc] = obtain_differentM(sgn,N_dim)
% sgn      = 2;
% N_dim    = 6;
%%
kv1     = kv(1:97);
kv2     = kv(98:end);
c_b1    = c_b(1:97);
c_b2    = c_b(98:end);
[R_kc_1,R_M1_1,R_M2_1] = obtain_longwaves(N_dim,kv1,c_b1,sgn);

R_kc_2  = 0.*c_b2;
R_M1_2  = R_kc_2;
R_M2_2  = R_kc_2;
for i = 1:length(kv2)
    [R_kc_2(i),R_M1_2(i),R_M2_2(i)] = obtain_shortwaves(N_dim,kv2(i),c_b2(i),sgn);
end
R_kc   = [R_kc_1 R_kc_2];
R_M1   = [R_M1_1 R_M1_2];
R_M2   = [R_M2_1 R_M2_2];

end
