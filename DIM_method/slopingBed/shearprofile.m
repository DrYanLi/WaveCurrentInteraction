function [I,z] = shearprofile(cst_str,dm_nz)
%% define the shear profile 
g     = cst_str.g;
h     = cst_str.h;
Fr    = cst_str.Fr;
alpha = cst_str.alpha; 

% h_0   = ; %
h_2   = 1;
%%
        uxf      = @(z) Fr*sqrt(g*h_2).*exp(alpha*z/h_2)-Fr*sqrt(g*h_2);        
        duxf     = @(z) Fr*sqrt(g/h_2)*alpha.*exp(alpha*z/h_2);       
        dduxf    = @(z) Fr*sqrt(g/h_2^3)*alpha^2.*exp(alpha*z/h_2); 
         
%         uxf      = @(z) Fr*sqrt(g*h).*exp(alpha*z/h)+sqrt(9.8/h)*z;        
%         duxf     = @(z) Fr*sqrt(g/h)*alpha.*exp(alpha*z/h)+sqrt(9.8/h);       
%         dduxf    = @(z) Fr*sqrt(g/h^3)*alpha^2.*exp(alpha*z/h);
              
        uyf        = @(z) 0.*z;  
        duyf       = @(z) 0.*z;  
        dduyf      = @(z) 0.*z;  
        
%% 
        I.U_i        = [];
        I.dU_i       = [];
        I.ddU_i      = [];    
        I.U0(1)      = uxf(0);
        I.U0(2)      = uyf(0);
        I.dU0(1)     = duxf(0);
        I.dU0(2)     = duyf(0);        
        
        %% place holders
        I.g       = 9.8;
        I.h       = [];
        I.N       = []; % intervals, N+1 nodes
        z.n_dz    = []; %1*nz 
        z.dz      = []; % h/N
       
        
        %% values
        z.dm_nz  = []; % vector 1*nz;
        z.z_k    = dm_nz;
        
        %% size = nz*nk
        [I.Ukx_i]      = uxf(dm_nz);
        [I.Uky_i]      = uyf(dm_nz);
        [I.dUkx_i]     = duxf(dm_nz);
        [I.dUky_i]     = duyf(dm_nz);
        [I.ddUkx_i]    = dduxf(dm_nz);
        [I.ddUky_i]    = dduyf(dm_nz);
            
end

