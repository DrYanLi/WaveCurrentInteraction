function [I,z] = shearprofile(pxinput,pyinput,dm_nz)
%% define the shear profile   
   %  pxinput is the coefficients vector of Ux(z)
   %  pyinput is the coefficients vector of Uy(z)
        p1       = pxinput(1);
        p2       = pxinput(2);
        p3       = pxinput(3);
        p4       = pxinput(4);
        p5       = pxinput(5);
        p6       = pxinput(6);
%         p7       = pxinput(7);
        p1y      = pyinput(1);
        p2y      = pyinput(2);
        p3y      = pyinput(3);
        p4y      = pyinput(4);
        p5y      = pyinput(5);
        p6y      = pyinput(6);
%         p7y      = pyinput(7);
        
        uxf      = @(z) polyval(pxinput,z);
        pcofdx   = [6*p1 5*p2 4*p3 3*p4 2*p5 p6];
        duxf     = @(z) polyval(pcofdx,z);       
        pcofddx  = [30*p1 20*p2 12*p3 6*p4 2*p5];
        dduxf    = @(z) polyval(pcofddx,z);
        
        uyf      = @(z) polyval(pyinput,z);       
        pcofdy   = [6*p1y 5*p2y 4*p3y 3*p4y 2*p5y p6y];
        duyf     = @(z) polyval(pcofdy,z);        
        pcofddy  = [30*p1y 20*p2y 12*p3y 6*p4y 2*p5y];
        dduyf    = @(z) polyval(pcofddy,z);
        
        %% VALUES
%         I.U_i(1,:)   = uxf(dm_nz);
%         I.U_i(2,:)   = uyf(dm_nz);
%         I.dU_i(1,:)  = duxf(dm_nz);
%         I.dU_i(2,:)  = duyf(dm_nz);
%         I.ddU_i(1,:) = dduxf(dm_nz);
%         I.ddU_i(2,:) = dduyf(dm_nz);        
        I.U_i        = [];
        I.dU_i       = [];
        I.ddU_i      = [];    
        I.U0(1)      = uxf(0);
        I.U0(2)      = uyf(0);
        I.dU0(1)     = duxf(0);
        I.dU0(2)     = duyf(0);        
        
        %% place holders
        I.g       = 9.8;
        I.h       = 1;
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

