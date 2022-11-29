function [nwz] = get_w(n,dz,kbcz)
%%      can develop to nk*nz matrix, so that the speed can be considerably improved
nwz      = 0.*kbcz;
nk       = length(kbcz(1,:));
% nwz      = zeros(1,n+1);
nwz(1,:) = 1;

%%  
%----------------------------------------%
%---- obtain the coefficient matrix  ----%
%----------------------------------------%
cfM_11   = dz(2,:).^2.*kbcz(2,:)+2;  %w2
cfM_12   = -1;

%%
rht      = zeros(n-1,nk);  
y        = rht;
gm       = y;

rht(1,:) = nwz(1,:);

%----------------------------------------%
%--- start the iteration for the matrix--%
%----------------------------------------%
%% b_i diagonal, a_i, lower half, c_i upper half
y(1,:)  = rht(1,:)./cfM_11;
gm(1,:) = cfM_12./cfM_11;
for  i  = 2:(n-2)
    %% from governing equation
        cfM_ii   = dz(i+1,:).^2.*kbcz(i+1,:)+2;  
        beta_i   = cfM_ii+gm(i-1,:);
        gm(i,:)  = -1./beta_i;
        y(i,:)   = (rht(i,:)+y(i-1,:))./beta_i;        
end

%% w_n
cfM_ii      = dz(n,:).^2.*kbcz(n,:)+2;  
i           = n-1;
beta_i      = cfM_ii+gm(i-1,:); 
y(i,:)      = (rht(i,:)+y(i-1,:))./beta_i;   
nwz(n,:)    = y(n-1,:);
for i = (n-2):-1:1
    nwz(i+1,:) = y(i,:)-gm(i,:).*nwz(i+2,:);
end

end