function [bp] = get_rhs(n,n_inlet,dx,dy,rho,u_star,v_star)

bp=zeros((n-1)*(n-1),1);    %vector of RHS for solving pressure corrections
stride = n-1;

% RHS is the same for all nodes except the p_prime(1,1)
% because p(1,1) is set to be zero, it has no pressure correction
for j=2:n
    for i=2:n
        position = (i-1) + (j-2)*stride; 
        bp(position) = rho * (u_star(i-1,j)*dy - u_star(i,j)*dy + v_star(i,j-1)*dx - v_star(i,j)*dx); 
        
    end
end

% modify for p_prime outlet
bp(((n-1)*(n-1))-(n_inlet-2):((n-1)*(n-1)),1) = 0;

return 
end
