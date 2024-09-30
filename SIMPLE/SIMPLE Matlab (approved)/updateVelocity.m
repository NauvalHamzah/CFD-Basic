function [u,v] = updateVelocity(n,n_inlet,u_star,v_star,p_prime,d_u,d_v,v1)
u = zeros(n,n+1);
v = zeros(n+1,n);


%update interior nodes of u and v
for i=2:n-1
    for j=2:n        
        u(i,j) = u_star(i,j) + d_u(i,j)*(p_prime(i,j)-p_prime(i+1,j));        
    end
end

for i=2:n
    for j=2:n-1        
        v(i,j) = v_star(i,j) + d_v(i,j)*(p_prime(i,j)-p_prime(i,j+1));        
    end
end

%update BCs
v(1,2:n-1) = -v(2,2:n-1); %left wall
v(n+1,2:n-1) = -v(n,2:n-1); %right wall
v(1:n+1, 1) = 0; %bottom wall
v(1:n+1, n) = 0; %top wall 

%%Inlet and Outlet
v(2*(n_inlet):3*(n_inlet)-2,1)=v1;   %Inlet
v((n+1)-(n_inlet-1):n,n)=v((n+1)-(n_inlet-1):n,n-1);               %Outlet


u(1,1:n+1) = 0; %left wall
u(n,1:n+1) = 0; %right wall
u(1:n, 1)  = -u(1:n, 2); %bottom wall
u(1:n, n+1) = -u(1:n, n); %top wall 

return
end
