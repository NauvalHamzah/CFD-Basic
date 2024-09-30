function [u_star,v_star] = vel_star(n,n_inlet,rhsU,rhsV,Au,Av,v1)
u_star = zeros(n,n+1);
v_star = zeros(n+1,n);

u_star_interior = pentaDiag_solve(Au,rhsU);
v_star_interior = pentaDiag_solve(Av,rhsV);

%convert velocity correction in to a matrix
%update velocity values
z=1; 
for j=2:n
    for i=2:n-1
        u_star(i,j)=u_star_interior(z); 
        z=z+1;
    end
end

%%Apply BCs
u_star(1,1:n+1) = 0; %left wall
u_star(n,1:n+1) = 0; %right wall
u_star(1:n, 1) = -u_star(1:n, 2); %bottom wall
u_star(1:n, n+1) = -u_star(1:n, n); %top wall 

z=1; 
for j=2:n-1
    for i=2:n
        v_star(i,j)=v_star_interior(z); 
        z=z+1;
    end
end

%%apply BCs
v_star(1,1:n) = -v_star(2,1:n); %left wall
v_star(n+1,1:n) = -v_star(n,1:n); %right wall
v_star(1:n+1, 1) = 0; %bottom wall
v_star(1:n+1, n) = 0; %top wall

%%Inlet and Outlet
v_star(2*(n_inlet):3*(n_inlet)-2,1)=v1;             %Inlet
v_star((n+1)-(n_inlet-1):n,n)=v_star((n+1)-(n_inlet-1):n,n-1);                   %Outlet

return
end


