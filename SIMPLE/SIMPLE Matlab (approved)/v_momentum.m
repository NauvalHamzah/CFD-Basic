function [v_star,d_v] = v_momentum_new(n,n_inlet,dx,dy,rho,mu,u,v,p,v1,alpha)

v_star=zeros(n+1,n);
d_v=zeros(n+1,n);

De  = mu*dy / dx;  %convective coefficients
Dw  = mu*dy / dx; 
Dn  = mu*dx / dy; 
Ds  = mu*dx / dy; 

A = @(F,D)( max(0, (1-0.1 * abs(F/D))^5 ) );

%%compute u_star
for i = 2:n
    for j = 2:n-1
        Fe  = .5*rho*dy*(u(i,j)+u(i,j+1));                                     
        Fw  = .5*rho*dy*(u(i-1,j)+u(i-1,j+1)); 
        Fn  = .5*rho*dx*(v(i,j)+v(i,j+1)); 
        Fs  = .5*rho*dx*(v(i,j-1)+v(i,j));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        
        d_v(i,j) = alpha * dx / aP;   %refer to Versteeg CFD book
        
        pressure_term = (p(i,j)-p(i,j+1)) * dx;
        
        v_star(i,j) = (alpha/aP) * ( (aE*v(i+1,j)+aW*v(i-1,j)+aN*v(i,j+1)+aS*v(i,j-1)) + pressure_term ) + (1-alpha)*v(i,j);
               
    end
end

%%set d_v for left and right BCs
%%they will be later used by the pressure correction equation 
%%they should not be zero, or BCs of pressure correction will get messed up
%%Apply BCs
i = 1;  %left BC
for j=2:n-1
    Fe  = .5*rho*dy*(u(i,j)+u(i,j+1));                                     
    Fw  = 0; 
    Fn  = .5*rho*dx*(v(i,j)+v(i,j+1)); 
    Fs  = .5*rho*dx*(v(i,j-1)+v(i,j));                                                       
        
    aE = De * A(Fe,De) + max(-Fe,0);
    aW = 0;
    aN = Dn * A(Fn,Dn) + max(-Fn,0);
    aS = Ds * A(Fs,Ds) + max(Fs,0);
    aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
    d_v(i,j) = alpha * dx / aP;
end

i = n+1;  %right BC
for j=2:n-1
    Fe  = 0;                                     
    Fw  = .5*rho*dy*(u(i-1,j)+u(i-1,j+1)); 
    Fn  = .5*rho*dx*(v(i,j)+v(i,j+1)); 
    Fs  = .5*rho*dx*(v(i,j-1)+v(i,j));
        
    aE = 0;
    aW = Dw * A(Fw,Dw) + max(Fw,0);
    aN = Dn * A(Fn,Dn) + max(-Fn,0);
    aS = Ds * A(Fs,Ds) + max(Fs,0);
    aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
    d_v(i,j) = alpha * dx / aP;
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