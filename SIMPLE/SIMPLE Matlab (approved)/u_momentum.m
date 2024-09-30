function [u_star,d_u] = u_momentum_new(n,dx,dy,rho,mu,u,v,p,alpha)

u_star=zeros(n,n+1);
d_u=zeros(n,n+1);

De  = mu*dy / dx;  %convective coefficients
Dw  = mu*dy / dx; 
Dn  = mu*dx / dy; 
Ds  = mu*dx / dy; 

A = @(F,D)( max(0, (1-0.1 * abs(F/D))^5 ) );

%%compute u_star
for i = 2:n-1
    for j = 2:n
        Fe  = .5*rho*dy*(u(i+1,j)+u(i,j));                                     
        Fw  = .5*rho*dy*(u(i-1,j)+u(i,j)); 
        Fn  = .5*rho*dx*(v(i+1,j)+v(i,j)); 
        Fs  = .5*rho*dx*(v(i+1,j-1)+v(i,j-1));
        
        aE = De * A(Fe,De) + max(-Fe,0);
        aW = Dw * A(Fw,Dw) + max(Fw,0);
        aN = Dn * A(Fn,Dn) + max(-Fn,0);
        aS = Ds * A(Fs,Ds) + max(Fs,0);
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        
        d_u(i,j) = alpha * dy / aP;   %refer to Versteeg CFD book
        
        pressure_term = (p(i,j)-p(i+1,j)) * dy;
        
        u_star(i,j) = (alpha/aP) * ( (aE*u(i+1,j)+aW*u(i-1,j)+aN*u(i,j+1)+aS*u(i,j-1)) + pressure_term ) + (1-alpha)*u(i,j);
                       
    end
end

%%set d_u for top and bottom BCs
%%they will be later used by the pressure correction equation 
%%they should not be zero, or BCs of pressure correction will get messed up
j = 1; %bottom
for i=2:n-1
    Fe  = .5*rho*dy*(u(i+1,j)+u(i,j));                                     
    Fw  = .5*rho*dy*(u(i-1,j)+u(i,j)); 
    Fn  = .5*rho*dx*(v(i+1,j)+v(i,j)); 
    Fs  = 0;
        
    aE = De * A(Fe,De) + max(-Fe,0);
    aW = Dw * A(Fw,Dw) + max(Fw,0);
    aN = Dn * A(Fn,Dn) + max(-Fn,0);
    aS = 0;
    aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
    d_u(i,j) = alpha * dy / aP;
end

j = n+1; %top
for i=2:n-1
    Fe  = .5*rho*dy*(u(i+1,j)+u(i,j));                                     
    Fw  = .5*rho*dy*(u(i-1,j)+u(i,j)); 
    Fn  = 0; 
    Fs  = .5*rho*dx*(v(i+1,j-1)+v(i,j-1));
        
    aE = De * A(Fe,De) + max(-Fe,0);
    aW = Dw * A(Fw,Dw) + max(Fw,0);
    aN = 0;
    aS = Ds * A(Fs,Ds) + max(Fs,0);
    aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
    d_u(i,j) = alpha * dy / aP;
end


%%Apply BCs
u_star(1,1:n+1) = 0; %left wall
u_star(n,1:n+1) = 0; %right wall
u_star(1:n, 1) = -u_star(1:n, 2); %bottom wall
u_star(1:n, n+1) = -u_star(1:n, n); %top wall 


return 
end