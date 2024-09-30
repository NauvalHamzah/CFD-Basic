function [bv] = get_rhsV(n,n_inlet,rho,dx,v,v_star,p,d_v,alpha,mu)

bv=zeros((n-1)*(n-2),1);    %vector of RHS for solving pressure corrections
stride = n-1;


for j=2:n-1
    for i=2:n
        position = (i-1) + (j-2)*stride;
        pressure_term = (p(i,j)-p(i,j+1)) * dx;
        bv(position,1) = pressure_term+(1-alpha)*(dx*v(i,j)/d_v(i,j));         
    end
end

Ds  = mu; 
A = @(F,D)( max(0, (1-0.1 * abs(F/D))^5 ) );


for i=1:(n_inlet-1)
    row = 2*(n_inlet)+(i-1);
    pos = 2*(n_inlet-1)+i;
    Fs  = .5*rho*dx*(v(row,1)+v(row,2));
    aS = Ds * A(Fs,Ds) + max(Fs,0);    
    bv(pos,1) = bv(pos)+aS*v_star(row,1);
end

Dn  = mu;
for i=1:(n_inlet-1)
    row = n-(i-1);
    pos = (n-1)*(n-2)-(i-1);
    Fn  = .5*rho*dx*(v(i,j)+v(i,j+1));
    aN = Dn * A(Fn,Dn) + max(-Fn,0);    
    bv(pos,1) = bv(pos)+aN*v_star(row,n);
end

return 
end
