function [Au,d_u] = get_Mat_U(n,dx,dy,rho,mu,u,v,alpha)

N = (n-2)*(n-1);
stride = n-2;
Au = zeros(N,N);

%U-star for the boundary nodes is set to zero
%the interior nodes (n-2)*(n-1) solved implicitly for U-star

d_u=zeros(n,n+1);

De  = mu*dy / dx;  %diffusive coefficients
Dw  = mu*dy / dx; 
Dn  = mu*dx / dy; 
Ds  = mu*dx / dy; 

A = @(F,D)( max(0, (1-0.1 * abs(F/D))^5 ) );

for j=2:n
    for i=2:n-1
        position = (i-1) + (j-2)*stride; 
        
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
        
        %set BSc for four corners
        if(i == 2 && j == 2)
            Au(position,position+1) = -aE;                                     
            Au(position,position+stride) = -aN;       
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Au(position,position) = aP/alpha+aS;
            continue;
        end
        if (i == n-1 && j == 2)
            Au(position,position-1) = -aW;                    
            Au(position,position+stride) = -aN;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Au(position,position) = aP/alpha+aS; 
            continue;
        end
        if (i == 2 && j == n)
            Au(position,position+1) = -aE;
            Au(position,position-stride) = -aS;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Au(position,position) = aP/alpha+aN;                   
            continue;
        end
        if (i == n-1 && j == n)
            Au(position,position-1) = -aW;
            Au(position,position-stride) = -aS;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Au(position,position) = aP/alpha+aN;
            continue;
        end
        %set four boundaries
        if (i == 2)
            Au(position,position+1) = -aE;
            Au(position,position+stride) = -aN;
            Au(position,position-stride) = -aS;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Au(position,position) = aP/alpha;
            continue;
        end
        if (j == 2)
            Au(position,position+1) = -aE;
            Au(position,position+stride) = -aN;
            Au(position,position-1) = -aW;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Au(position,position) = aP/alpha+aS;
            continue;
        end
        if (i == n-1)
            Au(position,position+stride) = -aN;
            Au(position,position-stride) = -aS;
            Au(position,position-1) = -aW;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Au(position,position) = aP/alpha;
            continue;
        end
        if (j == n )
            Au(position,position+1) = -aE;
            Au(position,position-stride) = -aS;
            Au(position,position-1) = -aW;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Au(position,position) = aP/alpha+aN;
            continue;
        end
        
        % interior nodes

        Au(position,position-1) = -aW;                      %sub diagonal
        Au(position,position+1) = -aE;                        %%upper diagonal
        Au(position,position-stride) = -aS;                 %%sub sub diagonal
        Au(position,position+stride) = -aN;                   %%upper upper diagonal       
        
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        Au(position,position) = aP/alpha;                                        %%main diagonal

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


return
end