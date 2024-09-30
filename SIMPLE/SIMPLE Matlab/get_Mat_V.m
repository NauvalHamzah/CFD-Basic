function [Av,d_v] = get_Mat_V(n,dx,dy,rho,mu,u,v,alpha)

N = (n-1)*(n-2);
stride = n-1;
Av = zeros(N,N);

%U-star for the boundary nodes is set to zero
%the interior nodes (n-2)*(n-1) solved implicitly for U-star

d_v=zeros(n+1,n);

De  = mu*dy / dx;  %diffusive coefficients
Dw  = mu*dy / dx; 
Dn  = mu*dx / dy; 
Ds  = mu*dx / dy; 

A = @(F,D)( max(0, (1-0.1 * abs(F/D))^5 ) );

for j=2:n-1
    for i=2:n
        position = (i-1) + (j-2)*stride;
        
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
        
        
        %set BSc for four corners
        if(i == 2 && j == 2)
            Av(position,position+1) = -aE;                                     
            Av(position,position+stride) = -aN;       
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Av(position,position) = aP/alpha+aW;
            continue;
        end
        if (i == n && j == 2)
            Av(position,position-1) = -aW;                    
            Av(position,position+stride) = -aN;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Av(position,position) = aP/alpha+aE; 
            continue;
        end
        if (i == 2 && j == n-1)
            Av(position,position+1) = -aE;
            Av(position,position-stride) = -aS;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Av(position,position) = aP/alpha+aW;                   
            continue;
        end
        if (i == n && j == n-1)
            Av(position,position-1) = -aW;
            Av(position,position-stride) = -aS;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Av(position,position) = aP/alpha+aE;
            continue;
        end
        %set four boundaries
        if (i == 2)
            Av(position,position+1) = -aE;
            Av(position,position+stride) = -aN;
            Av(position,position-stride) = -aS;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Av(position,position) = aP/alpha+aW;
            continue;
        end
        if (j == 2)  
            Av(position,position+1) = -aE;
            Av(position,position+stride) = -aN;
            Av(position,position-1) = -aW;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Av(position,position) = aP/alpha;
            continue;
        end
        if (i == n)
            Av(position,position+stride) = -aN;
            Av(position,position-stride) = -aS;
            Av(position,position-1) = -aW;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Av(position,position) = aP/alpha+aE;
            continue;
        end
        if (j == n-1 )
            Av(position,position+1) = -aE;
            Av(position,position-stride) = -aS;
            Av(position,position-1) = -aW;
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            Av(position,position) = aP/alpha;
            continue;
        end
        
        % interior nodes

        Av(position,position-1) = -aW;                      %sub diagonal
        Av(position,position+1) = -aE;                        %%upper diagonal
        Av(position,position-stride) = -aS;                 %%sub sub diagonal
        Av(position,position+stride) = -aN;                   %%upper upper diagonal       
        
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        Av(position,position) = aP/alpha;                                        %%main diagonal

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


return
end