function [Ap] = get_coeff_mat_modified(n,n_inlet,dx,dy,rho,d_u,d_v)

N = (n-1)*(n-1);
stride = n-1;
Ap = zeros(N,N);

%P-prime for the boundary nodes is set to zero
%the interior nodes (n-1)*(n-1) solved implicitly for P-prime

for j=2:n
    for i=2:n
        position = (i-1) + (j-2)*stride; 
        aE = 0;
        aW = 0;
        aN = 0;
        aS = 0;
        
        %set BSc for four corners
        if(j==n && i>(n-(n_inlet-1)))
            Ap(position,position) = 1;
            continue;
        end
        if(i == 2 && j == 2)
            Ap(position,position+1) = -rho*d_u(i,j)*dy;                    %pressure correction at the first node is zero                  
            aE = -Ap(position,position+1);
            Ap(position,position+stride) = -rho*d_v(i,j)*dx;       
            aN = -Ap(position,position+stride);
            aP = aE + aN + aW + aS;
            Ap(position,position) = aP;
            continue;
        end
        if (i == n && j == 2)
            Ap(position,position-1) = -rho*d_u(i-1,j)*dy;                    
            aW = -Ap(position,position-1);
            Ap(position,position+stride) = -rho*d_v(i,j)*dx;       
            aN = -Ap(position,position+stride);
            aP = aE + aN + aW + aS;
            Ap(position,position) = aP; 
            continue;
        end
        if (i == 2 && j == n)
            Ap(position,position+1) = -rho*d_u(i,j)*dy;                    
            aE = -Ap(position,position+1);
            Ap(position,position-stride) = -rho*d_v(i,j-1)*dx;  
            aS = -Ap(position,position-stride);
            aP = aE + aN + aW + aS;
            Ap(position,position) = aP;                   
            continue;
        end
        if (i == n && j == n)
            Ap(position,position-1) = -rho*d_u(i,j)*dy;                    
            aW = -Ap(position,position-1);
            Ap(position,position-stride) = -rho*d_v(i,j)*dx;       
            aS = -Ap(position,position-stride);
            aP = aE + aN + aW + aS;
            Ap(position,position) = aP;
            continue;
        end
        %set four boundaries
        if (i == 2)
            Ap(position,position+1) = -rho*d_u(i,j)*dy;                    
            aE = -Ap(position,position+1);
            Ap(position,position+stride) = -rho*d_v(i,j)*dx;       
            aN = -Ap(position,position+stride);
            Ap(position,position-stride) = -rho*d_v(i,j-1)*dx;       
            aS = -Ap(position,position-stride);
            aP = aE + aN + aW + aS;
            Ap(position,position) = aP;
            continue;
        end
        if (j == 2)
            Ap(position,position+1) = -rho*d_u(i,j)*dy;                    
            aE = -Ap(position,position+1);
            Ap(position,position+stride) = -rho*d_v(i,j)*dx;       
            aN = -Ap(position,position+stride);
            Ap(position,position-1) = -rho*d_u(i-1,j)*dy;                    
            aW = -Ap(position,position-1);
            aP = aE + aN + aW + aS;
            Ap(position,position) = aP;
            continue;
        end
        if (i == n)
            Ap(position,position+stride) = -rho*d_v(i,j)*dx;       
            aN = -Ap(position,position+stride);
            Ap(position,position-stride) = -rho*d_v(i,j-1)*dx;       
            aS = -Ap(position,position-stride);
            Ap(position,position-1) = -rho*d_u(i-1,j)*dy;                    
            aW = -Ap(position,position-1);
            aP = aE + aN + aW + aS;
            Ap(position,position) = aP;
            continue;
        end
        if (j == n )
            Ap(position,position+1) = -rho*d_u(i,j)*dy;                    
            aE = -Ap(position,position+1);
            Ap(position,position-stride) = -rho*d_v(i,j-1)*dx;       
            aS = -Ap(position,position-stride);
            Ap(position,position-1) = -rho*d_u(i-1,j)*dy;                    
            aW = -Ap(position,position-1);
            aP = aE + aN + aW + aS;
            Ap(position,position) = aP;
            continue;
        end
        
        % interior nodes

        Ap(position,position-1) = -rho*d_u(i-1,j)*dy;                      %sub diagonal
        aW = -Ap(position,position-1);

        Ap(position,position+1) = -rho*d_u(i,j)*dy;                        %%upper diagonal
        aE = -Ap(position,position+1);

        Ap(position,position-stride) = -rho*d_v(i,j-1)*dx;                 %%sub sub diagonal
        aS = -Ap(position,position-stride);

        Ap(position,position+stride) = -rho*d_v(i,j)*dx;                   %%upper upper diagonal
        aN = -Ap(position,position+stride);
        
        aP = aE + aN + aW + aS;
        Ap(position,position) = aP;                                        %%main diagonal
    end
end

return
end