%%LID DRIVEN CAVITY FLOW EXAMPLE
%%The code is based on the SIMPLE algorithm

clc
clear all
close all
format longG

%% GRID SIZE AND OTHER PARAMETERS
%i runs along x-direction and j runs along y-direction 

n_inlet=6;                        
n=(n_inlet-1)*5+1;              %grid size 
max_iteration=10000; 
maxRes = 1000;
iteration = 1;

rho = 1.225;                %density
v1=0.001;                   %velocity = lid velocity
dx=1/(n-1);                 %dx,dy cell sizes along x and y directions
dy=1/(n-1); 
mu = 1.7894e-5;             %viscosity
x=dx/2:dx:1-dx/2;
y=0:dy:1; 
alphaP = 0.1;               %pressure under-relaxation
alphaU = 0.9;               %velocity under-relaxation
tol = 1e-6;

%   u_star, v_star are Intermediate velocities
%   u and v = Final velocities

%Variable declaration
p   = zeros(n+1,n+1);             %   p = Pressure
p_star   = zeros(n+1,n+1);        
p_prime = zeros(n+1,n+1);         %   pressure correction 
p_final = zeros(n+1,n+1);
rhsp = zeros((n-1)*(n-1),1);            %   Right hand side vector of pressure correction equation
divergence = zeros(n,n); 

%Vertical velocity
v_star  = zeros(n+1,n);
vold    = zeros(n+1,n);
vRes    = zeros(n+1,n);
v       = zeros(n+1,n);
v_final = zeros(n,n);
d_v     = zeros(n+1,n);      %velocity correction coefficient

% Horizontal Velocity -----------
u_star = zeros(n,n+1);
uold   = zeros(n,n+1);
uRes   = zeros(n,n+1);
u      = zeros(n,n+1);
u_final = zeros(n,n);
d_u    = zeros(n,n+1);      %velocity correction coefficient

velocity_final = zeros(n,n);


%Boundary condition 
%Inlet and Outlet
v_star(2*(n_inlet):3*(n_inlet)-2,1)=v1;            %Inlet
v(2*(n_inlet):3*(n_inlet)-2,1)=v1;                 %Inlet
v_star((n+1)-(n_inlet-1):n,n)=v_star((n+1)-(n_inlet-1):n,n-1);                        %Outlet
v((n+1)-(n_inlet-1):n,n)=v((n+1)-(n_inlet-1):n,n-1);                             %Outlet


%% ---------- iterations -------------------
while ( (iteration <= max_iteration) && (maxRes > tol) ) 
    iteration = iteration + 1;
    [u_star,d_u] = u_momentum(n,dx,dy,rho,mu,u,v,p_star,alphaU);       %%Solve u-momentum equation for intermediate velocity u_star 
    [v_star,d_v] = v_momentum(n,n_inlet,dx,dy,rho,mu,u,v,p_star,v1,alphaU);                 %%Solve v-momentum equation for intermediate velocity v_star
    uold = u;
    vold = v; 
    [rhsp] = get_rhs(n,n_inlet,dx,dy,rho,u_star,v_star);                                 %%Calculate rhs vector of the Pressure Poisson matrix 
    [Ap] = get_coeff_mat_modified(n,n_inlet,dx,dy,rho,d_u,d_v);                          %%Form the Pressure Poisson coefficient matrix 
    [p,p_prime] = pres_correct(n,n_inlet,rhsp,Ap,p_star,alphaP);                         %%Solve pressure correction implicitly and update pressure
    [u,v] = updateVelocity(n,n_inlet,u_star,v_star,p_prime,d_u,d_v,v1);            %%Update velocity based on pressure correction
    %[divergence]=checkDivergenceFree(n,dx,dy,u,v);                               %%check if velocity field is divergence free
    p_star = p;                                                                          %%use p as p_star for the next iteration
    
    %find maximum residual in the domain
%     vRes = abs(v - vold);
%     uRes = abs(u - uold);
%     maxRes_u = max(max(uRes));
%     maxRes_v = max(max(vRes));
%     maxRes = max(maxRes_u, maxRes_v);
     error=0;
     for i=1:(n-1)*(n-1)
         error = error+abs(rhsp(i,1));
     end
     
     maxRes=error;                                                                        %%Check for convergence 
    
     iter(iteration)=iteration;
     All_error(iteration)=error;
     disp(['It = ',int2str(iteration),'; Res = ',num2str(maxRes)])
    if (maxRes > 2)
        disp('not going to converge!');
        break;
    end
end

%% plot

[u_final,v_final,p_final,velocity_final] = FinalMapping(u,v,p,n,n_inlet);

disp(['Total Iterations = ',int2str(iteration)])

figure (21)
contourf(y,y,velocity_final',50, 'edgecolor','none');colormap jet
colorbar;
axis([0 1 0 1]); 
title('Velocity'); 

[X,Y] = meshgrid(y,y);
figure (21)
hold on
quiver(X, Y, u_final', v_final', 5, 'k')

figure (22)
contourf(y,y,p_final',50, 'edgecolor','none');colormap jet
colorbar;
axis([0 1 0 1]); 
title('Pressure');

figure (23)
plot(iter,All_error)
grid on


