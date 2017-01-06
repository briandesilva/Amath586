function [h,k,error] = kdv_modified_CN(m,ax,bx)
%
% kdv_CN.m
%
% Solve u_t+6uu_x+ u_{xxx}=0 on [ax,bx] with Dirichlet boundary conditions,
% using fractional step with Crank-Nicolson and first order upwind methods
% with m interior points.
%
% Returns k, h, and the max-norm of the error.
% This routine can be embedded in a loop on m to test the accuracy,
% perhaps with calls to error_table and/or error_loglog.


% clf              % clear graphics
% hold on          % Put all plots on the same graph (comment out if desired)

alpha = 10;               % parameter in initial data
gamma = 1;             % parameter appearing in PDE: u_t + gamma u_{xxx}=0
tfinal = 0.02;                % final time

h = (bx-ax)/(m+1);         % h = delta x
x = linspace(ax,bx,m+2)';  % note x(1)=ax and x(m+2)=bx
                           % u(1)=g0 and u(m+2)=g1 are known from BC's
                           % With periodic BCs there are m+1 unknowns
I1 = 1:(m+1);        % Indices of unknowns (identify x_{m+2} with x_0)
I2 = 2:(m+2);         % Indices of unknowns with ghost cells
k = h/300;                  % time step

nsteps = round(tfinal / k);    % number of time steps
% nplot = 1;      % plot solution every nplot time steps
                 % (set nplot=2 to plot every 2 time steps, etc.)
nplot = nsteps;  % only plot at final time

if abs(k*nsteps - tfinal) > 1e-5
   % The last step won't go exactly to tfinal.
   disp(' ')
   disp(sprintf('WARNING *** k does not divide tfinal, k = %9.5e',k))
   disp(' ')
end

eta = @(x) (alpha^2) * (sech(alpha*x/2).^2) / 2;

% true solution for comparison:
% For sech initial data
utrue = @(x,t) eta(x-t*alpha^2);

% initial conditions:
u0 = eta(x);


% Each time step we first solve MOL system U' = AU using the Trapezoidal
% method, use this to solve u_t + u_xxx = 0

% Set up matrices for method (3):
e = ones(m+1,1);
A = spdiags([-e 3*e -3*e e], [-1 0 1 2], m+1, m+1);
A(end,1) =-3;
A(end,2) = 1;
A(end-1,1) = 1;
A(1,end) = -1;
A = -(gamma*k/(2*h^3)) * A;


% Precompute inv(I-A)*A
M = (eye(m+1)-A) \ (A+eye(m+1));

% Set up matrix for upwind method
B = spdiags([-e e], [-1 0], m+1, m+1);
B(1,end) = -1;
B = (-3*k/h) * B;


% initial data on fine grid for plotting:
xfine = linspace(ax,bx,1001);
% ufine = utrue(xfine,0);

% initialize u and plot:
tn = 0;
u = u0;

% periodic boundary conditions:
u(m+2) = u(1);   % copy value from ax to bx

% ghost cell on left
% u(1) = u(m+2);

% plot(x,u,'b.-', xfine,ufine,'r')
% legend('computed','true')
% title('Initial data at time = 0')
% 
% input('Hit <return> to continue  ');


% main time-stepping loop:

for n = 1:nsteps
     tnp = tn + k;   % = t_{n+1}

    % First solve u_t + u_{xxx} = 0 using method (3)
    % This gives u, intermediate solution
    
    % Just multiply by M to get next approximation 
    u(I1) = M*u(I1);
    
    % Periodic BCs
    u(m+2) = u(1);
    
    % Next solve u_t + 3(u^2)_x = 0 using upwind method
    u(I1) = u(I1) + B * (u(I1).^2);      % (u^2)_x
%     u(I1) = u(I1) + 2*u(I1).*(B*u(I1));     % 2uu_x
%     u(I2) = u(I2) - (3*k/h)*(u(I2).^2 - u(I2-1).^2);
    
    % Periodic BCs
     u(m+2) = u(1);

    % Ghost cell
%     u(1) = u(m+2);
    
     % plot results at desired times:
     if mod(n,nplot)==0 || n==nsteps
        ufine = utrue(xfine,tnp);
        plot(x,u,'b.-', xfine,ufine,'r')
        title(sprintf('t = %9.5e  after %4i time steps with %5i grid points',...
                       tnp,n,m+2))
        error = max(abs(u-utrue(x,tnp)));
%         disp(sprintf('at time t = %9.5e  max error =  %9.5e',tnp,error))
%         if n<nsteps
%             input('Hit <return> to continue  '); 
%         end
     end

     tn = tnp;   % for next time step
end

end
