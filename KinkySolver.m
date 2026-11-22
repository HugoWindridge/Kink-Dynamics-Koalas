% size of domain (goes from -L/2 to L/2)
L = 20;

%Number of grid points.
N = 250; 
dx = L/(N-1);

% need grid from -L to L

xgrid = linspace(-L/2, L/2, N);


%Taken from Demo_Shnack1D.m:
% Laplacian matrix
% Construct the Laplacian matrix using finite difference method
one_vec = ones(N,1);
lap = spdiags([one_vec, -2*one_vec, one_vec],[1,0,-1],N,N);
% Von Neumann Boundary Conditions
lap(1,1) = -1; lap(N,N) = -1;
lap = (1/dx^2)*lap;

% This is dU/dphi: (p is phi)
f = @(p)(2  * p .* (1 - p.^2));

% This is system to solve:
% p(1:N) is phi, p(2:N) is dphi/dt

% 1:N dphi/dt, N+1 : 2N d^2phi/dt
% this says: d/dt (phi) = dphi/dt (obviously), and
% d/dt (dphi/dt) = right hand side of eqn
F = @(t, p)[p(N+1:end); f(p(1:N)) + lap*p(1:N)];

initcond = zeros(2*N, 1); 
%initcond(1:N) = tanh(xgrid); % phi initial condition of tanh(x), speed dphi/dt = 0 initially (stationary)
initialcond(1:floorN/2)

odeparams = odeset('RelTol',1e-10,'AbsTol',1e-10);

T = linspace(0,100,1000)


%plot(xgrid, initcond(:, 1))

[T, phi] = ode45(F,T,initcond, odeparams);

size(phi)

plot(xgrid, phi(1, 1:N))
hold on
plot(xgrid, phi(end, 1:N), ":")

legend('Initial Conditions', 'Final Solution')
