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
% Von Neumann Boundary Conditions (Carl Neumann?) - Yeah...
lap(1,1) = -1; lap(N,N) = -1;
lap = (1/dx^2)*lap;

% This is -dU/dphi: (p is phi)
f = @(p)(2  * p .* (1 - p.^2));

% This is system to solve:
% p(1:N) is phi, p(N+1:2N) is dphi/dt

% 1:N dphi/dt, N+1 : 2N d^2phi/dt
% this says: d/dt (phi) = dphi/dt (obviously), and
% d/dt (dphi/dt) = right hand side of eqn
%F = @(t, p)[p(N+1:end); f(p(1:N)) + lap*p(1:N)];
%absorbing bc's:
F = @(t, p)[(p(2) - p(1))/dx; p(N+2:end-1); -(p(N) - p(N-1))/dx; f(p(1:N)) + lap*p(1:N)];


Av0 = 0.1; 
Ax0 = -5; %move kink left a bit so it has room to move right.
Bv0 = 0.1;
Bx0 = 5;
Agamma = 1 / sqrt(1 - Av0^2); %Lorentz factor for kink A
Bgamma = 1 / sqrt(1 - Bv0^2); %Lorentz factor for kink B

%this is a stationary solution.
initcond = zeros(2*N, 1); 
initcond(1:N) = tanh(Agamma * (xgrid - Ax0)) - tanh(Bgamma * (xgrid - Bx0)) - 1; % phi initial condition of tanh(x), speed dphi/dt = 0 initially (stationary)
%initcond(1:floor(N/2))

%Investigating a moving solution involves giving it an initial velocity.


initcond(N+1:2*N) = -Av0 * Agamma * sech((Agamma * (xgrid - Ax0))).^2 -Bv0 * Bgamma * sech((Bgamma * (xgrid - Bx0))).^2; %Applying an inital velocity in the form of a small bump. 


odeparams = odeset('RelTol',1e-10,'AbsTol',1e-10);

T = linspace(0,300,1000);


%plot(xgrid, initcond(:, 1))

[T, phi] = ode45(F,T,initcond, odeparams);

size(phi)

%plot(xgrid, phi(1, 1:N))
%hold on
%plot(xgrid, phi(end, 1:N), ":")

framelen = 5; %Displaying 5 time step updates per frame, crank it to speed up the animation and put it 1 to see every step.

figure;
for k = 1:framelen:length(T)
    plot(xgrid, phi(k, 1:N), 'LineWidth', 2);
    ylim([-1.5 1.5]);
    xlabel('Position x');
    ylabel('\phi(x, t)');
    title(['Time t = ', num2str(T(k), '%.1f')]);
    grid on;
    drawnow;
end

legend('Initial Conditions', 'Final Solution')
