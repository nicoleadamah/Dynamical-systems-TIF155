close all
clear all 
clc
% 3.3 Lyapunov exponents for the Lorenz model
% Lorentz attractor

% 3.3b
sigma = 10;
b = 8/3;
r = 28;
dt = 1e-3;


% 3.3c & 3.3d
% sigma = 10;
% b = 19/6;
% r = 28;
% dt = 1e-3;

% 3.3e & 3.3f
% sigma = 16;
% b = 5;
% r = 330;
% dt = 1e-4;

% Initialize and solve eq. 
t0 = 0;
t = 2;
xr1 = rand(1,1); 
xr2 = rand(1,1); 
xr3 = rand(1,1); 
x0=[xr1 xr2 xr3]; 
f = @(t, y)[sigma*(y(2)-y(1)); r*y(1)-y(2)-y(1)*y(3); y(1)*y(2)-b*y(3)];
[~, y] = ode45(f, linspace(t0, t, t_steps), x0);
Q = eye(3);
lambda = zeros(t_steps^2, 3);
t_steps = 1/dt;

for j = 1:t_steps
  for i = 1:t_steps

    % Compute new M(deformation matrix)
    jacobian = [-sigma, sigma, 0; r-y(i,3), -1, -y(i,1);y(i,2), y(i,1), -b];
    M = eye(3)+jacobian.*dt;
    % Q = Qnew
    [Q,R] = qr(M*Q);
    lambda(i+1 + t_steps*(j-1), :) = lambda(i+t_steps*(j-1),:) + log(abs(diag(R)))';

  end
  
  y0(1) = y(end - 1, 1);
  y0(2) = y(end - 1, 2);
  y0(3) = y(end - 1, 3);
  t0 = t;
  t = t + 1;
  
  [~, y] = ode45(f, linspace(t0, t, t_steps), y0);
  
  if rem(j, 100)== 0
    disp(j)
  end
    
end
%% PLOT AND DISPLATY LAMBDAS
L = lambda(end, :)./((t_steps^2+1)'.*dt);
disp(L)
figure;
semilogx(1:t_steps^2+1, lambda./((1:t_steps^2+1)'.*dt))
legend("{\lambda_1} = "+num2str(L(1)), "{\lambda_2} = "+num2str(L(2)), "{\lambda_3} = "+num2str(L(3)))
xlabel('Time (logarithmic scaling)'); 
title('Evolution of Lyapunov Exponents over Time with '+ "{\sigma} = "+ num2str(sigma) +" b = "+ num2str(b) +" r = "+ num2str(r));
grid on;