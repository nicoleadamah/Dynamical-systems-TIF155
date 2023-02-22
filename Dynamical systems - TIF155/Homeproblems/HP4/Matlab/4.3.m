%Nicole Adamah
%2022
%4.3a
close all
clear all 
clc
a = 1.4;
b = 0.3;
nrPoints = 1e6;
nrsteps = 1e5;
nonValid = true;
while nonValid
  x = zeros(nrPoints, 1);
  y = zeros(nrPoints, 1);
  x(1) = rand;
  y(1) = rand;
  % Initialize -> close to attractor
  for i = 2:nrsteps
    x(2) = y(1) + 1 - a*x(1)^2;
    y(2) = b*x(1);
    x(1) = x(2);
    y(1) = y(2);
  end
  
  % Step 3 - Check that x,t -> inf
  if ~isinf(x(1)) && ~isinf(y(1))
    nonValid = false;
    
    for i = 2:nrPoints
      x(i) = y(i-1) + 1 - a*x(i-1)^2;
      y(i) = b*x(i-1);
    end
  end
  
end
figure()
plot(x,y,'.','MarkerSize',4,'Color','b')
title('Approximation of the fractal attractor(Hénon map), a = 1.4, b = 0.3', 'FontSize', 15)
xlabel('x','FontSize', 15)
ylabel('y','FontSize', 15)

%% 4.3b
steps = 1000;
eps = linspace(2e-2, 1e-3, steps);
I0 = zeros(steps, 1);
I1 = zeros(steps, 1);
I2 = zeros(steps, 1);
xmin = round(min(x), 1);
ymin = round(min(y), 1);
xmax = round(max(x), 1);
ymax = round(max(y), 1);

for i = 1:steps
  e = eps(i);
  % Compute P(the occupancy matrix)
  P = zeros(ceil((xmax - xmin)/e), ceil((ymax - ymin)/e));
  for j = 1:nrPoints
    xPos = ceil((x(j) - xmin)/e);
    yPos = ceil((y(j) - ymin)/e);
    P(xPos, yPos) = P(xPos, yPos) + 1;
  end
  
  P = P/sum(sum(P));
  P = P(P ~= 0);
  I0(i, 1) = sum(P.^0);
  I1(i, 1) = sum((1./P).^P); 
  I2(i, 1) = sum(P.^2);
end
%%
figure;
subplot(1, 3, 1);
plot(log(1./eps), log(I0), 'b');
xlabel('ln(1/\epsilon)')
ylabel('ln(I_q)/(1-q)')
title('q = 0')

subplot(1, 3, 2);
plot(log(1./eps), log(I1), 'b');
xlabel('ln(1/\epsilon)')
ylabel('ln(I_q)')
title('q = 1')

subplot(1, 3, 3);
plot(log(1./eps), -log(I2), 'b');
xlabel('ln(1/\epsilon)')
ylabel('ln(I_q)/(1-q)')
title('q = 2')

%% 4.3c
% Dq for q = [0, 1, 2]-> the slope of the graphs;
D = zeros(3, 2);
Iq = [I0, I1, I2];

D(1, :) = polyfit(log(1./eps)', log(Iq(:, 1)), 1);
D(2, :) = polyfit(log(1./eps)', log(Iq(:, 2)), 1);
D(3, :) = polyfit(log(1./eps)', -log(Iq(:, 3)), 1);
slope = D(:, 1);
fprintf('D1 = %.2f\nD2 = %.2f\nD3 = %.2f\n', D(1, 1), D(2, 1), D(3, 1))
%% 4.3d)
nrq = 20;
eps2 = linspace(0.02,0.001,nrq);
q = linspace(0,4,nrq);
I = zeros(length(eps2),length(q));
for j=1:length(q)
    for i=1:length(eps2)
        [N,Xedges,Yedges] = histcounts2(x,y,'BinWidth',[eps2(i) eps2(i)]);
        P2 = N./nrPoints;
        P2 = P2(P2 ~= 0);
        P2 = P2.^q(j);
        I(i,j) = sum(sum(P2));
    end
end
q2 = linspace(0,4,nrq);
remidx = floor(0.2*nrq);
p = zeros(length(q2),2);
Dq = zeros(length(q2),1);

for i =1:length(q2)
    tempE = log(1./eps2(1:(end-remidx)));
    tempI = log(I(1:(end-remidx),i))';
    p(i,:) = polyfit(tempE, tempI, 1);
    Dq(i) = (1/(1-q2(i))).*p(i,1);
end
figure()
plot(q2,Dq,'*','MarkerSize',10);
title("D_{q} as a function of q with 20 q-values",'FontSize',15);
xlabel("q",'FontSize',15);
ylabel("D_{q}",'FontSize',15);
