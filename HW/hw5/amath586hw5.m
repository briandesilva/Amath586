%% Amath 586 HW 5 Matlab Notebook
% Brian de Silva


%% Problem 2

% Set the relevant parameters
a = 1;
h = 1/50;
k = 0.8*h;
p = [0:round(1/h)]';

shifted_unit_circle = [cos(0:0.001:2*pi) - 1; sin(0:0.001:2*pi) ].';


% Part (a)
epsilon = a*h/2;

% Find eigenvalues * k
kmu = k*[(-2*epsilon/(h^2))*(1-cos(2*pi*p*h)), (-a/h)*sin(2*pi*p*h)];

% Plot results
Fig1 = figure(1);
plot(shifted_unit_circle(:,1),shifted_unit_circle(:,2),'k-'), hold on
plot(kmu(:,1),kmu(:,2), 'bo');
axis([-2.5 0.5 -1.5 1.5]);
title('eps = 0.01 (Upwind) ')
grid on

 
% Part (b)
epsilon = -a*h/2;

% Find eigenvalues * k
kmu = k*[(-2*epsilon/(h^2))*(1-cos(2*pi*p*h)), (-a/h)*sin(2*pi*p*h)];

% Plot results
Fig2 = figure(2);
plot(shifted_unit_circle(:,1),shifted_unit_circle(:,2),'k-'), hold on
plot(kmu(:,1),kmu(:,2), 'bo');
axis([-2.1 2 -2.05 2.05]);
title('eps = -0.01 (Downwind) ')
grid on


%% Problem 5
close all

% Part (a) - Leapfrog
IC = 1;
numTrials = 8;

hVec = zeros(numTrials,1);
errorVec = zeros(numTrials,1);

% advection_Leap_pbc(239);

% Approximate PDE using increasing number of grid points (160,200,280,...)
% Note: we choose k=0.4h for these trials.
for j=1:numTrials
    [hVec(j),~,errorVec(j)] = advection_Leap_pbc(120+40*2^(j-1) - 1,IC);
end

% Create loglog plot of error as grid is refined
error_loglog(hVec,errorVec);

% Computed order of accuracy: 1.99
% Constant is O(10^3)


% Part (b) - New initial data
IC = 2;
numTrials = 8;

hVec = zeros(numTrials,1);
errorVec = zeros(numTrials,1);


% Approximate PDE using increasing number of grid points (200,240,320,...)
% Note: we choose k=0.4h for these trials.
for j=1:numTrials
    [hVec(j),~,errorVec(j)] = advection_Leap_pbc(160+40*2^(j-1) - 1,IC);
end

% Create loglog plot of error as grid is refined
error_loglog(hVec,errorVec);

% Computed order of accuracy: 1.96
% Constant is O(10^4)


% Part (c) 
IC = 3;

% Note: I had advection_Leap_pbc output multiple
% plots when I ran this line originally
[h,k,~]=advection_Leap_pbc(199,IC);















