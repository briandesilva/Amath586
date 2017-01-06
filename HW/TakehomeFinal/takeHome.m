%% Amath 586 Take-home Final Matlab Notebook
% Brian de Silva

%% Problem 2 (c)

% Throughout we use p=2

numTrials =8;

% Create vectors in which to store h values and error
hVec = zeros(numTrials,1);
errorVec = zeros(numTrials,1);

% Approximate PDE using an increasing number of grid points (200,220,260,320,...)
% with k=sqrt(h)/2
for j=1:numTrials
    [hVec(j),~,errorVec(j)] = dispersion_modified_CN(200+20*2^(j-1) - 1,0,2*pi);
end

% Create loglog plot of error as grid is refined
error_loglog(hVec,errorVec);

% Computed order of accuracy: 0.93
% Note: Choosing different relationships between h and k gave approximate
% orders of accuracy ranging from 0.88-1.63.
% Overall this suggests first order accuracy

% Aside
%---------------------------------------------------------
% Just to be sure we made the right choice, let's try the method
% corresponding to using equation (2) instead of (3):

% Create vectors in whic to store h values and error
hVec = zeros(numTrials,1);
errorVec = zeros(numTrials,1);

% Approximate PDE using an increasing number of grid points (200,220,260,320,...)
% with k=sqrt(h)/2
for j=1:numTrials
    [hVec(j),~,errorVec(j)] = dispersion_modified_CN_worse(200+20*2^(j-1) - 1,0,2*pi);
end

% Create loglog plot of error as grid is refined
error_loglog(hVec,errorVec);
%---------------------------------------------------------


%% Problem 3

clear all

numTrials = 6;

% Create vectors in which to store h values and error
hVec = zeros(numTrials,1);
errorVec = zeros(numTrials,1);
tic
% Approximate PDE using an increasing number of grid points (600,800,1200,... )
% with k=h/300
for j=1:numTrials
    [hVec(j),~,errorVec(j)] = kdv_modified_CN(600+200*2^(j-1) - 1,-1,3);
    if isnan(errorVec(j))
        disp('NaN encountered')
        break
    end
end
toc
% Create loglog plot of error as grid is refined
error_loglog(hVec,errorVec);

% Computed order of accuracy: 0.91

