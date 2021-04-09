% LayerOptimizer Tuning
%
% Jon Cooper

set(0, 'DefaultLineLineWidth', 2);
close all
rng(0);

%% Artificial Models

zd = (0:.002:1)';
obs = zeros(length(zd),2); obs(:,1) = zd;

% LOO
% pos = 0;
% res = 1;
% AOI
% pos = [0;.4;.55;.9];
% res = [1;2;1;2];
% Initial guess / Coordinate Descent
% pos = [0;.2;.26;.32;.38;.44;.50]; % CHECK IN FUTURE
% res = [1;2;1;2;1;2;1];
pos = [0;.4;.44;.48;.52;.56;.60];
res = [1;2;1;2;1;2;1];

qcTrue = LayerModelEval([pos;res],zd);
% plot(qcTrue,-zd);

%% Setup
% kern = normpdf(linspace(-3,3,65)); % use a basic bell curve for the psf
kern = chi2pdf(linspace(0,8,65),4);
kern = kern/sum(kern);
% N = 20; % initial number of layers

%% Force Blur Measures
M = size(obs,1); mask = eye(M); n = 0; z = 0;
mask = [zeros(M,floor(length(kern)/2)),mask,zeros(M,ceil(length(kern)/2)-1)];
mask = diag([zeros(1,z),linspace(0,1,n),ones(1,M-n-z)].*[ones(1,M-n-z),linspace(1,0,n),zeros(1,z)])*mask;
qcMeas = mask*conv(qcTrue,kern);
obs(:,2) = qcMeas;

%% Layer0 Choice
% LOO
% pos = [0;.2;.3;.6];
% res = [1;1.5;1;1.5];
% AOI
% pos = [0;.6];
% res = [1;1];
% Initial guess / Coordinate Descent
pos = []; res = [];

layers0 = [pos;res];
options = {'CD',1;'MIN',.04;'MAX',3;'Plt',3;'MaxIter',0};
options = [options;{'layers0',layers0}]; N = length(pos);
% N = 0;

%% Optimize
[layers,info] = LayerOptimizer(N,obs,kern,options);

%% Plot
layers0 = info.layers0; N = length(layers0)/2;
Nf = info.Nfinal;
IM = info.InitialMisfit;
FM = info.FinalMisfit;
mask = info.mask;

fprintf('Initial Misfit: %f\n',IM);
fprintf('Final Misfit: %f\n',FM);

figure
% subplot(1,2,1);
t0 = LayerModelEval(layers0,zd);
plot(t0,-zd,'-r'), hold on, plot(layers0(N+1:2*N),-layers0(1:N),'or'), plot(mask*conv(t0,kern),-zd,'--r');
plot(obs(:,2),-zd,'--b'), plot(qcTrue,-zd,'-b'), hold off;
legend('q_c^{inv}','Top of Layer','d^{sim}(q_c^{inv})','d^{meas}','q_c^{true}');
title("Initial Guess");
xlabel('$q_{c1n}$ Resistance','interpreter','LaTeX');
ylabel('Depth (m)','interpreter','LaTeX');
set(gca,'FontSize',12);
% subplot(1,2,2);
figure
t = LayerModelEval(layers,zd);
plot(t,-zd,'-r'), hold on, plot(layers(Nf+1:2*Nf),-layers(1:Nf),'or'), plot(mask*conv(t,kern),-zd,'--r');
plot(obs(:,2),-zd,'--b'), plot(qcTrue,-zd,'-b'), hold off;
legend('q_c^{inv}','Top of Layer','d^{sim}(q_c^{inv})','d^{meas}','q_c^{true}');
title("Standard Misfit");
xlabel('$q_{c1n}$ Resistance','interpreter','LaTeX');
ylabel('Depth (m)','interpreter','LaTeX');
set(gca,'FontSize',12);
