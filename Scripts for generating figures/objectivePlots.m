% script for produciing 3D objective function plots, comparing with log
% Jon Cooper

zd = 0:0.001:1;
trueLayers = [0;.5;.2;.8];
qcTrue = LayerModelEval(trueLayers,zd);
kern = chi2pdf(linspace(0,8,65),4);
kern = kern/sum(kern);
M = length(zd); mask = eye(M); n = 0; z = 0;
mask = [zeros(M,floor(length(kern)/2)),mask,zeros(M,ceil(length(kern)/2)-1)];
mask = diag([zeros(1,z),linspace(0,1,n),ones(1,M-n-z)].*[ones(1,M-n-z),linspace(1,0,n),zeros(1,z)])*mask;
blur = @(layer) mask*conv(layer,kern)';
qcMeas = blur(qcTrue);

scale = norm(ones(size(zd)),2);
misfit = @(l) norm(blur(LayerModelEval(l,zd))-qcMeas,2)/scale;
logMisfit = @(l) log(misfit(l));

layer = @(x) [0;x(1);.2;x(2)];
% layerPlot(zd,blur,trueLayers,qcMeas,"test");

%%
[X,Y] = meshgrid(0:.01:1,0:.01:1);
misfits = zeros(size(X));
for i = 1:size(X,1)
    for j = 1:size(X,2)
        misfits(i,j) = misfit(layer([X(i,j),Y(i,j)]));
    end
end

logMisfits = zeros(size(X));
for i = 1:size(X,1)
    for j = 1:size(X,2)
        logMisfits(i,j) = logMisfit(layer([X(i,j),Y(i,j)]));
    end
end

%%

figure
surfc(Y,-X,misfits), hold on
xlabel('Layer Resistance'); ylabel('Layer Location'); zlabel('Misfit Score');
plot3(qcTrue,-zd,zeros(size(zd))); shading interp

figure
surfc(Y,-X,logMisfits), hold on
xlabel('Layer Resistance'); ylabel('Layer Location'); zlabel('Misfit Score');
plot3(qcTrue,-zd,-5*ones(size(zd))); shading interp

% figure
% contour(Y,-X,log(misfits+1e-6)), hold on
% xlabel('Layer Resistance'); ylabel('Layer Location'); zlabel('Misfit Score');
% plot3(qcTrue,-zd,-5*ones(size(zd)));

%%


function layerPlot(zd,blur,model,meas,str)

figure
N = length(model)/2;
t = LayerModelEval(model,zd);
plot(t,-zd,'-r'), hold on, plot(model(N+1:2*N),-model(1:N),'or')
plot(blur(t),-zd,'--r'), plot(meas,-zd,'--b');
legend('$q_c$','Top of Layer','$\tilde{q_c}^{sim}$','$\tilde{q_c}^{meas}$','interpreter','latex');
title(str+", "+N+" Layers");
xlabel('$q_c$ Resistance (scaled)','interpreter','LaTeX');
ylabel('Depth (scaled)','interpreter','LaTeX');
set(gca,'FontSize',15);

end
