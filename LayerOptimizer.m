function [layers,info] = LayerOptimizer(obs,blur,options)
% LayerOptimizer uses the Particle Swarm optimizer to fit a soil profile
% model to observed data via the pointspread function kern. Requires the
% image processing toolbox and the global optimization toolbox. Note that
% "layers" in this program is equivalent to "m" in Cooper et al. (2021).
% 
% Jon Cooper
%
% Input:
%   obs     - A Mx2 matrix with M soil resistance observations. The first
%               column is the (positive) depth vector, the second contains 
%               the observed resistances.
%   blur    - A function handle that takes a true resistance profile and
%               outputs a blurred resistance profile of the same
%               dimensions.
% Options : To be specified in name-value cell array,
%               e.g. options = {'AOITolType','or';'LOOpTol',5;'AOIpTol',5};
%   layers0 - A 2N column vector defining a soil layer model with N layers.
%               The first N elements are the layer positions, given in
%               absolute positive depth, the second N elements are the
%               respective layer resistances.
%   HybridFcn - Used in PSO options. ({'fmincon'}, [], PSO HybridFcn options)
%   LOOTolType : ('iterative','total','absolute','and',{'or'})
%           iterative - Exit LOO if misfit increases by more than LOOpTol
%                           percent compared to the previous iteration.
%           total - Exit LOO if misfit increases by more than LOOpTol
%                       percent compared to the initial misfit.
%           absolute - Exit LOO if misfit score is greater than LOOTol.
%           and - Exit LOO if the conditions for both 'total' and
%                   'absolute' are met.
%           or - Exit LOO if either of the conditions for 'total' and
%                   'absolute' are met.
%   LOOpTol - Percent allowed change in misfit for LOO, if used. {5}
%   LOOTol - Maximum absolute misfit score. {1}
%   AOITolType : ('iterative','absolute',{'and'},'or')
%           iterative - Only accept adding a new layer if it decreases the
%                           misfit by more than AOIpTol percent, compared
%                           to the previous iteration.
%           absolute - Only accept adding a new layer if the misfit is
%                           below AOITol.
%           and - Only accept a new layer if the conditions for both
%                   'iterative' and 'absolute' are met.
%           or - Only accept a new layer if either of the the conditions for
%                   'iterative' and 'absolute' are met.
%   AOIpTol - Percent allowed change in misfit for AOI, if used. Should be
%               strictly greater than LOOpTol, if both are used. {10}
%   AOITol - Minimum absolute misfit score. Should be strictly less than
%               LOOTol, if both are used. {1}
%   MisfitNorm - The norm type to be used in calculating the misfit.
%                   ({2}, vector norm options)
%   ConvTol - The maximum allowed change in the layer model to stop the
%               main algorithm loop (in the inf-sense). {1e-3}
%   MIN - Smallest possible layer resistance. {0.01}
%   MAX - Largest possible layer resistance. {10}
%   Plt - Integer for controlling what is plotted. 0 for none, 1 for initial,
%           coordinate descent, and final, 2 for LOO, AOI, and particleswarm
%           outputs within each iteration, 3 for all (including LOO & AOI misfits). {0}
%   MaxIter - Maximum number of iterations for the main loop. {10}
%   CD - Boolean for using coordinate descent on the initial layer model. {1}
%   CDMaxIter - Maximum number iterations for coordinate descent. {20}
%   SocialWeight - SocialAdjustmentWeight for PSO. {.7}
%   SelfWeight - SelfAdjustmentWeight for PSO. {2}
%   MisfitTransform - Function for transforming the misfit scores. {@(x) x}
% Output:
%   layers - a soil layer model (see layers0 description)
%   info - a struct containing most internal parameters

%% Define default options & read user-defined options
layers0 = [];
HybridFcn = 'fmincon';
LOOTolType = 'iterative';
LOOpTol = 10;
LOOTol = 0.02;
AOITolType = 'iterative';
AOIpTol = 20;
AOITol = 0.01;
MisfitNorm = 2;
ConvTol = 1e-3;
MIN = .01;
MAX = 10;
Plt = 0;
MaxIter = 10;
CD = 1;
CDMaxIter = 20;
SocialWeight = .7;
SelfWeight = 2;
MisfitTransform = @(x) x;

% rewrite default parameters if needed
if nargin == nargin(mfilename)
  for j = 1:size(options,1), eval([options{j,1},'= options{j,2};']); end
end

% adjust the tolerances if using an objective function transform
LOOTol = MisfitTransform(LOOTol);
AOITol = MisfitTransform(AOITol);

if ~isempty(layers0), N0 = length(layers0)/2; end

%% Unpack variables and initialize layer model

% transform data into [0,1]^2 range
zd = obs(:,1); qcMeas = obs(:,2);
mzd = max(zd); mqc = MAX;
zd = zd/mzd; qcMeas = qcMeas/mqc;

% define the misfit function
scale = norm(ones(size(zd)),MisfitNorm); % add weighting to ensure between 0 and 1
Misfit=@(l) MisfitTransform(norm(blur(LayerModelEval(l,zd))-qcMeas,MisfitNorm)/scale);

% define the initial layer model, if not given, and set properties
if isempty(layers0)
    pos = (qcMeas(3:end)-qcMeas(2:end-1) >= 0) & (qcMeas(1:end-2)-qcMeas(2:end-1) > 0);
    neg = (qcMeas(3:end)-qcMeas(2:end-1) <= 0) & (qcMeas(1:end-2)-qcMeas(2:end-1) < 0);
    loc = nonzeros((pos | neg).*zd(2:end-1)); res = nonzeros((pos | neg).*qcMeas(2:end-1));
    loc = [0;loc]; res = [qcMeas(1);res];
    layers0 = [loc;res]; clear loc
else
    pos = layers0(1:N0); res = layers0(N0+1:end);
    [pos,order] = sort(pos); res = res(order);
    layers0 = [min(zd)+pos/mzd;res/mqc]; layers0(1) = 0;
end
N0 = length(layers0)/2; N = N0;
UB = ones(2*N,1);
if Plt > 0, layerPlot(zd,psf,mask,layers0,qcMeas,"Initial"); end

% coordinate descent for initial model
if CD
    CDoptions = optimoptions('fminunc','Display','off','Algorithm','quasi-newton',...
        'HessUpdate','steepdesc','MaxFunctionEvaluations',100); % backdoor gradient descent
    for j = 1:CDMaxIter
        layers = layers0;
        for i = 2*N:-1:2 % go through resistances first, then locations, skipping location 0
            layers0(i) = fminunc(@(l) Misfit([layers0(1:(i-1));l;layers0((i+1):(2*N))]),layers0(i),CDoptions);
            if layers0(i) > 1 % enforce bounds on both resistances and locations (both scaled to [0,1])
                layers0(i) = 1;
            elseif layers0(i) < 0
                layers0(i) = 0;
            end
        end
        if norm(layers-layers0) < 1e-2, break, end
    end
    
    % error checking / correcting
    pos = layers0(1:N); res = layers0(N+1:end);
    pos(pos > max(zd)) = NaN;
    res = res(~isnan(pos)); pos = pos(~isnan(pos));
    res(res < MIN) = MIN; res(res > MAX) = MAX;
    [pos,order] = sort(pos); res = res(order);
    layers0 = [pos;res];

    if Plt > 0, layerPlot(zd,psf,mask,layers0,qcMeas,"Coordinate Descent"); end
end


%% Main Loop

% setup PSO
nPop = 100;
PSOoptions = optimoptions('particleswarm','SwarmSize',nPop,'HybridFcn',HybridFcn,...
    'Display','off','SocialAdjustmentWeight',SocialWeight,'SelfAdjustmentWeight',SelfWeight);

% initial LOO application to help speed up main loop
layers = layers0;
layers = LOO(layers0,Misfit,LOOTol,LOOpTol,LOOTolType,Plt); Nc = length(layers)/2;
pos = layers(1:Nc); res = layers(Nc+1:end);
[pos,order] = sort(pos); res = res(order);
layers = [pos;res]; layers(1) = 0;
if Plt > 1, layerPlot(zd,psf,mask,layers,qcMeas,"LOO Initial"); end

% MAIN LOOP
iter = 0;
layersPrev = layers + 1;
while iter < MaxIter && (length(layersPrev) ~= length(layers) || norm(layersPrev-layers,inf) > ConvTol)
    layersPrev = layers;
    
    if iter > 0
        % LOO
        layers = LOO(layers,Misfit,LOOTol,LOOpTol,LOOTolType,Plt); Nc = length(layers)/2;
        pos = layers(1:Nc); res = layers(Nc+1:end);
        [pos,order] = sort(pos); res = res(order);
        layers = [pos;res]; layers(1) = 0;
        if Plt > 1, layerPlot(zd,psf,mask,layers,qcMeas,"LOO"); end
        
        %PSO applied to entire model
        PSOoptions.InitialSwarmMatrix = [layers';layers'];
        pts = lhsdesign(PSOoptions.SwarmSize-2,2*Nc); % latin hypercube initial particle selection
        pts(:,1:Nc) = sortrows(pts(:,1:Nc)')';
        PSOoptions.InitialSwarmMatrix = [PSOoptions.InitialSwarmMatrix; pts];
        LB = [zeros(Nc,1)+min(zd); zeros(Nc,1)+MIN];
        UB = ones(2*Nc,1);
        layers = particleswarm(Misfit,2*Nc,LB,UB,PSOoptions)';
        pos = layers(1:Nc); res = layers(Nc+1:end);
        [pos,order] = sort(pos); res = res(order);
        layers = [pos;res]; layers(1) = 0;
        if Plt > 1, layerPlot(zd,psf,mask,layers,qcMeas,"PSO"); end
    end
    
    % AOI
    layers = AOI(layers,Misfit,AOITol,AOIpTol,AOITolType,UB,PSOoptions,Plt);
    if Plt > 1, layerPlot(zd,psf,mask,layers,qcMeas,"AOI"); end
    Nc = length(layers)/2;
    
    % PSO again
    PSOoptions.InitialSwarmMatrix = [layers';layers'];
    PSOoptions.OutputFcn = @saveParticles;
    LB = [zeros(Nc,1)+min(zd); zeros(Nc,1)+MIN];
    UB = ones(2*Nc,1);
    layers = particleswarm(Misfit,2*Nc,LB,UB,PSOoptions)';
    PSOoptions.OutputFcn = [];
    pos = layers(1:Nc); res = layers(Nc+1:end);
    [pos,order] = sort(pos); res = res(order);
    layers = [pos;res]; layers(1) = 0;
    if Plt > 1, layerPlot(zd,psf,mask,layers,qcMeas,"PSO"); end
    
    iter = iter+1;
end
% extra LOO application to clean up extra layers from last AOI
layers = LOO(layers,Misfit,LOOTol,LOOpTol,LOOTolType,Plt); Nc = length(layers)/2;
pos = layers(1:Nc); res = layers(Nc+1:end);
[pos,order] = sort(pos); res = res(order);
layers = [pos;res]; layers(1) = 0;
if Plt > 1, layerPlot(zd,psf,mask,layers,qcMeas,"LOO Final"); end

%% Warnings
errorfunc = blur(LayerModelEval(layers,zd))-qcMeas;
error = norm(errorfunc,MisfitNorm)/scale;
if error > .001 % ARBITRARY FOR NOW
    if error > 0.02, warning("High misfit: layer(s) likely averaged"); end
    
    % moving average, local L2 approximations
    zs = 1:length(zd);
    posi = (errorfunc(2:end) < 0) & (errorfunc(1:end-1) >= 0);
    neg = (errorfunc(2:end) > 0) & (errorfunc(1:end-1) <= 0);
    zs = zs(posi | neg); zs = [1;zs';length(zd)];
    
    for i = 1:length(zs)-2
        start = zs(i); fin = zs(i+2); width = zd(fin)-zd(start);
        eval1 = abs(sum(errorfunc(start:fin)))/width;
        eval2 = sum(abs(errorfunc(start:fin)))/width;

        % warn if local abs misfit significantly greater than signed misfit
        if eval2/eval1 > 2 && eval2 > 10
            warning("Missing/incorrect layers likely between zd="+...
                string(zd(start)*mzd)+" and zd="+string(zd(fin)*mzd));
        end
    end
end


%% Post processing

% transform back from [0,1]^2 range
fmisfit = Misfit(layers); imisfit = Misfit(layers0);
if Plt > 0, layerPlot(zd,psf,mask,layers,qcMeas,"Final"); end
layers = [(pos-min(zd))*mzd;res*mqc]; layers(1) = 0;
pos = layers0(1:N); res = layers0(N+1:end);
layers0 = [(pos-min(zd))*mzd;res*mqc]; layers0(1) = 0;

% output struct for more information
info = struct('Ninitial',N0,'Nfinal',Nc,'layers0',layers0,'blur',blur,...
    'InitialMisfit',imisfit,'FinalMisfit',fmisfit);

end

% internal plotting function (on [0,1]^2 scale)
% zd - depth list
% blur - blurring function
% model - layer model
% meas - observed data (blurred)
% str - string for title of plot
function layerPlot(zd,psf,mask,model,meas,str)

    figure
    N = length(model)/2;
    t = LayerModelEval(model,zd);
    plot(t,-zd,'-r'), hold on, plot(model(N+1:2*N),-model(1:N),'or')
    plot(mask*conv(t,psf),-zd,'--r'), plot(meas,-zd,'--b');
    legend('$q_c$','Top of Layer','$\tilde{q_c}^{sim}$','$\tilde{q_c}^{meas}$','interpreter','latex');
    title(str+", "+N+" Layers");
    xlabel('$q_c$ Resistance (scaled)','interpreter','LaTeX');
    ylabel('Depth (scaled)','interpreter','LaTeX');
    set(gca,'FontSize',15);
end

% saveParticles is a bit of a hackey method for scraping the final state of
%   PSO, despite this not being readily accessible to the user. Saves the
%   final particles as the file FinalSwarm.mat, if it is used.
% NOTE: this function is meant to be passed to PSO as options.OutputFcn
% optimValues - particle states
% state - internal state of the optimizer
function stop = saveParticles(optimValues,state)
    if strcmp(state,'done')
        stop = true;
        swarm = optimValues.swarm;
        pos = swarm(:,1:size(swarm,2)/2); res = swarm(:,size(swarm,2)/2+1:end);
        for i = 1:size(pos,1)
            [tpos,order] = sort(pos(i,:)); res(i,:) = res(i,order); pos(i,:) = tpos;
        end
        swarm = [pos,res];
        swarmfvals = optimValues.swarmfvals;
        save('FinalSwarm','swarm','swarmfvals');
    else
        stop = false;
    end
end