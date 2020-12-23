function layers = AOILog(layers,Misfit,TOL,pTOL,tolType,UB,MaxAdded,options,parallel,p)
% AOILog (Add-One-In) repeatedly tries to add extra layers between the
% existing ones unless the decrease in misfit is only marginal (i.e. the
% tolerance is a lower bound). Each added layer IS optimized, and this is
% just a 2-variable optimization keeping all else fixed.
%
% Jon Cooper
%
% Input:
%   layers  - layer model vector
%   Misfit  - misfit function
%   TOL     - either a percent or a fixed value, depends on tolType
%   pTOL    - percent allowed change in misfit, if used
%   tolType - 'iterative', 'absolute', 'and', or 'or'; see LayerOptimizer
%               options for more details
%   UB      - upper bound vector used in PSO, just need first and last
%               elements
%   options - PSO options
%   parallel- boolean for whether to use parallel
%   p       - plots misfits if greater than 2

bottom = UB(1); maxRes = UB(end);
iter = 0; added = 0;

N = length(layers)/2; Nc = N;% number of current layers
while 1 % loop condition not actually necessary; breaks inside
    
    base = Misfit(layers);
    if strcmp(tolType,'iterative')
        TOL = base*(1-pTOL/100);
    elseif strcmp(tolType,'and')
        TOL = max([TOL,base*(1-pTOL/100)]);
    elseif strcmp(tolType,'or')
        TOL = min([TOL,base*(1-pTOL/100)]);
    end
    
    misfits = [];
    if (strcmp(tolType,'and') || strcmp(tolType,'absolute')) && base <= TOL
        if p > 2 && ~isempty(misfits) % plot final misfits
            figure, bar(misfits), hold on;
            yline(base,'--k'); yline(TOL,'--r');
            legend('Misfits','Starting Misfit','Tolerance','location','NorthWest');
            title('Final AOI Misfits');
            xlabel('Layer Number'); ylabel('Misfit After Adding a Layer Below');
            set(gca,'FontSize',15);
        end
        break
    end
    
    % compute add-1-in models
    layersTemp = zeros(2*(Nc+1),1);
    layersTemp(2:Nc+1) = layers(1:Nc);
    layersTemp(Nc+3:end) = layers(Nc+1:end);
    misfits = zeros(Nc,1);
    extraLayers = zeros(Nc,2);
    if parallel
        layer_bottom = layers(Nc);
        parfor i = 1:Nc
            if i == Nc % insert below the bottom layer
                if layer_bottom == bottom % if bottom layer is at the very bottom already (unlikely)
                    misfits(i) = base;
                    continue
                end
                extraLayer = [(layers(Nc)+bottom)/2,layers(end)];
                LB = [layers(Nc);0]; UB = [bottom;maxRes];
                extraLayern = [layers(i)+rand()*(bottom-layers(i)),rand()*maxRes]; % with noise

            else % insert between 2 layers
                extraLayer = [(layers(i)+layers(i+1))/2,(layers(Nc+i)+layers(Nc+i+1))/2];
                LB = [layers(i);0]; UB = [layers(i+1);maxRes];
                extraLayern = [layers(i)+rand()*(layers(i+1)-layers(i)),rand()*maxRes];
            end
            options_copy = optimoptions('particleswarm',options);
            options_copy.InitialSwarmMatrix = [extraLayer;extraLayern];
            extraLayer = particleswarm(@(l) Misfit(imputeLayer(l,layersTemp)),2,LB,UB,options_copy);
            extraLayers(i,:) = extraLayer(:);
            misfits(i) = Misfit(imputeLayer(extraLayer,layersTemp));
        end
    else
        for i = 1:Nc
            if i == Nc % insert below the bottom layer
                if layers(Nc) == bottom % if bottom layer is at the very bottom already (unlikely)
                    misfits(i) = base;
                    continue
                end
                extraLayer = [(layers(Nc)+bottom)/2,layers(end)];
                LB = [layers(Nc);0]; UB = [bottom;maxRes];
                extraLayern = [layers(i)+rand()*(bottom-layers(i)),rand()*maxRes]; % with noise

            else % insert between 2 layers
                extraLayer = [(layers(i)+layers(i+1))/2,(layers(Nc+i)+layers(Nc+i+1))/2];
                LB = [layers(i);0]; UB = [layers(i+1);maxRes];
                extraLayern = [layers(i)+rand()*(layers(i+1)-layers(i)),rand()*maxRes];
            end
            options.InitialSwarmMatrix = [extraLayer;extraLayern];
            extraLayer = particleswarm(@(l) Misfit(imputeLayer(l,layersTemp)),2,LB,UB,options);

            layersTemp = imputeLayer(extraLayer,layersTemp);
            extraLayers(i,1) = extraLayer(1);
            extraLayers(i,2) = extraLayer(2);
            misfits(i) = Misfit(layersTemp);
        end
    end
    if p > 2 && Nc == N % if first iteration, plot initial misfits
        figure, bar(misfits), hold on;
        yline(base,'--k'); yline(TOL,'--r');
        legend('Misfits','Starting Misfit','Tolerance','location','NorthWest');
        title('Initial AOI Misfits');
        xlabel('Layer Number'); ylabel('Misfit After Adding a Layer Below');
        set(gca,'FontSize',15);
    end
    
    if all(misfits > TOL) % if jump in misfit is too small, exit
        if p > 2% plot final misfits
            figure, bar(misfits), hold on;
            yline(base,'--k'); yline(TOL,'--r');
            legend('Misfits','Starting Misfit','Tolerance','location','NorthWest');
            title('Final AOI Misfits');
            xlabel('Layer Number'); ylabel('Misfit After Adding a Layer Below');
            set(gca,'FontSize',15);
        end
        break
    else
        [~,i] = min(misfits); % otherwise, keep the layer that decreased the misfit the most
        layers = imputeLayer(extraLayers(i,:),layersTemp);
        Nc = Nc+1; % number of current layers
        
        pos = layers(1:Nc); res = layers(Nc+1:end);
        [pos,order] = sort(pos); res = res(order);
        layers = [pos;res]; layers(1) = 0;
        
        added = added+1;
    end
    
    if added > MaxAdded
        break
    end
end
end

% helper function for AOI
function layers = imputeLayer(l,layers)
    Nc = length(layers)/2-1;
    layers(1) = l(1);
    layers(Nc+2) = l(2);
end