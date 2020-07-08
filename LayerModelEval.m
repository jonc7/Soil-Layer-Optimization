function eval = LayerModelEval(layers,zz)
% LayerModelEval is able to evaluate a soil layer model at all points zz.
% This function uses piecewise constant interpolation and returns a column
% vector of the evaluation points. If a point outside the range of the
% layer model is given, the resistance of the nearest layer is returned.
%
% Jon Cooper
%
% Input:
%   layers - A 2N column vector defining a soil layer model with N layers.
%            The first N elements are the layer positions, given in
%            absolute positive depth, the second N elements are the
%            respective layer resistances.
%   zz - A vector of positive depths to evaluate the model at.

N = length(layers)/2;
pos = layers(1:N); res = layers(N+1:end);
[pos,order] = sort(pos); res = res(order); % must be sorted here for mid-optimization evals
eval = zeros(size(zz));
for i = 1:length(eval)
    eval(i) = max([sum(pos <= zz(i)),1]);
    eval(i) = res(eval(i));
end

end