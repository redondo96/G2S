function errs = err_samples_prediction_stdev(gprMdl,samples_comp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ypred = predict(gprMdl,samples_comp(:,1:end-1));
% [ypred,ysd,~] = predict(gprMdl,samples_comp(:,1:end-1));

diffs = abs(samples_comp(:,end) - ypred);  % ./samples_comp(:,end)

errs = [diffs, zeros(size(diffs))];  % ysd

end

