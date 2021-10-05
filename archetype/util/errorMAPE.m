function err = errorMAPE(prediction,truth)
%% errorMAPE - Mean Absolute Percentage Error - cannot handle zero values
% err = errorMAPE(prediction,truth)

err = mean(abs(prediction-truth)./truth).*100;

end