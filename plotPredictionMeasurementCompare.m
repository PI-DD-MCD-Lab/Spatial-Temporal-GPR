function plotPredictionMeasurementCompare(Measurements,Predictions,Var,name)

% global name nonMeaVal TimeInterval beginTime_FixLoc;
% Measurements = Measurements(beginTime_FixLoc:end);
% Predictions = Predictions((beginTime_FixLoc-1)/TimeInterval+1:end);
% Var = Var((beginTime_FixLoc-1)/TimeInterval+1:end);

TimeInterval = 1;
T = length(Measurements);
idx = (1:1:T)'; 
idx1 = (1:TimeInterval:T)';
% temp = Measurements~=-1;
% Mea = NaN*ones(T,1);
% Mea(temp) = Measurements;

% temp = CI*Var;
% upper = Predictions + temp;
% lower = Predictions - temp;

figure;
plot(idx1,Predictions,'--'); 
hold on;
plot(idx(~isnan(Measurements)),Measurements(~isnan(Measurements)),'*');
% fill([Predictions fliplr(Predictions)], [fliplr(upper) fliplr(upper)],'g');
fill([idx1;flip(idx1)], [Predictions + 3*sqrt(abs(squeeze(Var))); ...
      flip(Predictions - 3*sqrt(abs(squeeze(Var))))], 'g', 'EdgeColor','g','FaceAlpha',0.2);

legend('Prediction','Measurements','Confidence Interval');
xlabel('Time Instants');
ylabel('Prediction');
title(name);
% XLabel = XLabelSet(TimeType,T);
% set(gca,'XTickLabel',XLabel);