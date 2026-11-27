function [MaxMin,estimation] = PredictionFit(testNumtimestd,testNumtimeL,testNumnoise,timerange,noiserange,spacestd,time,Params,Ytrain,space,Y,YtrainLength,type,Ytest,PredTime)

space.std = spacestd;
[logspace,logtime,lognoisevar,maxlog,GCVspace,GCVtime,GCVnoisevar,minGCV,SUREspace,SUREtime,SUREnoisevar,minSURE] = EstParamsWeather(space,time,Params.spaceLocs,Ytrain,testNumtimestd,testNumtimeL,testNumnoise,timerange,Params.samplingTime,Params.timeInstants,noiserange);
logspace.std = spacestd;
GCVspace.std = spacestd;
SUREspace.std = spacestd;

MaxMin = [maxlog,minGCV,minSURE];
log{1} = logspace;
log{2} = logtime;
log{3} = lognoisevar;
gcv{1} = GCVspace;
gcv{2} = GCVtime;
gcv{3} = GCVnoisevar;
sure{1} = SUREspace;
sure{2} = SUREtime;
sure{3} = SUREnoisevar;
estimation{1} = log;
estimation{2} = gcv;
estimation{3} = sure;

