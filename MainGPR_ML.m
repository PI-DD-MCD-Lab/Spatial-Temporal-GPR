close all
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testNum_space_lengthscale = 1;
testNum_time_lengthscale =1;
testNum_time_Amp = 1;
testNum_time_c = 1;
testNum_noise = 1;

time_Amp_Range = linspace(1,20000,testNum_time_Amp); %
time_lengthscale_Range = linspace(1,1000,testNum_time_lengthscale);
timec = linspace(0.01,100,testNum_time_c); %
noisevar_Range = linspace(1,2,testNum_noise); % sigma^2
space_lengthscale_Range = linspace(1,10,testNum_space_lengthscale); %sigma_s

for j = 1:testNum_time_Amp
    for k = 1:testNum_time_lengthscale
        for l = 1:testNum_noise
            for h = 1:testNum_time_c
            hyper = [time_Amp_Range(j),time_lengthscale_Range(k),noisevar_Range(l) timec(h)];
            Hyper{j,k,l,h} = hyper;
            end
        end
    end
end

Hyper = Hyper(:);
Load;
time.type = 'LocallyPer' ;%'exponentialTP4'
type = 'LocallyPer';
time.frequency=1/365;
tic;
for j = 1:testNum_space_lengthscale
    %for h =1:testNumspaceL
    space_lengthscale = space_lengthscale_Range(j);
    space.std = space_lengthscale;

    for i= 1:(testNum_time_lengthscale*testNum_time_Amp*testNum_noise*testNum_time_c)
        hyper = Hyper{i};
        time_Amp = hyper(1);
        time_lengthscale = hyper(2);
        noisevar = hyper(3);
        time.std =  time_lengthscale;
        time.scale = time_Amp;        
        time.c_t=hyper(4);
        timect =time.c_t;
        lb = eps*ones(1,5);
        ub = inf*ones(1,5);
        ub(1,5)=1; % when c = 1/time.c_t;
%         lb(1,5)=1; % when c = 1/time.c_t;
        ub(1,2) = 20;
        x0 = [time_Amp time_lengthscale noisevar space_lengthscale timect];
        
        ML = logMarginal(Y,space,time,noisevar,NM,Params)

        ff1 = @(x)stgp_loglklhd(x,Y,Params,space,type,time,NM);
        [x,fval,exitflag,output] = fmincon(ff1,x0,[],[],[],[],lb,ub,[]);
        
        NLog{j,i} =fval; 
        HyperML{j,i} = x;
    end
end
toc

Nlog = cell2mat(NLog);
[MAX_ML,I_ML] =min(Nlog(:));
[I_row1, I_column1] = ind2sub(size(Nlog),I_ML);

disp('optimal ML function value');
MAX_ML
disp('optimal ML hyper-parameter');
NLog{I_row1,I_column1}

time.scale = HyperML{I_row1,I_column1}(1);
time.std = HyperML{I_row1,I_column1}(2);
noisevar = HyperML{I_row1,I_column1}(3);
space.std = HyperML{I_row1,I_column1}(4);
time.c_t=HyperML{I_row1,I_column1}(5);
logspace = space;
logtime = time;
logtime.scale =HyperML{I_row1,I_column1}(1);
logtime.std = HyperML{I_row1,I_column1}(2);
lognoisevar = HyperML{I_row1,I_column1}(3);
logspace.std = HyperML{I_row1,I_column1}(4);
logtime.c_t = HyperML{I_row1,I_column1}(5);

noiseVar = noisevar*ones(M,N);
noiseVartrain = noiseVar(:,1:YtrainLength);
noiseVartest = noiseVar(:,YtrainLength+1:end);
%%%Hyperparameters for test
% load('SpaceTimeN.mat');
% logspace = space;
% logtime = time;

lognoiseVar = lognoisevar*ones(M,N);
lognoiseVartrain = lognoiseVar(:,1:YtrainLength);
[logPredGP_I,logVarGP_I] = GPR(Params.spaceLocsMeas,Ytrain,logspace,...
    logtime,Params.timeInstants(YtrainLength+1:end),Params.timeInstants(1:YtrainLength),lognoiseVartrain);

idxtest = ~isnan(Ytest);

Nfittestlog = 100*(1-norm(Ytest(idxtest)-logPredGP_I(idxtest))/norm(Ytest(idxtest)));
Nfittestlogvector = zeros(PredTime,1);
for i = 1:PredTime
    Nfittestlogvector(i) = 100*(1-norm(Ytest(idxtest(:,i),i)-logPredGP_I(idxtest(:,i),i))/norm(Ytest(idxtest(:,i),i)));
end
locationIdx = 18;
locationIdxstr = num2str(locationIdx);
MLPTitle = ['ML Prediction at Location ',locationIdxstr];
plotPredictionMeasurementCompare(Ytest(locationIdx,:)',logPredGP_I(locationIdx,:)',reshape(logVarGP_I(locationIdx,locationIdx,:),size(logPredGP_I(locationIdx,:)))',MLPTitle)

a = 1:PredTime;
figure
plot(a,Nfittestlogvector);
title('Prediction fit')


