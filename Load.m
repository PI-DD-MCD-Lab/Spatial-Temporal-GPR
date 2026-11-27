switch 'weather' % 'synthetic' 'weather'
    case 'synthetic'
        
        M = NumLocation*100 %%%%%%%max 1000 for exp 500 for matern
        Params.numLocs = M;
        Params.spaceLocsIdx = (1:Params.numLocs)';
        Params.spaceLocs = (1:Params.numLocs)';
        
        Params.samplingTime = 0.6;
        % Params.samplingTime = 0.05;
        % Params.samplingTime = 0.02;
        % Params.samplingTime = 0.01;
        % Params.samplingTime = 0.005;
        % Params.samplingTime = 0.002;
        Params.startTime = 0;
        Params.endTime = 12;
        
        Params.noiseStd = 0.1;
        
        space.type = 'gaussian';
        space.scale = 20;
        space.std = 1;
        Params.lossrate = 0.1;
        
        type = 'periodic';
        switch type
            case 'exponential'
                time.type = type;
                time.scale = 20;
                time.std = 1;
                load('Expo_simulateddata.mat')
                F = F(1:M,:);
                Y = Y(1:M,:);
                noiseVar = noiseVar(1:M,:);
            case 'Matern'
                time.type = type;
                time.l = 1;
                time.type = type;
                time.l = time.l;
                load('Mt_simulateddata.mat')
                F = F(1:M,:);
                Y = Y(1:M,:);
                noiseVar = noiseVar(1:M,:);
            case 'periodic'
                time.type = type;
                time.scale = 50;
                time.frequency = 1/12;
                time.std = 2;
                load('Per_simulateddata.mat')
                F = F(1:M,:);
                Y = Y(1:M,:);
                noiseVar = noiseVar(1:M,:);
        end
        Params.spaceLocsMeasIdx = sort(datasample(Params.spaceLocsIdx, round(1*Params.numLocs), 'Replace', false));
        Params.spaceLocsMeas = Params.spaceLocs(Params.spaceLocsMeasIdx,:);
        Params.spaceLocsPredIdx = setdiff(Params.spaceLocsIdx, Params.spaceLocsMeasIdx);
        Params.spaceLocsPred = Params.spaceLocs(Params.spaceLocsPredIdx,:);
        
        
        Params.timeInstants = (Params.startTime : Params.samplingTime : Params.endTime)';
        
        %         [F, Y, noiseVar] = datagenerate(Params,time,space);
        %         load Syn.mat;
        
        trainingrate = 0.99;
        YtrainLength = round(trainingrate*(Params.endTime-Params.startTime)/Params.samplingTime+1);
        YtestLength = size(Y,2) - YtrainLength;
        Ytrain = Y(:,1:YtrainLength);
        Ytest = Y(:,YtrainLength+1:end);
        Ftest = F(:,YtrainLength+1:end);
        noiseVartrain = noiseVar(:,1:YtrainLength);
        noiseVartest = noiseVar(:,YtrainLength+1:end);
        PredTime = YtestLength;
        [Num_Locations, Num_Time] = size(Ytrain);
        NM = sum(sum(~isnan(Ytrain)));
        Params.timeInstantstr = Params.timeInstants(1:end-YtestLength);
        Params.traintime=Params.timeInstantstr;
        
        space.type = 'gaussian';
        time.type = type;
        switch type
            case 'exponential'
                gcvspacerange = [15,25;0.5,1.5];
                gcvtimerange = [15,25;0.5,1.5];
                testNum = 10;
                %                 [logspace,logtime,GCVspace,GCVtime,SUREspace,SUREtime] = EstParams(space,time,Params.spaceLocs,Ytrain,testNum,gcvspacerange,gcvtimerange,noiseVartrain,Params.noiseStd,Params.timeInstants,Params.samplingTime);
            case 'Matern'
                gcvspacerange = [0.5,10;0.5,10];
                gcvtimerange = [0.1,1];
                testNum = 3;
                %                 [space,time] = GCVfilterMatern(space,time,Ytrain,Params.noiseStd^2,spacerange,timerange,Params.samplingTime,testNum,Params.spaceLocs);
                %                 [space,time,otherspace,othertime,otherspace1,othertime1] = EstParamsMatern(space,time,Params.spaceLocs,Ytrain,testNum,spacerange,timerange,noiseVartrain,Params.noiseStd,Params.timeInstants(1:YtrainLength));
        end
        
        
    case 'weather'
%                 load('TruncatedTemperature.mat');
%                 load('TruncatedTemperatureS12.mat');
%                 load('Temperature_North.mat');

%         load('TemperatureNorth12.mat');
%         load('Temperature20120S12.mat')
%         load('SmallTemperature20120.mat')
%         load('colorado.mat');
%         load('colorado20120.mat');
%         New_All_observation = measurementsPpt;
%         Location_3D = stationsLocations;
%         trainingrate = 0.9;
        load('TAVG20')
        Location_3D = Location(1:10,:);
        New_All_observation = TAVG(1:10,1:3652);
%         trainingrate = 6940/7305;
        trainingrate = 1- 365/3652;



        [M,N] = size(New_All_observation)
        idx = 1:M;
        Y = New_All_observation;
        Params.spaceLocsMeas = Location_3D;
        Params.noiseStd = 1;
        Params.samplingTime = 1;
        Params.startTime = 1;
        
        Params.spaceLocs = Location_3D(idx,:);
        Params.numLocs = size(Params.spaceLocs,1);
        Params.spaceLocsIdx = (1:Params.numLocs)';
        clear New_All_observation ;
        space.type = 'gaussian';
        %         space.type = 'gaussian_A';
%         time.type = 'PeriodicP';
%         time.type = 'exponentialTP3';
%         time.type = 'periodic';
        space.scale = 1;
        %         time.type = 'NewPer';
        %space.alt = 5;
        time.frequency = 1/12;
        time.scale = 17.8125;
        time.std = 3.1113;
        Effnoisevar =1.0185;
        space.std = 0.0283;
%         time.scale = 1400.31; %PeriodicP
%         time.std = 328.9711;
%         Effnoisevar = 1.6863;
%         space.std =  0.0248;
%         time.scale = 56075.3709; %PeriodicD1
%         time.std = 1491.1276;
%         Effnoisevar = 1;
%         space.std =  0.8155;
%         time.scale = 12.7455;
%         time.std = 1.5;
        
%         Effspace.type = 'gaussian';
%         %         Efftime.type = 'NewPer';
%         %         type = 'periodic';
%         Efftime.type = 'periodic';
%         type = 'periodic';
%         Effspace.scale = 1;
%         Effspace.std = 0.0892;
%         Efftime.frequency = 1/12;
%         Efftime.scale = 180.7455;
%         Efftime.std = 4.4292;
%         
        
        noiseVar = (Params.noiseStd .* ~isnan(Y)).^2;
        noiseVar(isnan(Y)) = Inf;
        noiseVar(Y==0) = Params.noiseStd^2;
        
        Params.spaceLocsMeasIdx = Params.spaceLocsIdx;
        Params.spaceLocsMeas = Params.spaceLocs(Params.spaceLocsMeasIdx,:);
        
        

        YtrainLength = floor(N*trainingrate); %round(trainingrate*(Params.endTime-Params.startTime)/Params.samplingTime+1);
        %         YtrainLength = round(trainingrate*newlength);
        YtestLength = N - YtrainLength;
        Ytrain = Y(:,1:YtrainLength);
        Ytest = Y(:,YtrainLength+1:end);
        Params.endTime = N;
        Params.timeInstants = (Params.startTime : Params.samplingTime : Params.endTime)';
        Params.traintime = (Params.startTime : Params.samplingTime :floor(trainingrate*Params.endTime))';
        Params.timeInstantstr =Params.traintime;
        PredTime = YtestLength;
        [Num_Locations, Num_Time] = size(Ytrain);
        NM = YtrainLength*M;
        %NM = sum(sum(~isnan(Ytrain)));
        
        %         noiseVar = (0.05 .* abs(Y)).^2;
        %         noiseVar(isnan(Y)) = Inf;
        %         noiseVar(Y==0) = 0.0025;
        noiseVar = Effnoisevar*noiseVar;
        noiseVartrain = noiseVar(:,1:YtrainLength);
        noiseVartest = noiseVar(:,YtrainLength+1:end);
        
end
