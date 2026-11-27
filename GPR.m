function [hat_F, Sigma_F] = GPR(Locations,Measurements,...
    spatial,time,timevectortest,timevectortrain,noiseVar)

[Num_Locations, ~] = size(Measurements);
Num_Time = length(timevectortest);

time_Covxx = CovMatrix(timevectortest,timevectortest,time);
time_Covxy = CovMatrix(timevectortest,timevectortrain,time);
time_Covyy = CovMatrix(timevectortrain,timevectortrain,time);
% time_Covxx = time.scale*ones(size(timevectortest,1));
% time_Covxy = time.scale*ones(size(timevectortest,1),size(timevectortrain,1));
% time_Covyy = time.scale*ones(size(timevectortrain,1),size(timevectortrain,1));

spactial_Cov = CovMatrix(Locations,Locations,spatial);

Kyy = kron(time_Covyy,spactial_Cov);
Kxy = kron(time_Covxy,spactial_Cov);
Kxx = kron(time_Covxx,spactial_Cov);
Measurements = Measurements(:);
noiseVar = sparse(diag(noiseVar(:)));

MeaIdx = ~isnan(Measurements);
I = eye(sum(MeaIdx));
temp = Kyy(MeaIdx,MeaIdx) + noiseVar(MeaIdx,MeaIdx);
[Q, R] = qr(temp);
inverse = (I/R)/Q;
F = Kxy(:,MeaIdx)*(inverse*Measurements(MeaIdx));
% NF = time_Covxy*inv(time_Covyy+ 3.3749*eye(108,108))*(Measurements(1,:))';
Sigma = Kxx - Kxy(:,MeaIdx)*inverse*Kxy(:,MeaIdx)';

hat_F = reshape(F,Num_Locations, Num_Time);
Sigma_F = zeros(Num_Locations,Num_Locations, Num_Time);
for t = 1:Num_Time
    idx = ((1:Num_Locations) + (t-1)*Num_Locations)';
    Sigma_F(:,:,t) = Sigma(idx,idx);
end
end