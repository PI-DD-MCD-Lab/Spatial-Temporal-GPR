function ML = logMarginal(Y,space,time,noisevar,NM,Params)
Locations = Params.spaceLocs;
% time_Covxx = CovMatrix(timevectortest,timevectortest,time);
% time_Covxy = CovMatrix(timevectortest,timevectortrain,time);
time_Covyy = CovMatrix(Params.traintime,Params.traintime,time);
spatial_Cov = CovMatrix(Locations,Locations,space);
N = size(Params.traintime,1);
Y = Y(:,1:N);
Y = Y(:);
I = eye(NM);
SigmaN = kron(time_Covyy,spatial_Cov) + noisevar * I;
[Q, R] = qr(SigmaN);
det_S = sum(log(abs(diag(R))));
inverse = (I/R)/Q;

% for i = 1:N
%     YinvY =YinvY+ Y'*inverse*Y;
% end
% ML = -0.5*(NM*log(2*pi) +det_S  + YinvY );

ML = -0.5*(NM*log(2*pi) +det_S  + Y'*inverse*Y );