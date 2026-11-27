function gcv = stgp_gcv(x,Y,Params,space,type,time,NM)

% timescale = x(1);
% timestd = x(2);
% noisevar = x(3);
% spacestd = x(4);
time.scale = x(1);
time.std = x(2);
noisevar = x(3);
space.std = x(4);
time.c_t = x(5);
time.type = type;

gcv = gcvfunction(Y,space,time,noisevar,NM,Params)
