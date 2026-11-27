function plotPredictionRegion(x,y,z)

minx = min(x);
maxx = max(x);
miny = min(y);
maxy = max(y);

tx = linspace(minx,maxx)'; 
ty = linspace(miny,maxy); 

[XI,YI] = meshgrid(tx,ty);
ZI = griddata(x,y,z,XI,YI);
imagesc(x,y,ZI)