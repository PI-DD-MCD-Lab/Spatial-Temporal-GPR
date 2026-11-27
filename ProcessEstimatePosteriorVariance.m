function ProcessEstimatePosteriorVariance(TestLocations,TrainingLocations,Measurements,PreLoc,hat_f,Var_f)

Num_TrainingLoc = length(TrainingLocations);
I_Loc = eye(Num_TrainingLoc);
MeaIdx = (Measurements~=-1);
Ik = I_Loc(MeaIdx,:);
y = Ik*Measurements;
TrainingLoc = Ik*TrainingLocations;

figure;
plotPredictionRegion(PreLoc(:,1),PreLoc(:,2),hat_f); hold
plot(TrainingLoc(:,1),TrainingLoc(:,2),'k*');
plot(TestLocations(:,1),TestLocations(:,2),'bo');
xlabel( 'longitude' );
ylabel( 'latitude' );
zlabel( 'Predictions' );
grid on
colorbar

figure;
plotPredictionRegion(PreLoc(:,1),PreLoc(:,2),Var_f);hold
plot(TrainingLoc(:,1),TrainingLoc(:,2),'k*');
plot(TestLocations(:,1),TestLocations(:,2),'bo');
xlabel( 'longitude' );
ylabel( 'latitude' );
zlabel( 'Posterior Variance' );
grid on
colorbar