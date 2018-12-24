function rint = DistInteractionStrengthMT_PB(Nc,Nm,ri0,fp)

% Interaction matrix based on strength probability distribution B
% Nc: number of interacting species
% ri0: maximum interaction strength
% fp: fraction for positive interactions
rint = ri0*rand(Nc,Nm).*sign(rand(Nc,Nm)-(1-fp));

return;