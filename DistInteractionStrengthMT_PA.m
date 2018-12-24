function rint = DistInteractionStrengthMT_PA(Nc,Nm,ri0)

% Interaction matrix based on strength probability distribution A
% Nc: number of interacting species
% ri0: maximum interaction strength
% rpn: ratio of positive to negative interactions
rint = ri0*(2*rand(Nc,Nm)-1);

return;