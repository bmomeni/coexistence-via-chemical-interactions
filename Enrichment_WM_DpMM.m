%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
% ExMT2: corrected the error in ExMT, now rIntMat only includes links in R

clear
rndseed0 = 8712;
rand('twister',rndseed0);

nSample = 5000; % # of samples being screened
rndseed = round(100*nSample*rand(1,nSample));

nGen = 200;
nInitialCell = 1e4; % total initial cells
dilTh = 1e10; % coculture dilution threshold
GenPerRound = log(dilTh/nInitialCell)/log(2);
nRound = round(nGen/GenPerRound); % number of rounds of propagation

nCellType = 20; % # of cell types in the initial pool
nMediator = 15; % # of mediators
kSatLevel = 1e4; % interaction strength saturation level of each population
extTh = 0.1; % population extinction threshold
ri0 = 0.2; % maximum interaction strength, 1/hr
posIntRatio = 0.1; % fraction of interactions that are positive
tau0 = 0; % in hours
tauf = 250; % in hours
dtau = 0.01; % in hours, cell growth update and uptake timescale
del = 0; % avg. decay rate (per hour) 
at = 1; % avg. consumption values (fmole per cell); alpha_ij: population i, resource j
bt = 0.1; % avg. production rates (fmole per cell per hour); beta_ij: population i, resource j
qp = 0.7; % probability of production link per population
qc = 0.7; % probability of influence link per population

r0T = zeros(nCellType,nSample); %matrix of zeros: 4 rows because 4 species and 1000 columns because 1000 samples being screened
SiT = zeros(nMediator,nSample); 
AT = zeros(nMediator,nCellType,nSample);
BT = zeros(nMediator,nCellType,nSample);
DT = zeros(nMediator,nSample);
CmpADT = zeros(nCellType,nSample);
CmpBDT = zeros(nCellType,nSample);
CmpEDT = zeros(nCellType,nSample);
CmpART = zeros(nCellType,nSample);
CmpBRT = zeros(nCellType,nSample);
CmpERT = zeros(nCellType,nSample);
rintAT = zeros(nMediator,nCellType,nSample);
rintBT = zeros(nMediator,nCellType,nSample);
rintET = zeros(nMediator,nCellType,nSample);
V0DT = zeros(3,nCellType,nSample);
VDT = zeros(3,nCellType,nSample);
V0RT = zeros(3,nCellType,nSample);
VRT = zeros(3,nCellType,nSample);

DAD = zeros(1,nSample); % d is depletable, r is reusable, ABE are distrubutions
DBD = zeros(1,nSample);
DED = zeros(1,nSample);
DAR = zeros(1,nSample);
DBR = zeros(1,nSample);
DER = zeros(1,nSample);

NE0D = zeros(3,nSample);
NED = zeros(3,nSample);
NE0R = zeros(3,nSample);
NER = zeros(3,nSample);

for ns = 1:nSample,
    disp(ns)
    tic
    
    rand('twister',rndseed(ns));
    r0 = 0.08+0.04*rand(nCellType,1); % population reproduction rates, per hour
    kSatVector = kSatLevel * (0.5 + rand(nMediator, 1)); % population levels for influence saturation

    %% Parameters
    % Network configuration
    % NetworkConfig_Balanced(Nc,Nm,q): link between Nc and Nm present with
    % a probability q
    R = NetworkConfig_Binomial(nCellType,nMediator,qc);
    P = NetworkConfig_Binomial(nCellType,nMediator,qp);

    % interaction matrix
    alpha = at * (0.5+rand(nCellType,nMediator)); % consumption rates
    beta = bt * (0.5+rand(nCellType,nMediator)); % mediator release rates
    delta = del * (0.5+rand(nMediator,1)); % decay rates, reusable chemicals only
    A = (R.*alpha)';
    B = (P.*beta)';
    D = delta;

    rIntMatA = R .* DistInteractionStrengthMT_PA(nCellType, nMediator, ri0); % matrix of interaction coefficients, 50/50
    rIntMatB = R .* DistInteractionStrengthMT_PB(nCellType, nMediator, ri0, posIntRatio); % matrix of interaction coefficients, more negative
    rIntMatE = R .* DistInteractionStrengthMT_PB(nCellType, nMediator, ri0, 1-posIntRatio); % matrix of interaction coefficients, more positive

    cellRatioArray = 1 / nCellType * ones(1,nCellType); % cell distribution; population ratios
    
    %% Simulating dynamics, Dp, depletable
    [NeAD, CmpAD] = WellmixedInteraction_DpMM_ExMT4(nRound,r0,cellRatioArray,rIntMatA,nInitialCell,kSatVector,A,B,kSatLevel,extTh,dilTh,tauf,dtau);
    [NeBD, CmpBD] = WellmixedInteraction_DpMM_ExMT4(nRound,r0,cellRatioArray,rIntMatB,nInitialCell,kSatVector,A,B,kSatLevel,extTh,dilTh,tauf,dtau);
    [NeED, CmpED] = WellmixedInteraction_DpMM_ExMT4(nRound,r0,cellRatioArray,rIntMatE,nInitialCell,kSatVector,A,B,kSatLevel,extTh,dilTh,tauf,dtau);

    %% Simulating dynamics, Ru, reusable
    [NeAR, CmpAR] = WellmixedInteraction_Ru_ExMT4(nRound,r0,cellRatioArray,rIntMatA,nInitialCell,kSatVector,D,B,extTh,dilTh,tauf,dtau);
    [NeBR, CmpBR] = WellmixedInteraction_Ru_ExMT4(nRound,r0,cellRatioArray,rIntMatB,nInitialCell,kSatVector,D,B,extTh,dilTh,tauf,dtau);
    [NeER, CmpER] = WellmixedInteraction_Ru_ExMT4(nRound,r0,cellRatioArray,rIntMatE,nInitialCell,kSatVector,D,B,extTh,dilTh,tauf,dtau);

    V0AD = zeros(1,nCellType);
    V0BD = zeros(1,nCellType);
    V0ED = zeros(1,nCellType);
    V0AD(NeAD) = 1;
    V0BD(NeBD) = 1;
    V0ED(NeED) = 1;

    V0AR = zeros(1,nCellType);
    V0BR = zeros(1,nCellType);
    V0ER = zeros(1,nCellType);
    V0AR(NeAR) = 1;
    V0BR(NeBR) = 1;
    V0ER(NeER) = 1;

    NE0D(:,ns) = sum([V0AD; V0BD; V0ED],2);
    NE0R(:,ns) = sum([V0AR; V0BR; V0ER],2);

    Cmp0AD = zeros(1,nCellType);
    Cmp0BD = zeros(1,nCellType);
    Cmp0ED = zeros(1,nCellType);
    Cmp0AD(NeAD) = CmpAD;
    Cmp0BD(NeBD) = CmpBD;
    Cmp0ED(NeED) = CmpED;
    
	Cmp0AR = zeros(1,nCellType);
    Cmp0BR = zeros(1,nCellType);
    Cmp0ER = zeros(1,nCellType);
    Cmp0AR(NeAR) = CmpAR;
    Cmp0BR(NeBR) = CmpBR;
    Cmp0ER(NeER) = CmpER;
	
	CmpADT(:,ns) = Cmp0AD;
	CmpBDT(:,ns) = Cmp0BD;
	CmpEDT(:,ns) = Cmp0ED;
	CmpART(:,ns) = Cmp0AR;
	CmpBRT(:,ns) = Cmp0BR;
	CmpERT(:,ns) = Cmp0ER;
	
    r0T(:,ns) = r0;
    SiT(:,ns) = kSatVector;
    AT(:,:,ns) = A;
    BT(:,:,ns) = B;
    DT(:,ns) = D;
    rintAT(:,:,ns) = rIntMatA';
    rintBT(:,:,ns) = rIntMatB';
    rintET(:,:,ns) = rIntMatE';
    V0DT(:,:,ns) = [V0AD; V0BD; V0ED];
    V0RT(:,:,ns) = [V0AR; V0BR; V0ER];
    toc
end

save(strcat('StabilityScreenCmp_Dp_vs_Ru_ExMT4_ABE_fp',num2str(round(100*posIntRatio)),'_CSD',num2str(nInitialCell,2),'_DilTh',num2str(dilTh,2),'_ExtTh',num2str(extTh),'_Ksat',num2str(kSatLevel),'_ri',num2str(round(100*ri0)),'_bt',num2str(round(100*bt)),'_at',num2str(round(100*at)),'_A_Nc',num2str(nCellType),'_Nm',num2str(nMediator),'_qp',num2str(round(100*qp)),'_qc',num2str(round(100*qc)),'_Nr',num2str(nRound),'_Ns',num2str(nSample),'_rndseed',num2str(rndseed0),'.mat'))

