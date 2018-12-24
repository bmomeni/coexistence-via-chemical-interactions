%% This program shows the dynamics leading to final stable community during
%% enrichment for a specific sample among cases simulated in the screen

%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators

clear

inpath = 'Z:\Members\Ian\Experiments\ExMT4 PowerLaw vs Binomial\Data\Nc (old parameters)\PowerLaw\';
infile = 'StabilityScreenCmp_Dp_vs_Ru_ExMT4_ABE_fp20_CSD1e+04_DilTh1e+10_ExtTh0.1_Ksat10000000_ri10_bt100_at10_A_Nc25_Nm15_qp140_qc140_Nr15_Ns10000_rndseed6035.mat';

load(strcat(inpath,infile))

Ncxt = 2;
slc = 1; % select network: 1:A, 2:B, 3:E

Nc = nCellType;
Nm = nMediator;
Nr = nRound;
CSD = nInitialCell;

NstD = shiftdim(sum(V0DT(slc,:,:),2),1);
pck = find(NstD>=Ncxt);
Npck = length(pck);
VT = zeros(10*Npck,12);
CXT = zeros(10*Npck,Nc);
VTE = zeros(10*Npck,12);
CXTE = zeros(10*Npck,Nc);

cntinst = 0;
cntinstE = 0;
for cntp = 1:Npck
    disp([num2str(cntp),' out of ',num2str(Npck)])
    sc = pck(cntp); % which instance in the screen

    if slc == 1
        rint = rintAT(:,:,sc);
    end
    if slc == 2
        rint = rintBT(:,:,sc);
    end
    if slc == 3
        rint = rintET(:,:,sc);
    end
    A = AT(:,:,sc);
    B = BT(:,:,sc);
    r0 = r0T(:,sc);
    Si = SiT(:,sc);

    %% Dp enriched
    rith = ri0/3000;
    rintEn = zeros(Nm,Nc);
    AtEn = zeros(Nm,Nc);
    BtEn = zeros(Nm,Nc);
    Ncrng = 1:Nc;
    NeD = Ncrng(V0DT(slc,:,sc)==1);
    rintEn(:,NeD) = rint(:,NeD);
    rintEn = rintEn.*((sum(BT(:,NeD,sc),2)>0)*ones(1,Nc));
    rintEn(abs(rintEn)<rith) = 0; % remove weak interactions
    AtEn(:,NeD) = AT(:,NeD,sc);
    AtEn = AtEn.*((sum(BT(:,NeD,sc),2)>0)*ones(1,Nc));
    BtEn(:,NeD) = BT(:,NeD,sc);
    BtEn = BtEn.*((sum(AT(:,NeD,sc),2)>0)*ones(1,Nc));

    NeDr = [];
    cntst = 0;
    for scnt = 1:length(NeD)
        if sum(AtEn(:,NeD(scnt)))>0
            cntst = cntst+1;
            NeDr = [NeDr NeD(scnt)];
        end
    end

    [Nc, Nm] = size(rint');
    %% Initial state
    rp = zeros(1,Nc); % cell distrbution
    rp(NeDr) = 1/length(NeDr); % cell distrbution
    NDs = WellmixedInteraction_DpMM_ExMT4(Nr,r0,rp,rintEn',CSD,Si,AtEn,BtEn,kSatLevel,extTh,dilTh,tauf,dtau);
    rintEn = zeros(Nm,Nc);
    AtEn = zeros(Nm,Nc);
    BtEn = zeros(Nm,Nc);
    rintEn(:,NDs) = rint(:,NDs);
    rintEn = rintEn.*((sum(BT(:,NDs,sc),2)>0)*ones(1,Nc));
    % rintEn(abs(rintEn)<rith) = 0; % remove weak interactions
    AtEn(:,NDs) = AT(:,NDs,sc);
    AtEn = AtEn.*((sum(BT(:,NDs,sc),2)>0)*ones(1,Nc));
    BtEn(:,NDs) = BT(:,NDs,sc);
    BtEn = BtEn.*((sum(AT(:,NDs,sc),2)>0)*ones(1,Nc));
    Nchc = find(rintEn); % index of influences that can change
    V = zeros(length(Nchc)+1,Nc);
    V(length(Nchc)+1,NDs) = 1;
    VE = zeros(length(Nchc)+1,Nc);
    VE(length(Nchc)+1,NDs) = 1;
    Ninst = sum(V(length(Nchc)+1,:),2);

    if Ninst > 1
        % examine stable community
        for nst = 1:length(Nchc)
            %% Initial state
            rp = zeros(1,Nc); % cell distrbution
            rp(NeDr) = 1/length(NeDr); % cell distrbution
            rintEch = rintEn;
            rintEch(Nchc(nst)) = 0;
            NDRmv = WellmixedInteraction_DpMM_ExMT4(Nr,r0,rp,rintEch',CSD,Si,AtEn,BtEn,kSatLevel,extTh,dilTh,tauf,dtau);
            V(nst,NDRmv) = 1;

            cntinst = cntinst+1;
            disp(cntinst)
            VT(cntinst,1) = Ninst; % number of species in the examined instance
            VT(cntinst,2) = sum(abs(V(nst,:)-V(length(Nchc)+1,:)),2); % change in community structure
            VT(cntinst,3) = sign(rintEn(Nchc(nst)))*(1+(BtEn(Nchc(nst))>0)); % category
            VT(cntinst,4) = rintEn(Nchc(nst)); % interaction coefficient
            VT(cntinst,5) = sum(sum(rintEn(:,NDs).*AtEn(:,NDs)>0)); % total number of positive influence links
            VT(cntinst,6) = sum(sum(rintEn(:,NDs).*AtEn(:,NDs)<0)); % total number of negative influence links
            VT(cntinst,7) = sum(sum(rintEn(:,NDs).*AtEn(:,NDs).*BtEn(:,NDs)>0)); % number of positive influences on producers
            VT(cntinst,8) = sum(sum(rintEn(:,NDs).*AtEn(:,NDs).*BtEn(:,NDs)<0)); % number of positive influences on producers
            VT(cntinst,9) = sum(sum(rintEn(:,NDs).*AtEn(:,NDs).*(1-(BtEn(:,NDs)>0))>0)); % number of positive influences on nonproducers
            VT(cntinst,10) = sum(sum(rintEn(:,NDs).*AtEn(:,NDs).*(1-(BtEn(:,NDs)>0))<0)); % number of positive influences on nonproducers
            VT(cntinst,11) = sc; % index in the original screen
            VT(cntinst,12) = Nchc(nst); % index of removed link

            CXT(cntinst,NDRmv) = 1; % list of coexisting species
        end
        % examine enrichment
        for nst = 1:length(Nchc)
            %% Initial state
            rp = 1/Nc*ones(1,Nc); % cell distrbution
            rintch = rint;
            rintch(Nchc(nst)) = 0;
            NDRmv = WellmixedInteraction_DpMM_ExMT4(Nr,r0,rp,rintch',CSD,Si,A,B,kSatLevel,extTh,dilTh,tauf,dtau);
            VE(nst,NDRmv) = 1;

            cntinstE = cntinstE+1;
            disp(cntinstE)
            VTE(cntinstE,1) = Ninst; % number of species in the examined instance
            VTE(cntinstE,2) = sum(abs(VE(nst,:)-VE(length(Nchc)+1,:)),2);
            VTE(cntinstE,3) = sign(rint(Nchc(nst)))*(1+(B(Nchc(nst))>0)); % category
            VTE(cntinstE,4) = rint(Nchc(nst)); % interaction coefficient
            VTE(cntinstE,5) = sum(sum(rint.*A>0)); % total number of positive influence links
            VTE(cntinstE,6) = sum(sum(rint.*A<0)); % total number of negative influence links
            VTE(cntinstE,7) = sum(sum(rint.*A.*B>0)); % number of positive influences on producers
            VTE(cntinstE,8) = sum(sum(rint.*A.*B<0)); % number of positive influences on producers
            VTE(cntinstE,9) = sum(sum(rint.*A.*(1-(B>0))>0)); % number of positive influences on nonproducers
            VTE(cntinstE,10) = sum(sum(rint.*A.*(1-(B>0))<0)); % number of positive influences on nonproducers
            VTE(cntinstE,11) = sc; % index in the original screen
            VTE(cntinstE,12) = Nchc(nst); % index of removed link

            CXTE(cntinstE,NDRmv) = 1; % list of coexisting species
        end
        disp([cntinst VT(cntinst,2) cntinstE VTE(cntinstE,2)])
    end

end

IntType = [-2 -1 1 2];
winf = find(VT(1:cntinst,2)>0);
woinf = find(VT(1:cntinst,2)==0);
figure
bar(1/cntinst*[hist(VT(winf,3),IntType); hist(VT(woinf,3),IntType)]',0.3,'stacked')
set(gca,'XTickLabel',{'- on self','- on others','+ on others','+ on self'})
legend('Influencing coexistence','Not influencing coexistence')
xlabel('Chemical influence type')
ylabel('Frequency')
title('Removing influence links - Stable-State')

IntType = [-2 -1 1 2];
winf = find(VTE(1:cntinstE,2)<0);
woinf = find(VTE(1:cntinstE,2)==0);
figure
bar(1/cntinstE*[hist(VTE(winf,3),IntType); hist(VTE(woinf,3),IntType)]',0.3,'stacked')
set(gca,'XTickLabel',{'- on self','- on others','+ on others','+ on self'})
legend('Influencing coexistence','Not influencing coexistence')
xlabel('Chemical influence type')
ylabel('Frequency')
title('Removing influence links - Enrichment')

% positive influence on others
pioSS = find(VT(1:cntinst,3)==1);
pochSS(1) = 1/length(pioSS)*sum(VT(pioSS,2)>0); % fewer coexisting
pochSS(2) = 1/length(pioSS)*sum(VT(pioSS,2)==0); % same coexisting

% positive influence on self
pisSS = find(VT(1:cntinst,3)==2);
pschSS(1) = 1/length(pisSS)*sum(VT(pisSS,2)>0); % fewer coexisting
pschSS(2) = 1/length(pisSS)*sum(VT(pisSS,2)==0); % same coexisting

% negative influence on others
nioSS = find(VT(1:cntinst,3)==-1);
nochSS(1) = 1/length(nioSS)*sum(VT(nioSS,2)>0); % fewer coexisting
nochSS(2) = 1/length(nioSS)*sum(VT(nioSS,2)==0); % same coexisting

% negative influence on self
nisSS = find(VT(1:cntinst,3)==-2);
nschSS(1) = 1/length(nisSS)*sum(VT(nisSS,2)>0); % fewer coexisting
nschSS(2) = 1/length(nisSS)*sum(VT(nisSS,2)==0); % same coexisting

figure
bar([nschSS' nochSS' pochSS' pschSS']','stack')
set(gca,'XTickLabel',{'- on self','- on others','+ on others','+ on self'})
ylabel('Frequency')
xlabel('Category of influence removed')
text(0.7,1.1,num2str(length(nisSS)))
text(1.5,1.1,num2str(length(nioSS)))
text(2.7,1.1,num2str(length(pioSS)))
text(3.7,1.1,num2str(length(pisSS)))
text(0.5,0.8,'same')
text(0.5,0.2,'fewer')
ylim([0 1.2])

sch(1) = 1/cntinstE*sum(sum(CXTE(1:cntinstE,:),2)<VTE(1:cntinstE,1)); % fewer coexisting
sch(2) = 1/cntinstE*sum(sum(CXTE(1:cntinstE,:),2)==VTE(1:cntinstE,1)); % same coexisting
sch(3) = 1/cntinstE*sum(sum(CXTE(1:cntinstE,:),2)>VTE(1:cntinstE,1)); % more coexisting

pio = find(VTE(1:cntinstE,3)==1);
poch(1) = 1/length(pio)*sum(sum(CXTE(pio,:),2)<VTE(pio,1)); % fewer coexisting
poch(2) = 1/length(pio)*sum(sum(CXTE(pio,:),2)==VTE(pio,1)); % same coexisting
poch(3) = 1/length(pio)*sum(sum(CXTE(pio,:),2)>VTE(pio,1)); % more coexisting

pis = find(VTE(1:cntinstE,3)==2);
psch(1) = 1/length(pis)*sum(sum(CXTE(pis,:),2)<VTE(pis,1)); % fewer coexisting
psch(2) = 1/length(pis)*sum(sum(CXTE(pis,:),2)==VTE(pis,1)); % same coexisting
psch(3) = 1/length(pis)*sum(sum(CXTE(pis,:),2)>VTE(pis,1)); % more coexisting

nio = find(VTE(1:cntinstE,3)==-1);
noch(1) = 1/length(nio)*sum(sum(CXTE(nio,:),2)<VTE(nio,1)); % fewer coexisting
noch(2) = 1/length(nio)*sum(sum(CXTE(nio,:),2)==VTE(nio,1)); % same coexisting
noch(3) = 1/length(nio)*sum(sum(CXTE(nio,:),2)>VTE(nio,1)); % more coexisting

nis = find(VTE(1:cntinstE,3)==-2);
nsch(1) = 1/length(nis)*sum(sum(CXTE(nis,:),2)<VTE(nis,1)); % fewer coexisting
nsch(2) = 1/length(nis)*sum(sum(CXTE(nis,:),2)==VTE(nis,1)); % same coexisting
nsch(3) = 1/length(nis)*sum(sum(CXTE(nis,:),2)>VTE(nis,1)); % more coexisting

figure
bar([sch' nsch' noch' poch' psch']','stack')
set(gca,'XTickLabel',{'any','- on self','- on others','+ on others','+ on self'})
ylabel('Frequency')
xlabel('Category of influence removed')
text(0.5,1.1,num2str(cntinstE))
text(1.7,1.1,num2str(length(nis)))
text(2.5,1.1,num2str(length(nio)))
text(3.7,1.1,num2str(length(pio)))
text(4.7,1.1,num2str(length(pis)))
text(0.5,0.9,'more')
text(0.5,0.5,'same')
text(0.5,0.1,'fewer')
ylim([0 1])

save(strcat(infile,'_CheckInflRemoval_SSvsEn_A'))
