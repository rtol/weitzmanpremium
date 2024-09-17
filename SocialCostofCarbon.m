
%SocialCostofCarbon
%The Climate Framework for Uncertainty, Negotiation and Distribution,
%version 4.1-matlab-global
%
%This script is part of FUND 4.1 MG
%It computes the social cost of carbon
%
%Richard Tol, 28 March 2018
%This code is protected by the MIT License

dimpact = impactd - impact;

for i=1:NImpact,
    dimpabs(i,:,:) = squeeze(dimpact(i,:,:)).*Y;
end

gy = zeros(NYear,NScen);
for t=2:NYear,
    for s=1:NScen,
        gy(t,s) = (YpC(t,s)-YpC(t-1,s))/YpC(t-1,s);
    end
end

gp = zeros(NYear,NScen);
for t=2:NYear,
    for s=1:NScen,
        gp(t,s) = (Population(t,s)-Population(t-1,s))/Population(t-1,s);
    end
end

NDR = 8;
PRTP = [0.001 0.003 0.010 0.015 0.020 0.030 0.040 0.050];
%NDR = 4;
%PRTP = [0.00125 0.005 0.00875 0.0125];
NRA = 6;
RA = [0.5 1.0 1.3 1.5 2.0 2.5];
%NRA = 4;
%RA = [1.25 1 0.75 0.5];
NCDR = 6;
DR = [0.072823464 0.070 0.0545 0.030 0.025 0.02];

DDR = [2.24800683 0.28982023 0.06113661];

df = zeros(NYear,NScen,NDR,NRA);
df2 = zeros(NYear,NCDR);
dfw = zeros(NYear,1);
for s=1:NScen,
    for p=1:NDR,
        for r=1:NRA,
            df(SCCYear,s,p,r)=1;
        end
    end
end
for r=1:NCDR,
    df2(SCCYear,r)=1;
end
dfw(SCCYear)=1;

for t=SCCYear+1:NYear,
    for s=1:NScen,
        for p=1:NDR,
            for r=1:NRA
                df(t,s,p,r) = df(t-1,s,p,r)/(1+PRTP(p)+gp(t-1,s)+RA(r)*gy(t-1,s));
            end
        end
    end
    for r=1:NCDR,
        df2(t,r) = df2(t-1,r)/(1+DR(r) + 0.00857974);
    end
    draux = (1-DDR(3))*DDR(1)/(DDR(2) + t-1-SCCYear) + 0.857974; %0.86 is the average rate of population growth 2014-23
    dfw(t) = dfw(t-1)/(1+0.01*draux); 
end

SCC = zeros(NImpact, NScen, NDR, NRA);
SCC2 = zeros(NImpact, NScen, NCDR);
SCCWeitzman = zeros(NImpact, NScen);

for i=1:NImpact,
    for s=1:NScen,
        for p=1:NDR,
            for r=1:NRA,
                SCC(i,s,p,r) = squeeze(df(:,s,p,r))'*dimpabs(i,:,s)';
            end
        end
        for p=1:NCDR,
            SCC2(i,s,p) = df2(:,p)'*dimpabs(i,:,s)';
        end
        SCCWeitzman(i,s) = dfw'*dimpabs(i,:,s)';
    end
end

SCC = -0.01*SCC/1000000;
SCC2 = -0.01*SCC2/1000000;
SCCWeitzman = -0.01*SCCWeitzman/1000000;

s = PrintTable;
s.addRow('Time pref \ Risk aversion', RA(1), RA(2), RA(3), RA(4), RA(5), RA(6));
s.addRow(num2str(PRTP(1),5),num2str(SCC(1,1,1,1),7),num2str(SCC(1,1,1,2),7),num2str(SCC(1,1,1,3),7),num2str(SCC(1,1,1,4),7),num2str(SCC(1,1,1,5),7),num2str(SCC(1,1,1,6),7));
s.addRow(num2str(PRTP(2),5),num2str(SCC(1,1,2,1),7),num2str(SCC(1,1,2,2),7),num2str(SCC(1,1,2,3),7),num2str(SCC(1,1,2,4),7),num2str(SCC(1,1,2,5),7),num2str(SCC(1,1,2,5),7));
s.addRow(num2str(PRTP(3),5),num2str(SCC(1,1,3,1),7),num2str(SCC(1,1,3,2),7),num2str(SCC(1,1,3,3),7),num2str(SCC(1,1,3,4),7),num2str(SCC(1,1,3,5),7),num2str(SCC(1,1,3,5),7));
s.addRow(num2str(PRTP(4),5),num2str(SCC(1,1,4,1),7),num2str(SCC(1,1,4,2),7),num2str(SCC(1,1,4,3),7),num2str(SCC(1,1,4,4),7),num2str(SCC(1,1,4,5),7),num2str(SCC(1,1,4,5),7));
s.addRow(num2str(PRTP(5),5),num2str(SCC(1,1,5,1),7),num2str(SCC(1,1,5,2),7),num2str(SCC(1,1,5,3),7),num2str(SCC(1,1,5,4),7),num2str(SCC(1,1,5,5),7),num2str(SCC(1,1,5,5),7));
s.addRow(num2str(PRTP(6),5),num2str(SCC(1,1,6,1),7),num2str(SCC(1,1,6,2),7),num2str(SCC(1,1,6,3),7),num2str(SCC(1,1,6,4),7),num2str(SCC(1,1,6,5),7),num2str(SCC(1,1,6,5),7));
s.addRow(num2str(PRTP(7),5),num2str(SCC(1,1,7,1),7),num2str(SCC(1,1,7,2),7),num2str(SCC(1,1,7,3),7),num2str(SCC(1,1,7,4),7),num2str(SCC(1,1,7,5),7),num2str(SCC(1,1,7,5),7));
s.addRow(num2str(PRTP(8),5),num2str(SCC(1,1,8,1),7),num2str(SCC(1,1,8,2),7),num2str(SCC(1,1,8,3),7),num2str(SCC(1,1,8,4),7),num2str(SCC(1,1,8,5),7),num2str(SCC(1,1,8,5),7));
disp('Social cost of carbon ($/tC)')
disp('alternative rates of time preference (rows) and risk aversion (columns)')
str = ['model = Tol, scenario = SRES A1, year = ' num2str(SCCYear+StartYear)];
disp(str)
s.display

line = sprintf('\n');
disp(line)

t = PrintTable;
t.addRow('Model \ Scenario', 'SRES A1', 'SRES A2', 'SRES B1', 'SRES B2', 'SSP1', 'SSP2', 'SSP3', 'SSP4', 'SSP5');
t.addRow('Tol parabola',num2str(SCC(1,1,3,2),7),num2str(SCC(1,2,3,2),7),num2str(SCC(1,3,3,2),7),num2str(SCC(1,4,3,2),7),num2str(SCC(1,5,3,2),7),num2str(SCC(1,6,3,2),7),num2str(SCC(1,7,3,2),7),num2str(SCC(1,8,3,2),7),num2str(SCC(1,9,3,2),7));
t.addRow('Weitzman (6)',num2str(SCC(2,1,3,2),7),num2str(SCC(2,2,3,2),7),num2str(SCC(2,3,3,2),7),num2str(SCC(2,4,3,2),7),num2str(SCC(2,5,3,2),7),num2str(SCC(2,6,3,2),7),num2str(SCC(2,7,3,2),7),num2str(SCC(2,8,3,2),7),num2str(SCC(2,9,3,2),7));
t.addRow('Weitzman (7)',num2str(SCC(3,1,3,2),7),num2str(SCC(3,2,3,2),7),num2str(SCC(3,3,3,2),7),num2str(SCC(3,4,3,2),7),num2str(SCC(3,5,3,2),7),num2str(SCC(3,6,3,2),7),num2str(SCC(3,7,3,2),7),num2str(SCC(3,8,3,2),7),num2str(SCC(3,9,3,2),7));
t.addRow('Nordhaus',num2str(SCC(4,1,3,2),7),num2str(SCC(4,2,3,2),7),num2str(SCC(4,3,3,2),7),num2str(SCC(4,4,3,2),7),num2str(SCC(4,5,3,2),7),num2str(SCC(4,6,3,2),7),num2str(SCC(4,7,3,2),7),num2str(SCC(4,8,3,2),7),num2str(SCC(4,9,3,2),7));
t.addRow('Hope',num2str(SCC(5,1,3,2),7),num2str(SCC(5,2,3,2),7),num2str(SCC(5,3,3,2),7),num2str(SCC(5,4,3,2),7),num2str(SCC(5,5,3,2),7),num2str(SCC(5,6,3,2),7),num2str(SCC(5,7,3,2),7),num2str(SCC(5,8,3,2),7),num2str(SCC(5,9,3,2),7));
t.addRow('Ploeg',num2str(SCC(6,1,3,2),7),num2str(SCC(6,2,3,2),7),num2str(SCC(6,3,3,2),7),num2str(SCC(6,4,3,2),7),num2str(SCC(6,5,3,2),7),num2str(SCC(6,6,3,2),7),num2str(SCC(6,7,3,2),7),num2str(SCC(6,8,3,2),7),num2str(SCC(6,9,3,2),7));
t.addRow('Golosov',num2str(SCC(7,1,3,2),7),num2str(SCC(7,2,3,2),7),num2str(SCC(7,3,3,2),7),num2str(SCC(7,4,3,2),7),num2str(SCC(7,5,3,2),7),num2str(SCC(7,6,3,2),7),num2str(SCC(7,7,3,2),7),num2str(SCC(7,8,3,2),7),num2str(SCC(7,9,3,2),7));
t.addRow('Tol piecewise linear',num2str(SCC(8,1,3,2),7),num2str(SCC(8,2,3,2),7),num2str(SCC(8,3,3,2),7),num2str(SCC(8,4,3,2),7),num2str(SCC(8,5,3,2),7),num2str(SCC(8,6,3,2),7),num2str(SCC(8,7,3,2),7),num2str(SCC(8,8,3,2),7),num2str(SCC(8,9,3,2),7));
disp('Social cost of carbon ($/tC)')
disp('alternative impact models (rows) and scenarios (columns)')
str =['pure rate of time preference = 0.01; rate of risk aversion = 1, year = ' num2str(SCCYear+StartYear)]; 
disp(str)
t.display

%%
Drupp = csvread('Drupp.csv');
Drupp = Drupp(~isnan(Drupp(:,1))&~isnan(Drupp(:,2)),:);
ND = size(Drupp,1);

Falk = csvread('Falk-GPS.csv');
NH = 65;
NF = (size(Falk,1)-NH)/5;

FalkInd = readtable('FalkInd.csv');
%FalkInd = readtable('FalkPop.csv');
%FalkInd = readtable('FalkNAWE.csv');
FalkInd = table2array(FalkInd);
NI = size(FalkInd,1);
FalkCountry = csvread('FalkCountry.csv');
NC = size(FalkCountry,1);
FalkAge = csvread('FalkAge.csv');
NA = size(FalkAge,1);

Pop = csvread('population.csv');

prefs = [Drupp; Falk];

NPref = ND + NF + NF + NF + NF + NF + NH;

df3 = zeros(NYear,NScen,ND);
for s=1:NScen,
    for p=1:NPref,
        df3(SCCYear,s,p)=1;
    end
end

for t=SCCYear+1:NYear,
    for s=1:NScen,
        for p=1:NPref,
            df3(t,s,p) = df3(t-1,s,p)/(1+0.01*prefs(p,1)+gp(t-1,s)+prefs(p,2)*gy(t-1,s));
        end
    end
end

dfi = zeros(NYear,NScen,NI);
for s=1:NScen,
    for p=1:NI,
        dfi(SCCYear,s,p)=1;
    end
end

for t=SCCYear+1:NYear,
    for s=1:NScen,
        for p=1:NI,
            dfi(t,s,p) = dfi(t-1,s,p)/(1+0.01*FalkInd(p,1)+gp(t-1,s)+FalkInd(p,2)*gy(t-1,s));
        end
    end
end

dfc = zeros(NYear,NScen,NC);
for s=1:NScen,
    for p=1:NC,
        dfc(SCCYear,s,p)=1;
    end
end

for t=SCCYear+1:NYear,
    for s=1:NScen,
        for p=1:NC,
            dfc(t,s,p) = dfc(t-1,s,p)/(1+0.01*FalkCountry(p,1)+gp(t-1,s)+FalkCountry(p,2)*gy(t-1,s));
        end
    end
end

dfa = zeros(NYear,NScen,NA);
for s=1:NScen,
    for p=1:NA,
        dfa(SCCYear,s,p)=1;
    end
end

for t=SCCYear+1:NYear,
    for s=1:NScen,
        for p=1:NA,
            dfa(t,s,p) = dfa(t-1,s,p)/(1+0.01*FalkAge(p,1)+gp(t-1,s)+FalkAge(p,2)*gy(t-1,s));
        end
    end
end

SCC3 = zeros(NImpact, NScen, NPref);
SCCI = zeros(NImpact, NScen, NI);
SCCC = zeros(NImpact, NScen, NC);
SCCA = zeros(NImpact, NScen, NA);

for i=1:NImpact,
    for s=1:NScen,
        for p=1:NPref,
            SCC3(i,s,p) = squeeze(df3(:,s,p))'*dimpabs(i,:,s)';
        end
	for p=1:NI,
            SCCI(i,s,p) = squeeze(dfi(:,s,p))'*dimpabs(i,:,s)';
        end
	for p=1:NC,
            SCCC(i,s,p) = squeeze(dfc(:,s,p))'*dimpabs(i,:,s)';
        end
	for p=1:NA,
            SCCA(i,s,p) = squeeze(dfa(:,s,p))'*dimpabs(i,:,s)';
        end
    end
end

SCC3 = -0.01*SCC3/1000000;
SCCI = -0.01*SCCI/1000000;
SCCC = -0.01*SCCC/1000000;
SCCA = -0.01*SCCA/1000000;

%%
for i=1:5,
    SCCDrupp = squeeze(SCC3(9,i+4,ND+1:ND+NF-2));
    SCCWeighted = squeeze(SCC3(9,i+4,ND+NF+1:ND+2*NF-2));
    SCCLimited = squeeze(SCC3(9,i+4,ND+2*NF+1:ND+3*NF-2));
    SCCObserved = squeeze(SCC3(9,i+4,ND+3*NF+1:ND+4*NF-2));
    SCCImputed = squeeze(SCC3(9,i+4,ND+4*NF+1:ND+5*NF-2));
    SCCHofstede = squeeze(SCC3(9,i+4,ND+5*NF+1:ND+5*NF+NH-2));

    SCCD(i) = sum(SCCDrupp)/(NF-2);
    SCCW(i) = Pop'*SCCWeighted/sum(Pop);
    SCCL(i) = sum(SCCLimited)/(NF-2);
    SCCO(i) = sum(SCCObserved)/(NF-2);
    SCCI(i) = sum(SCCImputed)/(NF-2);
    SCCH(i) = sum(SCCHofstede)/(NH-2);
end
%% repeat for impact function
for i=1:10,
    SCCDrupp = squeeze(SCC3(i,6,ND+1:ND+NF-2));
    SCCWeighted = squeeze(SCC3(i,6,ND+NF+1:ND+2*NF-2));
    SCCLimited = squeeze(SCC3(i,6,ND+2*NF+1:ND+3*NF-2));
    SCCImputed = squeeze(SCC3(i,6,ND+3*NF+1:ND+4*NF-2));
    SCCObserved = squeeze(SCC3(i,6,ND+4*NF+1:ND+5*NF-2));
    SCCHofstede = squeeze(SCC3(i,6,ND+5*NF+1:ND+5*NF+NH-2));

    SCCD2(i) = sum(SCCDrupp)/(NF-2);
    SCCW2(i) = Pop'*SCCWeighted/sum(Pop);
    SCCL2(i) = sum(SCCLimited)/(NF-2);
    SCCO2(i) = sum(SCCObserved)/(NF-2);
    SCCI2(i) = sum(SCCImputed)/(NF-2);
    SCCH2(i) = sum(SCCHofstede)/(NH-2);
end