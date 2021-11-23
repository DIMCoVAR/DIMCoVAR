%DIMCoVAR Consortium

% REFERENCE (Please cite):
%Wastewater and marine bioindicators surveillance to anticipate COVID-19 prevalence and to explore SARS-CoV-2 diversity by next generation sequencing: one-year study

% Code authorship
% Manuel Pájaro [1,2] (manuel.pajaro@udc.es). Irene Otero-Muras [1]. Noelia Fajar [1]. Antonio A. Alonso [1].
% [1] CSIC Spanish National Research Council
% [2] UDC University of A Coruña


% Main exe to generate figure predictions
% Here parameters fixed for the following locations/date:
% BAIONA  20-1-2021
% MELIDE  12-9-2020
% ARES SIRO 10-9-2020
% (FIGURE IN MAIN TEXT)



load DIMCoVARdata.mat %load data

figure
%% BAIONA
hab = 30000;
%Data
f= [2021 01 20];
[~, k]=min(abs(datetime(Fecha)-datetime(f))); 
B=Fecha(k:k+6)'; %I WTP 
C=Fecha(k+7:k+13)'; %O Xunta
WTP_Baiona(WTP_Baiona == -999) = NaN;
WTP=WTP_Baiona(k:k+6)';
dateWTP=B;
RealData=I_Baiona(k+7:k+13)';%Real Data Infected detected by Xunta
ICum0=ICum14_Baiona(k+7:k+13)'; %Incidencia Acumulada a 14 días
cumI0 =ICum0(1); %é o valor inicial de Incidencia Acumulada do día que empezas a simular
dateO =C ;

% Parameters
alpha=1/14; % being 14 the cumulative incidence fixed as constant
beta=0.03; 
gamma = 0.04;
% Baiona 
hab = 30000; %% inhabitants, but should be transformed to susceptibles....
% Time
nDays = length(RealData);
Tgrid=0:1:nDays-1; % days 
% Initial condition % [S,I,R,O,OR] 
x0 = [hab-max(cumI0,WTP(1)),max(+cumI0,WTP(1)),0, cumI0,0]; 
% Number of realizations
nsimula = 1000;
% Maximum Infected number threshold
Imax = sum(x0(1:3));%restrinxir poboación total para que non haaxa mais I que reais
% call the function to obtain the realizations
simulation = SIRO_ssa(alpha,beta,gamma,x0,Imax,Tgrid,nsimula);
II = zeros(nsimula,7);
OO = II;
newO = zeros(nsimula,7);
%Inew(:,1) = x0(2);
newO(:,1) = x0(4) - cumI0;
ApproxError = zeros(nsimula,1);
ApproxErrorI = zeros(nsimula,1);

for i=1:nsimula
    II(i,:) = simulation{i}(2,:);
    OO(i,:) = simulation{i}(4,:);
    %InewT = simulation{i}(1,1:end-1)-simulation{i}(1,2:end);
    newOT = simulation{i}(4,2:end)-simulation{i}(4,1:end-1)+...
              simulation{i}(5,2:end)-simulation{i}(5,1:end-1);
    %Inew(i,2:end) = InewT;
    newO(i,2:end) = newOT;
    ApproxError(i) = sum(abs(newO(i,:)-RealData));
    indO = find(WTP>=0);
    if indO > 0
        ApproxErrorI(i) = sum(abs(II(i,indO)-WTP(indO)));
    else
        ApproxErrorI(i)=0;
    end
end

acumulados =ICum0; 

subplot (3,2,6) 
% Xunta data Observed people
hold on 
plot(C,acumulados,'ks','MarkerSize',8,'LineWidth',3)
hold on 
plot(C,mean(OO,1),'k-','LineWidth',2)
hold on 
plot(C,mean(OO,1)+std(OO),'k--','LineWidth',1.5)
hold on 
plot(C,mean(OO,1)-std(OO),'k--','LineWidth',1.5)
box on
xlabel('Day','Fontsize',14,'Interpreter','latex')
ylabel('Infected','Fontsize',14,'Interpreter','latex')
title('Baiona Io 2021','Fontsize',14,'Interpreter','latex')
set(gca,'FontName', 'Times New Roman','FontSize',14,'TickLabelInterpreter','latex','TitleFontWeight', 'normal')

%WTP data Infected waters 
subplot (3,2,5) 
hold on 
plot(B,WTP,'bo','MarkerSize',8,'LineWidth',3)
hold on 
plot(B,mean(II,1),'b-','LineWidth',2)
hold on 
plot(B,mean(II,1)+std(II),'b--','LineWidth',1.5)
hold on 
plot(B,mean(II,1)-std(II),'b--','LineWidth',1.5)
box on
xlabel('Day','Fontsize',14,'Interpreter','latex')
ylabel('Infected','Fontsize',14,'Interpreter','latex')
title('Baiona I 2021','Fontsize',14,'Interpreter','latex')
set(gca,'FontName', 'Times New Roman','FontSize',14,'TickLabelInterpreter','latex','TitleFontWeight', 'normal')
%legend('Data O','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}',...
       %'Data WTP','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}')




%% MELIDE
hab = 14756;
%Data
f= [2020 09 12];
[~, k]=min(abs(datetime(Fecha)-datetime(f))); 
B=Fecha(k:k+6)'; %I WTP 
C=Fecha(k+7:k+13)'; %O Xunta
WTP_Melide(WTP_Melide == -999) = NaN;
WTP=WTP_Melide(k:k+6)';
dateWTP=B;
RealData=I_Melide(k+7:k+13)';%Real Data Infected detected by Xunta
ICum0=ICum14_Melide(k+7:k+13)'; %Incidencia Acumulada a 14 días
cumI0 =ICum0(1); %é o valor inicial de Incidencia Acumulada do día que empezas a simular
dateO =C ;

% Parameters
alpha=1/14; % being 14 the cumulative incidence fixed as constant
beta=0.005; 
gamma = 0.005;
% Time
nDays = length(RealData);
Tgrid=0:1:nDays-1; % days 
% Initial condition % [S,I,R,O,OR] 
x0 = [hab-max(cumI0,WTP(1)),max(+cumI0,WTP(1)),0, cumI0,0]; 
% Number of realizations
nsimula = 1000;
% Maximum Infected number threshold
Imax = sum(x0(1:3));%restrinxir poboación total para que non haaxa mais I que reais
% call the function to obtain the realizations
simulation = SIRO_ssa(alpha,beta,gamma,x0,Imax,Tgrid,nsimula);
II = zeros(nsimula,7);
OO = II;
newO = zeros(nsimula,7);
%Inew(:,1) = x0(2);
newO(:,1) = x0(4) - cumI0;
ApproxError = zeros(nsimula,1);
ApproxErrorI = zeros(nsimula,1);

for i=1:nsimula
    II(i,:) = simulation{i}(2,:);
    OO(i,:) = simulation{i}(4,:);
    %InewT = simulation{i}(1,1:end-1)-simulation{i}(1,2:end);
    newOT = simulation{i}(4,2:end)-simulation{i}(4,1:end-1)+...
              simulation{i}(5,2:end)-simulation{i}(5,1:end-1);
    %Inew(i,2:end) = InewT;
    newO(i,2:end) = newOT;
    ApproxError(i) = sum(abs(newO(i,:)-RealData));
    indO = find(WTP>=0);
    if indO > 0
        ApproxErrorI(i) = sum(abs(II(i,indO)-WTP(indO)));
    else
        ApproxErrorI(i)=0;
    end
end

acumulados =ICum0; 

subplot (3,2,4) 
% Xunta data Observed people
hold on 
plot(C,acumulados,'ks','MarkerSize',8,'LineWidth',3)
hold on 
plot(C,mean(OO,1),'k-','LineWidth',2)
hold on 
plot(C,mean(OO,1)+std(OO),'k--','LineWidth',1.5)
hold on 
plot(C,mean(OO,1)-std(OO),'k--','LineWidth',1.5)
box on
% xlabel('Day','Fontsize',14,'Interpreter','latex')
ylabel('Infected','Fontsize',14,'Interpreter','latex')
title('Melide Io 2020','Interpreter','latex')
%legend('Data O','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}',...
       %'Data WTP','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}')
set(gca,'FontName', 'Times New Roman','FontSize',14,'TickLabelInterpreter','latex','TitleFontWeight', 'normal')


%WTP data Infected waters 
subplot (3,2,3) 
hold on 
plot(B,WTP,'bo','MarkerSize',8,'LineWidth',3)
hold on 
plot(B,mean(II,1),'b-','LineWidth',2)
hold on 
plot(B,mean(II,1)+std(II),'b--','LineWidth',1.5)
hold on 
plot(B,mean(II,1)-std(II),'b--','LineWidth',1.5)
box on
% xlabel('Day','Fontsize',14,'Interpreter','latex')
ylabel('Infected','Fontsize',14,'Interpreter','latex')
%legend('Data O','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}',...
       %'Data WTP','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}')
title('Melide I 2020','Interpreter','latex')
set(gca,'FontName', 'Times New Roman','FontSize',14,'TickLabelInterpreter','latex','TitleFontWeight', 'normal')





%%ARES
hab = 14493;
%Data
f= [2020 09 10];
[~, k]=min(abs(datetime(Fecha)-datetime(f))); 
B=Fecha(k:k+6)'; %I WTP 
C=Fecha(k+7:k+13)'; %O Xunta
WTP_Ares(WTP_Ares == -999) = NaN;
WTP=WTP_Ares(k:k+6)';
dateWTP=B;
RealData=I_Ares(k+7:k+13)';%Real Data Infected detected by Xunta
ICum0=ICum14_Ares(k+7:k+13)'; %Incidencia Acumulada a 14 días
cumI0 =ICum0(1); %é o valor inicial de Incidencia Acumulada do día que empezas a simular
dateO =C ;

% Parameters
alpha=1/5; % being 14 the cumulative incidence fixed as constant
beta=0.005; 
gamma = 0.06;
% Ares 
hab = 14493; %% inhabitants, but should be transformed to susceptibles....
% Time
nDays = length(RealData);
Tgrid=0:1:nDays-1; % days 
% Initial condition % [S,I,R,O,OR] 
x0 = [hab-max(cumI0,WTP(1)),max(+cumI0,WTP(1)),0, cumI0,0]; 
% Number of realizations
nsimula = 1000;
% Maximum Infected number threshold
Imax = sum(x0(1:3));%restrinxir poboación total para que non haaxa mais I que reais
% call the function to obtain the realizations
simulation = SIRO_ssa(alpha,beta,gamma,x0,Imax,Tgrid,nsimula);
II = zeros(nsimula,7);
OO = II;
newO = zeros(nsimula,7);
%Inew(:,1) = x0(2);
newO(:,1) = x0(4) - cumI0;
ApproxError = zeros(nsimula,1);
ApproxErrorI = zeros(nsimula,1);

for i=1:nsimula
    II(i,:) = simulation{i}(2,:);
    OO(i,:) = simulation{i}(4,:);
    %InewT = simulation{i}(1,1:end-1)-simulation{i}(1,2:end);
    newOT = simulation{i}(4,2:end)-simulation{i}(4,1:end-1)+...
              simulation{i}(5,2:end)-simulation{i}(5,1:end-1);
    %Inew(i,2:end) = InewT;
    newO(i,2:end) = newOT;
    ApproxError(i) = sum(abs(newO(i,:)-RealData));
    indO = find(WTP>=0);
    if indO > 0
        ApproxErrorI(i) = sum(abs(II(i,indO)-WTP(indO)));
    else
        ApproxErrorI(i)=0;
    end
end

acumulados =ICum0; 


subplot (3,2,2) 
% Xunta data Observed people
hold on 
plot(C,acumulados,'ks','MarkerSize',8,'LineWidth',3)
hold on 
plot(C,mean(OO,1),'k-','LineWidth',2)
hold on 
plot(C,mean(OO,1)+std(OO),'k--','LineWidth',1.5)
hold on 
plot(C,mean(OO,1)-std(OO),'k--','LineWidth',1.5)
box on
%xlabel('Day','Fontsize',14,'Interpreter','latex')
ylabel('Infected','Fontsize',14,'Interpreter','latex')
title('Ares Io 2020','Interpreter','latex')
%legend('Data O','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}',...
       %'Data WTP','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}')
set(gca,'FontName', 'Times New Roman','FontSize',14,'TickLabelInterpreter','latex','TitleFontWeight', 'normal')


%WTP data Infected waters 
subplot (3,2,1) 
hold on 
plot(B,WTP,'bo','MarkerSize',8,'LineWidth',3)
hold on 
plot(B,mean(II,1),'b-','LineWidth',2)
hold on 
plot(B,mean(II,1)+std(II),'b--','LineWidth',1.5)
hold on 
plot(B,mean(II,1)-std(II),'b--','LineWidth',1.5)
box on
%xlabel('Day','Fontsize',14,'Interpreter','latex')
ylabel('Infected','Fontsize',14,'Interpreter','latex')
%legend('Data O','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}',...
%        'Data WTP','mean_{SSA}','mean_{SSA} + std_{SSA}','mean_{SSA} - std_{SSA}')
title('Ares I 2020','Interpreter','latex')
set(gca,'FontName', 'Times New Roman','FontSize',14,'TickLabelInterpreter','latex','TitleFontWeight', 'normal')

