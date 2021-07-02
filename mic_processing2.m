%%% Determining Minimal Inhibitory Concentration (MIC)
clear;close all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Printing options
printfigures=1;% 1 - "print"; 0 - "do not print";
frmt='-dpng';
res=200;%resolution in dpi
FS=19;%Font size for figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reading experimental data
% Open file dialog form
[fname,fpath,~]=uigetfile('*.xlsx');
% fname='Izoniazid  H37Rv  07.03.2017.xlsx';
% fname='DO-214  7074  15.04.2021.xlsx';
% fpath='C:\Users\Eugene\Dropbox\Summer2020\MIC\micRawData\';
% Read Excel file
[xlsdata,xlstext]=xlsread([fpath,fname]);
outName=['PROCESSED',fname];
% Test description
T=cell(6,1);
sz=size(xlstext);
if sz(1)>9
if length(xlstext([3,5:9],1))==6
    T=xlstext([3,5:9],1);
end
end
% Extracting fluorescence data
data=xlsdata(2:end-1,1:end-1);
[n,m]=size(data);
% Reshaping: 1% suspension without a drug should be first
data=[data(:,m),data(:,1:m-1)];
% Check whether 1% suspension do not grow and correct if needed
if mean(data(:,1))>=(mean(data(:,end))/2)
   data(:,1)=0.1*data(:,end);
end
% Extracting concentrations
conc=xlsdata(end,1:m-2);
% Converting array of concentration to labels
xlab(2:m-1)=strsplit(num2str(conc));
xlab{m}='Ctr';
xlab{1}='Ctr(1%)';
%% Data processing

figure('NumberTitle', 'off', 'Name', 'Raw experimental data assessment');
boxplot(data)
hold on
for j=1:m
    plot(j,data(:,j),'o','color','black','MarkerSize',9);
end
xlim([0.5 m+0.5])
xlabel('Concentration, mg/ml')
ylabel('Recorded fluorescence, RFU')
set(gca,'XTick',[1:m],'XTickLabel',xlab,'FontSIze',FS);
title([T{4},' ',T{5},' ',T{6}],'FontSize',FS)


figure('NumberTitle', 'off', 'Name', 'Raw experimental data assessment');
subplot(3,1,[1 2])
data_clear=[];
g=[];
for j=1:m
     plot(j,data(:,j),'o','color','magenta')
    hold on
    % 75th percentiles of the sample data 
    q3=quantile(data(:,j),0.75);
    % 25th percentiles of the sample data
    q1=quantile(data(:,j),0.25);
    % Upper outlier boundary
    outUp=q3+1.5*(q3-q1);
    % Lower outlier boundary
    outDown=q1-1.5*(q3-q1);
    notOut=(data(:,j)>outDown)&(data(:,j)<outUp);
    % Data without outliers
    data_notOut=data(notOut,j);
    % Mann-Whitney U-test
    if j==1
        dataRef=data_notOut;
    else
       [p(j-1),h(j-1)]=ranksum(dataRef,data_notOut);
    end
    
    
    % Cleared data's median and quartiles
    DataMedian(j)=median(data_notOut);
    q3=quantile(data_notOut,0.75);
    q1=quantile(data_notOut,0.25);
    % Length of data without outliers
    len=length(data_notOut);
    % Confidence intervals
    ConfUpDown(:,j)=[DataMedian(j)+1.57*(q3-q1)/sqrt(len);...
                     DataMedian(j)-1.57*(q3-q1)/sqrt(len)];
    data_clear=[data_clear;data_notOut];
    gj=repmat({xlab{j}},len,1);
    g=[g;gj];
%     plot([j;j],ConfUpDown(:,j),'.-')
end
hold on
%boxplot(data)
boxplot(data_clear,g,'Notch','marker')

MaxCtr1=repmat(ConfUpDown(1,1),1,m-2);
MinCtr1=repmat(ConfUpDown(2,1),1,m-2);

maxCI=max(MinCtr1,ConfUpDown(2,2:end-1));
minCI=min(MaxCtr1,ConfUpDown(1,2:end-1));

overlap_ind=(minCI>=maxCI);

lower_ind=DataMedian(2:end-1)<DataMedian(1);
indCtr1=sign(overlap_ind+lower_ind);
ind=find(indCtr1==1);
indMIC=max(ind);

indMWMIC=find((h==1)&DataMedian(2:end)>DataMedian(1))-1;
if length(indMWMIC>0)
    indMWMIC=indMWMIC(1);
end

plot([1,m],[ConfUpDown(2,1),ConfUpDown(2,1)],'-.','color','black')
plot([1,m],[ConfUpDown(1,1),ConfUpDown(1,1)],'-.','color','black')
xlabel('Concentration, mg/ml')
ylabel('Recorded fluorescence, RFU')
xlim([0.5 m+0.5])

set(gca,'XTick',[1:m],'XTickLabel',xlab,'FontSIze',FS);
title([T{4},' ',T{5},' ',T{6}],'FontSize',FS)

subplot(3,1,3)
xlim([0.5 m+0.5])
ylim([-2 2])


for j=1:m-2
    if overlap_ind(j)==1
        %plot(j+1,0.5,'.','Marker','none')
        text(j+1,0.5,'=','FontSize',FS)
        hold on
    else
        if lower_ind(j)==1
            %plot(j+1,0.5,'.','Marker','none')
            text(j+1,0.5,'<','FontSize',FS)
            hold on
        else
            %plot(j+1,0.5,'.','Marker','none')
            text(j+1,0.5,'>','FontSize',FS)
            hold on
        end
    end
    text(j+1,-1.5,['H_',num2str(h(j))],'FontSize',FS)
    text(j+0.75,-2.5,['p=',num2str(p(j),1)],'FontSize',FS-5)
end
axis off

if length(indMIC)>0
    text(0.5,1.5,['Median - notches overlapping test:'...
                  ' MIC=',num2str(conc(indMIC)),' mg/ml'],'FontSize',FS)
else
    text(0.1,1.5,['Median - notches overlapping test:'...
                  'MIC>',num2str(conc(1)),' mg/ml'],'FontSize',FS)
end

if indMWMIC>0
    text(0.5,-0.5,['Mann-Whitney U-test:'...
                  ' MIC=',num2str(conc(indMWMIC)),' mg/ml'],'FontSize',FS)
else
    text(0.1,-0.5,['Mann-Whitney U-test:'...
                  'MIC>',num2str(conc(1)),' mg/ml'],'FontSize',FS)
end

%% Test of Hill's kinetics

% Normed fluorescence: 0 for the 1% control, 1 for the normal control
Mn=(DataMedian-DataMedian(1))/(DataMedian(end)-DataMedian(1));

figure('NumberTitle', 'off', 'Name', 'Test of Hill`s kinetics');
MLog=1./Mn(2:end-1)-1;
ind=(MLog>0);
logc=2:10;
plot(logc(ind),log2(MLog(ind)),'o','color','red','LineWidth',1,'MarkerSize',9)
hold on
pMn=polyfit(logc(ind),log2(MLog(ind)),1);
plot(logc(ind),polyval(pMn,logc(ind)),'LineWidth',1.5,'color','blue')
xlim([0.5 m+0.5])
xlabel('Concentration (log_2 scale), mg/ml')
ylabel('log_2(M_{normed}^{-1}-1)')

set(gca,'XTick',[1:m],'XTickLabel',xlab,'FontSIze',FS);
title([T{4},' ',T{5},' ',T{6}],'FontSize',FS)

% % Normed fluorescence: 1 for the normal control
Mn=DataMedian/DataMedian(end);
modelfun = @(b,x) b(1)+(b(2)-b(1))./(1+exp(b(3)*(x-b(4))));
beta0=[0;1;2;mean(log(conc))];
% Using robust fitting options
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
[beta,R,J,CovB,MSE]=nlinfit(log(conc),Mn(2:end-1),modelfun,beta0);
ka=log(conc(1)/conc(end))/(m-3);
kb=log(conc(end))-2*ka;
nodes=linspace(1,m,m*100+1);
lcfit=ka*nodes+kb;
[Mpred,delta] = nlpredci(modelfun,lcfit,beta,R,'Covar',CovB,...
                         'MSE',MSE,'SimOpt','on');
lower = Mpred - delta;
upper = Mpred + delta;

% Determination of MIC
if (beta(2)-Mn(1))/(Mn(1)-beta(1))>0
    logMIC=beta(4)+log((beta(2)-Mn(1))/(Mn(1)-beta(1)))/beta(3);
    MIC=exp(logMIC);
else
    MIC=[];
end
logMICplat=beta(4)+log(0.95/0.01)/beta(3);
MICplat=exp(logMICplat);

figure('NumberTitle', 'off', 'Name', 'MIC from the continuous approximation');

subplot(1,4,[1:3])
plot([m-1:-1:2],Mn(2:end-1)*DataMedian(end),'x','color','red','LineWidth',1)
hold on
boxplot(flipud(data_clear),flipud(g),'Notch','marker')
plot(nodes,DataMedian(end)*Mpred,'-','color','red','LineWidth',1.5)
plot(nodes,DataMedian(end)*[lower;upper],':','color','red','LineWidth',1)

plot([1,m],[DataMedian(1),DataMedian(1)],'--','color','black')
plot([1,m],[ConfUpDown(2,1),ConfUpDown(2,1)],'-.','color','black')
plot([1,m],[ConfUpDown(1,1),ConfUpDown(1,1)],'-.','color','black')

if length(MIC)>0
    plot((logMIC-kb)/ka,modelfun(beta,logMIC)*DataMedian(end),'*','color','green','LineWidth',1)
end
plot((logMICplat-kb)/ka,modelfun(beta,logMICplat)*DataMedian(end),'*','color','magenta','LineWidth',1)

xlim([0.8*min(conc),1.2*max(conc)])
xlabel('Concentration, mg/ml')
ylabel('Recorded fluorescence, RFU')
xlim([0.5 m+0.5])

set(gca,'XTick',[1:m],'XTickLabel',fliplr(xlab),'FontSIze',FS);


subplot(1,4,4)
for j=1:6
text(0.1,1.1-0.1*j,num2str(cell2mat(T(j))),'FontSize',FS)
end
axis off
text(0.1,1.2-0.75,['-----------------------------------------'],'FontSize',FS)
%text(0.1,1.2-0.8,['Fl(c)=b_1-b_2tanh\{r_{IC}[ln(c/IC_{0.5})]\}'],'FontSize',FS)
text(0.1,1.2-0.9,['C_{0.5}=',num2str(exp(beta(4)),3),' mg/ml'],'FontSize',FS)
text(0.1,1.2-1,['\alpha=',num2str(beta(3),3)],'FontSize',FS)
if length(MIC)>0
    text(0.1,0.1,['MIC_{fit}=',num2str(MIC,3),' mg/ml'],'FontSize',FS)
else
    text(0.1,0.1,['MIC_{fit}>',num2str(max(conc)),' mg/ml'],'FontSize',FS)
end
text(0.1,0,['MIC_{EC95}=',num2str(MICplat,3),' mg/ml'],'FontSize',FS)


%% Printing figures
if printfigures==1
    figure(1)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 8],'PaperSize',[12 8])
    nameplate=['plateRaw',outName(1:end-4),'.',frmt(3:end)];
    print(nameplate,frmt,['-r',num2str(res)])
    figure(2)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 8],'PaperSize',[12 8])
    nameplate=['plate',outName(1:end-4),'.',frmt(3:end)];
    print(nameplate,frmt,['-r',num2str(res)])  
    figure(3)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 8],'PaperSize',[12 8])
    nameplate=['plateHill',outName(1:end-4),'.',frmt(3:end)];
    print(nameplate,frmt,['-r',num2str(res)])  
    figure(4)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 6],'PaperSize',[14 6])
    namefit=['fit',outName(1:end-4),'.',frmt(3:end)];
    print(namefit,frmt,['-r',num2str(res)]) 
end
