%%% Determining Minimal Inhibitory Concentration (MIC)
clear;close all;clc;
Rodbard=@(b,x)  b(1)+(b(2)-b(1))./(1+(x/b(3)).^b(4));
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
[xlsdata,xlstext]=xlsread([fpath,fname]);
%[xlsdata,xlstext]=xlsread('Fluor16table.xlsx');

% Extracting  data
data=xlsdata(1:end-1,1:9);
conc=xlsdata(end,1:9);
%% Date without visible outliers
ind=[1:3,5:9]; 
% Median and median average deviation
cM=median(data);
% Data fitting
beta0=[min(cM),max(cM),mean(conc),1];
beta=nlinfit(conc(ind),cM(ind),Rodbard,beta0);
Mc=(cM(1:9)-beta(1))/(beta(2)-beta(1));
Mmadc=mad(data)/(beta(2)-beta(1));
errorbar([1:9],Mc,Mmadc(1:9),'s','color','blue','LineWidth',1)
hold on
x=linspace(1,9,1001);
p=polyfit([1:9],log(conc),1);
plot(x,(Rodbard(beta,exp(polyval(p,x)))-beta(1))/(beta(2)-beta(1)),...
     '--','color','blue','LineWidth',1.5);
ylim([min(Mc-Mmadc) max(Mc+Mmadc)])
xlabel('Concentration, \mug/ml')
ylabel('Normed drug response')
set(gca,'XTick',[1:9],'XTickLabel',conc)
set(gca,'FontSize',FS,'FontName','Times');

grid on

set(gca,'FontSize',FS,'FontName','Times');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 12 8],'PaperSize',[12 8])

if printfigures==1
    print -dpdf -r300 'fig_MIC'
end

% %%
% pconc=polyfit([1:9],log10(conc),1)
% %
% IC95clr=10^(polyval(pconc,1.024))
% IC95fluo=10^(polyval(pconc,1.448))
% %
% IC90clr=10^(polyval(pconc,1.584))
% IC90fluo=10^(polyval(pconc,1.896))
% %
% IC75clr=10^(polyval(pconc,2.408))
% IC75fluo=10^(polyval(pconc,2.544))
% %
% IC50clr=10^(polyval(pconc,3.232))
% IC50fluo=10^(polyval(pconc,3.2))