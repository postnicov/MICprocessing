clear; close all;clc;FS=18;

TR=readtable('Red.csv');
TG=readtable('Green.csv');
TB=readtable('Blue.csv');
%
TR.R=10.^(-TR.Acorr);
TG.G=10.^(-TG.Acorr);
TB.B=10.^(-TB.Acorr);
%
M=max([TR.R,TG.G,TB.B]);
Norm=max([1,max(M)]);
TR.R=TR.R/Norm;
TG.G=TG.G/Norm;
TB.B=TB.B/Norm;
%
for j=1:8
    TB.Y(12*(j-1)+1:12*j,1)=(9-j)*ones(12,1);
end
for j=1:96
    plot(TB.Column(j), TB.Y(j),'.','MarkerSize',90,'color',[TR.R(j),TG.G(j),TB.B(j)])
    hold on
end
xlim([0 13])
ylim([0 9])
set(gca,'XTick',[1:12],'YTick',[1:8],'YTickLabel',{'H','G','F','E','D','C','B','A'},'FontSize',FS)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 18 12],'PaperSize',[18 12])
print('FigPlate','-dpng','-r150')
%%
Lab=rgb2lab([TR.R,TG.G,TB.B]);
Row=TB.Y;
Column=TB.Column;
L=Lab(:,1);
a=Lab(:,2);
b=Lab(:,3);
TLab=table(Row,Column,L,a,b);

%
for j=1:12;
    median_ab(j,1)=median(TLab.a(TLab.Column==j));
    amad(j)=mad(TLab.a(TLab.Column==j),1);
    median_ab(j,2)=median(TLab.b(TLab.Column==j));
    bmad(j)=mad(TLab.b(TLab.Column==j),1);
end
%
figure
subplot(2,1,1)
for j=1:12
    aj=TLab.a(TLab.Column==j);
    plot(j,aj,'o','color','magenta','MarkerSize',4,'LineWidth',0.5);
    hold on
end    
errorbar([1:12],median_ab(:,1),1.5*amad,'*:',...
    'color','magenta','MarkerSize',8,'LineWidth',1)
xlim([0.5 12.5])
xlabel('Column''s number')
ylabel('a*')
set(gca,'FontSize',FS,'FontName','Times');
subplot(2,1,2)
for j=1:12
    bj=TLab.b(TLab.Column==j);
    plot(j,bj,'o','color','blue','MarkerSize',4,'LineWidth',0.5);
    hold on
end    
errorbar([1:12],median_ab(:,2),1.5*bmad,'*:',...
    'color','blue','MarkerSize',8,'LineWidth',1)
xlim([0.5 12.5])
xlabel('Column''s number')
ylabel('b*')
set(gca,'FontSize',FS,'FontName','Times');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 12],'PaperSize',[16 12])
print('Figab','-dpdf','-r150')

save matlabLab

%%
conc=readtable('Concentrations.csv');
for j=1:12
    ind=find(TLab.Column==j);
    Tout(:,j)=TLab.a(ind);
end
Tout=flipud(Tout);
Tout(9,1:9)=table2array(conc)';
Tout=array2table(Tout);
writetable(Tout,'Tout.xlsx')