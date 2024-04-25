clearvars
close all
clc

% READ MATRIX FOR ITALY
opts=detectImportOptions('MATRICE PENDOLARISMO 2011/matrix_pendo2011_10112014.txt');
opts.SelectedVariableNames=[1 3 8 15];
t=readtable('MATRICE PENDOLARISMO 2011/matrix_pendo2011_10112014.txt',opts);

tbr=(cell2mat(t{:,1})=='L' | t{:,3}==0);
t(tbr,:)=[];

from_prov=t{:,2};
to_prov=t{:,3};
trips=str2double(t{:,4});

pr = readtable('pendo_2011/prov_reg.csv');
[~,idx]=unique(pr,'rows','first');
pr=pr(idx,:);

OD=sparse(from_prov,to_prov,trips)';

OI = zeros(20,20); %matrix for regions
for R1 = 1:20
    prov_in_R1 = table2array(pr(pr.Var1==R1,2));
    for R2 = 1:20
        prov_in_R2 = table2array(pr(pr.Var1==R2,2));
        OI(R1,R2) = sum(OD(prov_in_R1,prov_in_R2),'all');
    end
end

P=bsxfun(@rdivide,OI,sum(OI,1));
    
% for plotting
regions = ["Piemonte",...
    "Valle d'Aosta",...
    "Lombardia",...
    "Trentino Alto Adige",...
    "Veneto",...
    "Friuli-Venezia Giulia",...
    "Liguria",...
    "Emilia Romagna",...
    "Toscana",...
    "Umbria",...
    "Marche",...
    "Lazio",...
    "Abruzzo",...
    "Molise",...
    "Campania",...
    "Puglia",...
    "Basilicata",...
    "Calabria",...
    "Sicilia",...
    "Sardegna"];


save('mobility.mat', 'OI', 'P',"regions")

np = 20;



r_vect = 1:np;
[xxx,yyy]=meshgrid(1:np,1:np);
mask=find(OI);


figure()
subplot(121)
scatter(xxx(mask),yyy(mask),250,log10(OI(mask)),'filled')
axis([0 np+1 0 np+1])
axis ij; axis square; box on
cbh=colorbar(gca,'YTick',0:1:6,'YTickLabel',10.^(0:1:6));
set(get(cbh,'Title'),'String',{'People moving';'from region l to region j'})
set(gca,'XTick',r_vect,'YTick',r_vect,'TickDir','out','XTickLabels',regions,'YTickLabels',regions)
xlabel('Region of origin'); ylabel('Region of destination')
subplot(122)
scatter(xxx(mask),yyy(mask),250,log10(P(mask)),'filled')
axis([0 np+1 0 np+1])
axis ij; axis square; box on
cbh=colorbar(gca,'YTick',-6:1:0,'YTickLabel',10.^(-6:1:0));
set(get(cbh,'Title'),'String',{'Movement probability';'from region l to region j'})
set(gca,'XTick',r_vect,'YTick',r_vect,'TickDir','out','XTickLabels',regions,'YTickLabels',regions)
xlabel('Region of origin'); ylabel('Region of destination')
set(findall(gcf,'-property','FontSize'),'FontSize',10)

