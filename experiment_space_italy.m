clear all
close all
clc

%% LOAD DATA
% 1. Load incidence
cases = readtable("data/cases_region.csv",'ReadVariableNames',true); 
Time = cases.data';

% 2. Treat incidence data
count_temp = table2array(cases(:,2:end)); count_temp(isnan(count_temp))=0;
% Reorganise Trentino Alto Adige
count = zeros(size(count_temp,1),20);
count(:,1:3) = count_temp(:,1:3);
count(:,4) = count_temp(:,20)+count_temp(:,21);
count(:,5:end) = count_temp(:,4:19);
count = diff([zeros(1,20); count],1);
cases_region = count'; clear count count_temp cases;
lim = 600; cases_region = cases_region(:,1:lim); Time = Time(1:lim);
cases_region(cases_region<0)=0;
Fp = smoothdata(cases_region, 2, 'movmean', [13 0]); Fp(Fp<0)=0;
Ft = sum(Fp,1);

% 3. Load Mobility and Population
load data/mobility P
population = readtable("data/population_regions.csv");
ResPop = table2array(population(:,2));

% Make matrices
C = full(P); clear P;
x=(1-diag(C)); % percentage of moving pop
Q=(C-diag(diag(C)))*diag(1./x); % extradiagonal fluxes


% 4. Load Google Mobility Data
load data/google-data.mat
a = find(Time_GMD == Time(1)); b = find(Time_GMD == Time(end));
GMD = GMD(:,a:b);
gmd_fill = fillmissing(GMD,'linear',2);
gmd_fill_smooth = smoothdata(gmd_fill, 2, 'movmean', 14);

% 5. Make exposure matrix
csi = x.*gmd_fill_smooth;
Z = compute_matrix_Z(ResPop,Q,csi);

%% Define parameters
% Max age of infection
q = 21; 
% Build generation time distribution
mean_GD = 5.2;
std_GD = 1.72;
a_beta = (mean_GD/std_GD)^2;
b_beta = mean_GD/a_beta;
beta = @(x) gampdf(x,a_beta,b_beta);
beta_x = beta(1:q)/sum(beta(1:q));      % generation times
% Survival function (on infectious compartment)
p = @(x) exp(-0.068*x);
p_x = p(0:q-1);       % survival probability    
sigma = exp(-0.068);  % fraction of reaching next dayphi = zeros(5,21);
% Ratio
phi = beta_x./p_x;

%% Compute convolution
alpha = compute_alpha(Fp,beta_x);
alpha_sumS = squeeze(sum(alpha,1))';

%% Compute effective reproduction numbers
% Parameters
par.delay = 0; 
par.init = 10; 
par.Np = 20000;
par.cv_r_0 = 1;
par.low_cv_r=0.75;
par.alpha_min = 0;
par.delta = 0.95;
par.lik = 'V1';

% Simulate
[R2,di,mo] = pf(Fp,par,Q,ResPop,csi,squeeze(alpha));

% Compute Leslie Projection matrix
[L2,K2] = compute_leslie_matrix(phi, p_x, R2.Q50, Z, 1);

% Compute transformation matrices;
X1 = zeros(20,20*21);
X2 = zeros(20,20*21);

for i = 1:20
    for j = 1:21
        X2(i,i+20*(j-1)) = 1/ResPop(i);
        X1(i,i+20*(j-1)) = 1;
    end
end

% Compute epidemicity and global effective RN
ES = compute_epidemic_subset(L2,20,21,X1,X2);
ES2 = compute_epidemic_subset(L2,20,21,X2,X2);
E2 = compute_epidemicity(L2,X1,X2);
E3 = compute_epidemicity(L2,X2,X2);
[RG2] = compute_global_RN(K2);

% Compute amplification envelope
Amax = zeros(3,size(L2,3));
for t = 1:size(L2,3)
    temp1 = zeros(1,50);
    temp2 = zeros(1,50);
    temp3 = zeros(1,50);
    Ltemp = squeeze(L2(:,:,t));
    for i = 1:length(temp1)
        ee = compute_epidemicity(Ltemp^i,X1,X2);
        temp1(i) = ee.E1;
        temp2(i) = ee.E2;
        ee = compute_epidemicity(Ltemp^i,X2,X2);
        temp3(i) = ee.E1;
    end
    Amax(1,t) = max(temp1,[],'omitmissing');
    Amax(2,t) = max(temp2,[],'omitmissing');
    Amax(3,t) = max(temp3,[],'omitmissing');
end
Amax_save = Amax;
Amax(:,RG2>1)=NaN;


%% PLOTS BELOW
% Render it nicely on Matlab (vectorial)
set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')

%% RGDRAW

% When the global RN is >1, the amplification envelope has no upper bound
RGdraw = RG2-1;
TRANS = find(diff(sign(RGdraw)));
first_non_nan = find(RGdraw>0,1,'first');
TRANS(first_non_nan-1)=first_non_nan;
TRANS(2:first_non_nan-2)=[];
if RGdraw(TRANS(2))>0
    TRANS(2)=[];
end
if mod(length(TRANS),2)
    TRANS(end+1) = length(RGdraw);
end

X1patch = TRANS(1:2:end);
X2patch = TRANS(2:2:end);

Xpatch = [X1patch; X2patch; X2patch; X1patch];


%% SPACE
regions = ["Piedmont",...
    "Aosta Valley",...
    "Lombardy",...
    "Trentino-Sudtirol",...
    "Veneto",...
    "Friuli-Venezia Giulia",...
    "Liguria",...
    "Emilia-Romagna",...
    "Tuscany",...
    "Umbria",...
    "Marche",...
    "Lazio",...
    "Abruzzo",...
    "Molise",...
    "Campania",...
    "Apulia",...
    "Basilicata",...
    "Calabria",...
    "Sicily",...
    "Sardinia"];

low_lim_A = 0.9; up_lim_A = 3;

color_list_2 = ["#e9c46a","#264653","#f4a261","#e76f51", "#2a9d8f"];

figure()
subplot(3,1,1)
yyaxis left
plot(Time, R2.Q50','linewidth',.25,'color',[0 0 0 0.4],'marker','none','linestyle','-');
hold on
plot(Time,ones(1,length(Time)),'--r','LineWidth',0.5)
ylabel('$\mathcal{R}_j$','interpreter','latex')
xlim([Time(1) Time(end)])
ylim([0 1.1*max(R2.Q50,[],'all')])
box off

yyaxis right
p1 = plot(Time,RG2,'Color',color_list_2(1),'LineWidth',1.25);
xlim([Time(1) Time(end)])
ylim([0 1.1*max(R2.Q50,[],'all')])
ylabel('$\mathcal{R}^G$','interpreter','latex')
box off
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = color_list_2(1);

subplot(3,1,2)
yyaxis left
hold on
p61 = plot(Time,E2.E1,'Color',color_list_2(2),'LineWidth',1);
p62 = plot(Time,E3.E1,'Color',color_list_2(2),'LineWidth',1,'LineStyle',':');
ylabel('$\mathcal{E}^{\ell^1}_t$','interpreter','latex')
plot(Time,ones(1,length(Time)),'--r','LineWidth',0.5)
Ypatch = repmat([low_lim_A low_lim_A up_lim_A up_lim_A]', 1,size(Xpatch,2));
p=patch(Xpatch,Ypatch,'k');
p.FaceAlpha=.2;p.EdgeAlpha=0;
xlim([Time(1) Time(end)])
ylim([low_lim_A up_lim_A])
hold off

yyaxis right
hold on
p63 = plot(Time,Amax(1,:),'Color',color_list_2(3),'LineWidth',1);
p64 = plot(Time,Amax(3,:),'Color',color_list_2(3),'LineWidth',1,'LineStyle',':');
legend([p61,p62,p63,p64],'Incidence','Prevalence','INC','PRE')
ylabel('$\mathcal{A}^{\mathrm{max},\ell^1}_t$','interpreter','latex')
xlim([Time(1) Time(end)])
ylim([low_lim_A up_lim_A])
box off
ax = gca;
ax.YAxis(1).Color = color_list_2(2);
ax.YAxis(2).Color = color_list_2(3);

subplot(3,1,3)
yyaxis left
hold on
plot(Time,E2.E2,'Color',color_list_2(4),'LineWidth',1)
ylabel('$\mathcal{E}^{\ell^2}_t$','interpreter','latex')
plot(Time,ones(1,length(Time)),'--r','LineWidth',0.5)
Ypatch = repmat([low_lim_A low_lim_A up_lim_A up_lim_A]', 1,size(Xpatch,2));
p=patch(Xpatch,Ypatch,'k');
p.FaceAlpha=.2;p.EdgeAlpha=0;
xlim([Time(1) Time(end)])
ylim([low_lim_A up_lim_A])

yyaxis right
hold on
plot(Time,Amax(2,:),'Color',color_list_2(5),'LineWidth',1)
ylabel('$\mathcal{A}^{\mathrm{max},\ell^2}_t$','interpreter','latex')
xlim([Time(1) Time(end)])
ylim([low_lim_A up_lim_A])
box off
ax = gca;
ax.YAxis(1).Color = color_list_2(4);
ax.YAxis(2).Color = color_list_2(5);

%% CORRELATIVE ANALYSIS
R = R2.Q50;


CORR_COEFS_RG = zeros(20,1);
CORR_COEFS_E2 = zeros(20,1);
DOMINANT_RL = zeros(20,1); DOMINANT_RR = zeros(20,1);
subindex = @(A, r, c) A(r, c);     % An anonymous function for 2-D indexing

%compute rho
rho = zeros(20,600);
for i = 1:size(Z,3)
    Zt = squeeze(Z(:,:,i));
    temp = Zt./repmat(ResPop,1,20);
    rho(:,i) = sum(temp,1)'.*ResPop;
end

for i = 1:20
    CORR_COEFS_E2(i)=subindex(corrcoef(R2.Q50(i,~isnan(R(i,:)))',E2.E2(~isnan(E2.E2))'),1,2);
    CORR_COEFS_E1(i)=subindex(corrcoef(R2.Q50(i,~isnan(R(i,:)))',E2.E1(~isnan(E2.E1))'),1,2);
    CORR_COEFS_RG(i)=subindex(corrcoef(R2.Q50(i,~isnan(R(i,:)))',RG2(~isnan(RG2))'),1,2);
    DOMINANT_RL(i) = sum(R(i,:)==max(R))/sum(~isnan(max(R)))*100;
    DOMINANT_RR(i) = sum((R(i,:).*rho(i,:))==max(R.*rho))/sum(~isnan(max(R)))*100;
end
CORR_COEFS_ER = subindex(corrcoef(RG2(~isnan(RG2))',E2.E2(~isnan(E2.E2))'),1,2);
CORR_COEFS_EM = subindex(corrcoef(max(R(:,~isnan(RG2)))',E2.E2(~isnan(E2.E2))'),1,2);

CORR_COEFS_E1R = subindex(corrcoef(RG2(~isnan(RG2))',E2.E1(~isnan(E2.E1))'),1,2);
CORR_COEFS_E1M = subindex(corrcoef(max(R(:,~isnan(RG2)))',E2.E1(~isnan(E2.E1))'),1,2);

figure;
subplot(3,1,1)
plot(1:20,CORR_COEFS_RG,'.','MarkerSize',20)
ylim([0.8 1])
xlim([0 21])
set(gca,'XTick',1:20)
set(gca,'xticklabels',[])
ylabel('$\mathrm{correl}(\mathcal{R}_j, \mathcal{R}^G)$','Interpreter','latex')
subplot(3,1,2)
plot(1:20,CORR_COEFS_E2,'.','MarkerSize',20)
ylim([0.8 1])
xlim([0 21])
set(gca,'XTick',1:20)
set(gca,'xticklabels',[])
ylabel('$\mathrm{correl}(\mathcal{R}_j, \mathcal{E}^{\ell^2})$','Interpreter','latex')
subplot(3,1,3)
bar(1:20,DOMINANT_RL*100)
ylabel('$\mathrm{Dominant} \mathcal{R}_j \ [\%]$','Interpreter','latex')
xlim([0 21])
set(gca,'XTick',1:20)
xticklabels(regions)

%% SENSITIVITY ANALYSIS
RGv = zeros(20,2);
E1v = zeros(20,2);
E2v = zeros(20,2);

E1_bsl = E3.E1;
E2_bsl = E2.E2;
RG_bsl = RG2;
for i = 1:20
    R(i,:) = R(i,:)*1.2;
    [L,K] = compute_leslie_matrix(phi,p,R,Z,1);

    temp = compute_epidemicity(L,X2,X2);
    E1_tmp = temp.E1;
    E2_tmp = temp.E2;
    RG_tmp = compute_global_RN(K);

    E1v(i,1) = -sum((E1_bsl-E1_tmp)./E1_bsl,'omitmissing')/591*100;
    E2v(i,1) = -sum((E2_bsl-E2_tmp)./E2_bsl,'omitmissing')/591*100;
    RGv(i,1) = -sum((RG_bsl-RG_tmp)./RG_bsl,'omitmissing')/591*100;

    R(i,:) = R(i,:)/1.2*0.8;
    [L,K] = compute_leslie_matrix(phi,p,R,Z,1);
    
    temp = compute_epidemicity(L,X2,X2);
    E1_tmp = temp.E1;
    E2_tmp = temp.E2;
    RG_tmp = compute_global_RN(K);

    E1v(i,2) = -sum((E1_bsl-E1_tmp)./E1_bsl,'omitmissing')/591*100;
    E2v(i,2) = -sum((E2_bsl-E2_tmp)./E2_bsl,'omitmissing')/591*100;
    RGv(i,2) = -sum((RG_bsl-RG_tmp)./RG_bsl,'omitmissing')/591*100;

    R(i,:) = R(i,:)/0.8;
end

DIFF = abs(RGv(:,1)-RGv(:,2));
[~,IDX] = sort(DIFF,'descend');
E1v_TEMP = E1v;
E2v_TEMP = E2v;
RGv_TEMP = RGv;
fields = regions;


for rk = 1:length(DIFF)
        E1v(rk,:) = E1v_TEMP(IDX(rk),:);
        E2v(rk,:) = E2v_TEMP(IDX(rk),:);
        RGv(rk,:) = RGv_TEMP(IDX(rk),:);

        fields(rk) = regions(IDX(rk));
end


fig=figure;
tiledlayout(1,3)

nexttile
    hold on    
    d = barh(1:length(fields),RGv(:,2));
    d.FaceColor = 'blue'; d.EdgeColor = 'none';
    u = barh(1:length(fields),squeeze(RGv(:,1)));
    u.FaceColor = 'red'; u.EdgeColor = 'none'; 
    xlabel('$\delta \mathcal{R}^{G} \ [\%]$','Interpreter','latex')
    box off
    ax = gca;
    ax.YColor = 'w';
    ax.YAxis.Label.Color='k';
    yticks(1:length(fields))
    yticklabels(fields)
    set(gca,'FontSize',8)
    set(gca,'YDir','reverse')
    set(gca,'FontSize',9)
    xlim([-max(RGv(:,1),[],'all') max(RGv(:,1),[],'all')])

nexttile
    hold on    
    d = barh(1:length(fields),E1v(:,2));
    d.FaceColor = 'blue'; d.EdgeColor = 'none';
    u = barh(1:length(fields),squeeze(E1v(:,1)));
    u.FaceColor = 'red'; u.EdgeColor = 'none'; 
    xlabel('$\delta \mathcal{E}^{\ell^1} \ [\%]$','Interpreter','latex')
    box off
    ax = gca;
    ax.YColor = 'w';
    ax.YAxis.Label.Color='k';
    yticks(1:length(fields))
    yticklabels(fields)
    set(gca,'FontSize',8)
    set(gca,'YDir','reverse')
    set(gca,'FontSize',9)
    xlim([-max(E1v(:,1),[],'all') max(E1v(:,1),[],'all')])

nexttile
    hold on    
    d = barh(1:length(fields),E2v(:,2));
    d.FaceColor = 'blue'; d.EdgeColor = 'none';
    u = barh(1:length(fields),squeeze(E2v(:,1)));
    u.FaceColor = 'red'; u.EdgeColor = 'none'; 
    xlabel('$\delta \mathcal{E}^{\ell^2} \ [\%]$','Interpreter','latex')
    box off
    ax = gca;
    ax.YColor = 'w';
    ax.YAxis.Label.Color='k';
    yticks(1:length(fields))
    yticklabels(fields)
    set(gca,'FontSize',8)
    set(gca,'YDir','reverse')
    set(gca,'FontSize',9)
    xlim([-max(E2v(:,1),[],'all') max(E2v(:,1),[],'all')])
