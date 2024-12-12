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
lim = 965; cases_region = cases_region(:,1:lim); Time = Time(1:lim);
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

% 6. Load variant data
load data/variants_italy.mat
variants(lim+1:end,:) = []; % has the same start as the cases

vIT = table2array(variants);

% Build generation time distribution
mean_GD = [5.20 5.20 5.18 4.43 3.30 2.90 2.80 2.30];
U95 =     [6.00 5.47 5.43 4.49 3.70 3.10 6.70 3.10];
D95 =     [4.40 4.87 4.93 4.36 2.80 2.70 1.50 1.60];


std_GD = (U95-D95)/3.92*sqrt(21);

no_of_var = length(mean_GD);
beta_x = zeros(no_of_var,q);
beta_plot = zeros(no_of_var,q*100+1);
% Survival function (on infectious compartment)
p = @(x) exp(-0.068*x);
p_x = p(0:q-1);       % survival probability    
sigma = exp(-0.068);  % fraction of reaching next day
phi = zeros(no_of_var,21);
for i = 1:no_of_var
    a_beta = (mean_GD(i)/std_GD(i))^2;
    b_beta = mean_GD(i)/a_beta;
    beta = @(x) gampdf(x,a_beta,b_beta);

    beta_x(i,:) = beta(1:q)/sum(beta(1:q));      % generation times
    beta_plot(i,:) = beta(0:0.01:q)/sum(beta(0:0.01:q));  
    % Ratio
    phi(i,:) = beta_x(i,:)./p_x;

end

% computation of limit RT
Rstar = (1-sigma)./max(phi,[],2);

% Estimate incidence per variant
Fv = Ft.*vIT';
Fw = zeros(size(csi,1),no_of_var,lim);

for i = 1:no_of_var
    Fw(:,i,:) = Fp.*vIT(:,i)';
end

%% Compute convolution
alpha = compute_alpha(Fw,beta_x);
alpha_sumV = squeeze(sum(alpha,2));
alpha_sumS = squeeze(sum(alpha,1));
alpha_sumVS = sum(alpha_sumS,1);

% Compute average Phi (variant-unstratified scenario)

phi_avg = zeros(lim,q);
for i = 1:lim
    for j = 1:q
        phi_avg(i,j)=sum(phi(:,j).*alpha_sumS(:,i))/alpha_sumVS(i);
    end
end
phi_avg(1,:) = phi_avg(2,:);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE MATRIX X (OUTPUT);
W1_1V = zeros(20,20*21); % for one strain
W = zeros(20*5,20*21*5); % for all strains


for i = 1:20
    for j = 1:21
        W1_1V(i,i+20*(j-1)) = 1;
    end
end

for i = 1:20
    for j = 1:21
        for v = 1:no_of_var
            W(i+20*(v-1),i+20*(j-1)+20*21*(v-1)) = 1;
        end
    end
end

RG_WT = zeros(200,1); E_WT = zeros(200,1); 
RG_AL = zeros(200,1); E_AL = zeros(200,1);
%%

%A.  Variant unstratified (VUS)
R_VUS = pf(Fp,par,Q,ResPop,csi,alpha_sumV);
for i = 1:lim
    [L_VUSt,K_VUSt] = compute_leslie_matrix(phi_avg(i,:), p_x, R_VUS.Q50(:,i), Z(:,:,i), 1);
    L_VUS(:,:,i) = L_VUSt; K_VUS(:,:,i) = K_VUSt;
end
E_VUS = compute_epidemicity(L_VUS,W1_1V);
RG_VUS = compute_global_RN(K_VUS);

%B. Variant specific (VS)
for i = 1:no_of_var
    i
    par.init = find(vIT(:,i)>0,1,'first')+7;
    Rtemp = pf(squeeze(Fw(:,i,:)),par,Q,ResPop,csi,squeeze(alpha(:,i,:)));
    R_VS(i,:,:) = Rtemp.Q50;
    R_VS(i,:,vIT(:,i)==0) = 0;
end

E_VS.E1 = zeros(1,lim);
RG_VS = zeros(1,lim);

%%
start = find(squeeze(sum(R_VS,[1 2]))>1,1,'first');
for t = start:lim
    t
L_VS = zeros(21*20*no_of_var,21*20*no_of_var);
K_VS = zeros(21*20*no_of_var,21*20*no_of_var);
for i = 1:no_of_var
    if vIT(t,i) >0
        [L_VSt,K_VSt] = compute_leslie_matrix(phi(i,:), p_x, squeeze(R_VS(i,:,t))', Z(:,:,t),1);
    else % if variant is absent, both matrices are null
        L_VSt = zeros(20*21); K_VSt = zeros(20*21);
    end
    if t >=216 && t <= 415
        if i == 1
            E_VStemp = compute_epidemicity(L_VSt,W1_1V,W1_1V);
            E_WT(t-215) = E_VStemp.E1; 
            RG_WT(t-215) = compute_global_RN(K_VSt);
        elseif i == 2
            E_VStemp = compute_epidemicity(L_VSt,W1_1V,W1_1V);
            E_AL(t-215) = E_VStemp.E1; 
            RG_AL(t-215) = compute_global_RN(K_VSt);
        end
    end

    L_VS(1+21*20*(i-1):21*20*i,1+21*20*(i-1):21*20*i) = L_VSt;
    K_VS(1+21*20*(i-1):21*20*i,1+21*20*(i-1):21*20*i) = K_VSt;
    
end
    % Simplify matrix
    sum_col = sum(L_VS,1); sum_row = sum(L_VS,2)';
    idx = find(sum_col == 0 & sum_row == 0);
    L_VS(idx,:) = []; L_VS(:,idx) = [];
    Wtemp = W; Wtemp(:,idx) = []; Wtemp(sum(Wtemp,2)==0,:) = [];
    E_VStemp = compute_epidemicity(L_VS,Wtemp);
    E_VS.E1(t) = E_VStemp.E1;
    RG_VS(t) = compute_global_RN(K_VS);
end


% Make zeros NaNs for easier plotting
RG_VUS(RG_VUS==0) = NaN;
RG_VS(RG_VS==0) = NaN;
RG_WT(RG_WT==0)=NaN; 
RG_AL(RG_AL==0)=NaN;

E_VUS.E1(E_VUS.E1==0) = NaN;
E_VS.E1(E_VS.E1==0) = NaN;
E_WT(E_WT==0)=NaN; 
E_AL(E_AL==0)=NaN;


%% PLOTS BELOW
% Render it nicely on Matlab (vectorial)
set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')

color_variants = ["#277da1","#4d908e","#43aa8b","#90be6d",...
    "#f9c74f","#f8961e","#f3722c","#f94144"];

title_strings = ["Wild Type"; "B.1.1.117 (\alpha)"; "P.1 (\gamma)";...
    "B.1.157.2 (\delta)"; "BA.1 (o)"; "BA.2 (o)";...
    "BA.4 (o)"; "BA.5+BQ.1 (o)"];



% FIGURE 1
f1 = figure;
subplot(2,3,[1 2])
b = bar(Time,vIT'*100,'stacked');
for i = 1:no_of_var
    b(i).FaceColor = color_variants(i);
    b(i).FaceAlpha = .5;
    b(i).EdgeColor = "none";
end
ylim([0 100])
xlim([Time(1) Time(end)])
ylabel("[%]")
subplot(2,3,[4 5])
a = area(Time,Fv');
for i = 1:no_of_var
    a(i).FaceColor = color_variants(i);
    a(i).FaceAlpha = .5;
    a(i).EdgeColor = "none";
end
xlim([Time(1) Time(end)])
ylabel("New cases")
box off

subplot(2,3,6)
hold on
for i = 1:no_of_var
    area(0:0.01:q,beta_plot(i,:), 'FaceColor',color_variants(i),...
        'EdgeColor',color_variants(i),'FaceAlpha',.5)
end
legend(title_strings)
legend("boxoff")
xlabel('Days since infection')
ylabel('\beta')
f1.Units = 'points';
set(findall(gcf,'-property','FontSize'),'FontSize',9)


f1c = figure;
hold on
for i=1:2
    plot(Time,Fv(i,:),'LineWidth',1,'Color',color_variants(i));
end
xlim([Time(216) Time(415)])
box off
f1c.Units = 'points';
set(findall(gcf,'-property','FontSize'),'FontSize',9)


% FIGURE 2

f2 = figure;
tiledlayout(5*2,2*2);

for i = 1:no_of_var
    nexttile([2 2])
    plot(Time, squeeze(R_VS(i,:,:))','linewidth',.4,'color',color_variants(i));
    hold on
    plot(Time, ones(length(Time),1),'--r')
    xlim([Time(1) Time(end)])
    limsy=get(gca,'YLim');
    ylim([0 limsy(2)])
    title(title_strings(i))
    if mod(i,2)~=0
        ylabel('$\mathcal{R}^v_j(t)$','interpreter','latex') 
    end
end
nexttile([2 4])
plot(Time, R_VUS.Q50,'linewidth',.4,'color','#353535');
hold on
plot(Time, ones(length(Time),1),'--r')
xlim([Time(1) Time(end)])
limsy=get(gca,'YLim');
ylim([0 limsy(2)])
title('Not variant-stratified')
ylabel('$\mathcal{R}_j(t)$','interpreter','latex')
set(findall(gcf,'-property','FontSize'),'FontSize',9)

% FIGURE 3

color_indices = ["#219ebc","#023047","#ffb703","#fb8500"];
f3 = figure;
subplot(2*2,3*2,[1 2 3 4 7 8 9 10])
hold on
plot(Time,RG_VUS,'linewidth',1,'Color',color_indices(2))
plot(Time,RG_VS,'linewidth',1,'Color',color_indices(4))
plot(Time,ones(length(Time),1),'--r')
xlim([Time(1) Time(end)])
ylim([.5 3.5])
ylabel("$\mathcal{R}^G(t)$",'Interpreter','latex')
line([Time(216) Time(216)], [.5 3.5], 'color','#353535','linestyle','--')
line([Time(415) Time(415)], [.5 3.5], 'color','#353535','linestyle','--')
for i = 1:8
    a = find(vIT(:,i)>0,1,'first');
    plot(Time(a),.5,'.','Marker','diamond','MarkerSize',10,'color',color_variants(i),'MarkerFaceColor',color_variants(i));
end
subplot(2*2,3*2,[13 14 15 16 19 20 21 22])
hold on
plot(Time,E_VUS.E1,'linewidth',1,'Color',color_indices(2))
plot(Time,E_VS.E1,'linewidth',1,'Color',color_indices(4))
plot(Time,ones(length(Time),1),'--r')
xlim([Time(1) Time(end)])
ylim([.5 4])
line([Time(216) Time(216)], [.5 4], 'color','#353535','linestyle','--')
line([Time(415) Time(415)], [.5 4], 'color','#353535','linestyle','--')
ylabel("$\mathcal{E}(t)$",'Interpreter','latex')
for i = 1:8
    a = find(vIT(:,i)>0,1,'first');
    plot(Time(a),.5,'.','Marker','diamond','MarkerSize',10,'color',color_variants(i),'MarkerFaceColor',color_variants(i));
end

subplot(2*2,3*2,[5 6 11 12])
hold on
plot(Time(216:415), RG_WT,'linewidth',1,'Color', color_variants(1))
plot(Time(216:415), RG_AL,'linewidth',1,'Color', color_variants(6))
plot(Time,RG_VUS,'linewidth',1,'Color','#353535','linestyle','--')
plot(Time,ones(length(Time),1),'--r')
xlim([Time(216) Time(415)])
ylim([.5 3.5])
set(gca,'yticklabels',[])

subplot(2*2,3*2,[17 18 23 24])
hold on
plot(Time(216:415), E_WT,'linewidth',1,'Color', color_variants(1))
plot(Time(216:415), E_AL,'linewidth',1,'Color', color_variants(6))
plot(Time,E_VUS.E1,'linewidth',1,'Color','#353535','linestyle','--')
plot(Time,ones(length(Time),1),'--r')
xlim([Time(216) Time(415)])
ylim([.5 4])
set(gca,'yticklabels',[])

set(findall(gcf,'-property','FontSize'),'FontSize',9)
