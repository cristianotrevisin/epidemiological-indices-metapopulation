clear all
close all
clc

% Build generation time distribution
mean_GD = 5.2; 
std_GD = 1.7; 
k = 21; 
a_beta = (mean_GD/std_GD)^2;
b_beta = mean_GD/a_beta;
beta = @(x) gampdf(x,a_beta,b_beta);
beta_x = beta(1:k)/sum(beta(1:k));      % generation times
% Survival function
p = @(x) exp(-0.068*x);
p_x = p(0:k-1);                           % survival probability
sigma = p(2:k)./p(1:k-1);               % fraction of reaching next day
phi = beta_x./p_x;
sx = p(1);
sigma0 = 1;

% Tested outgoing mobility
xis = [0.001; 0.002; 0.005; 0.01; 0.02; 0.05; 0.1; 0.2; 0.5];

% Initialise
N_metapopulations = 50;
T=50;
E2con = zeros(N_metapopulations,length(xis),T);
E2dis = zeros(N_metapopulations,length(xis),T);

RGcon = zeros(N_metapopulations,length(xis),T);
RGdis = zeros(N_metapopulations,length(xis),T);

% Distribution for effective reproduction numbers
mu_ln = 1.05;
va_ln = 0.33;

sigma = sqrt(log(va_ln/(mu_ln^2)+1));
mu = log(mu_ln) - sigma^2/2;

%% Simulate

for tt = 1:N_metapopulations
    N = 1+randi(9); % size of the metapopulation (no of connected nodes)
    ranks = randperm(N); 
    ResPop = zipf(ranks, 2, 10000)';
    
    C = rand(N); C = C./sum(C,1);
    x=1-diag(C); % percentage of moving pop
    Q=(C-diag(diag(C)))*diag(1./x); % extradiagonal fluxes

    R = lognrnd(mu,sigma, [N, T]);
    store_N(tt) = N;
    store_R{tt} = R;

    % Get transformation matrix
    X2 = zeros(N,N*21); 

    for i = 1:N
        for j = 1:21
            X2(i,i+N*(j-1)) = 1/ResPop(i);
        end
    end

    % Do computation
    for j = 1:length(xis)
        j
    
        % Connected scenario
        csi = normrnd(xis(j),xis(j),[N, T]);
        csi(csi<0) = -csi(csi<0); csi(csi>1) = 1;

        Z = compute_matrix_Z(ResPop,Q,csi);
        [L,K] = compute_leslie_matrix(phi,p,R,Z,1);
        E2con(tt,j,:) = compute_epidemicity(L,X2,X2).E2;

        % Disconnected scenario
        [L0,K0] = compute_leslie_matrix(phi,p,R,repmat(eye(N),1,1,T),1);
        E2dis(tt,j,:) = compute_epidemicity(L0,X2,X2).E2;
        
        E21 = squeeze(E2con(tt,j,:)); E22 = squeeze(E2dis(tt,j,:));

        diff_E2(tt,j) = sum(abs(E21-E22)./E22/(length(E22)-sum(isnan(E22))),'omitnan');

      
    end
end


%% GENERATE THE FIGURE

color_list = ["#03071e";
    "#370617";
    "#6a040f";
    "#9d0208";
    "#d00000";
    "#dc2f02";
    "#e85d04";
    "#f48c06";
    "#faa307";
    "#ffba08"];


bounds_x = [0.9*prctile(E2dis,1,'all') 1.1*prctile(E2dis,99,'all')];
bounds_y = [0.9*prctile(E2con,1,'all') 1.1*prctile(E2con,99,'all')];

figure()
iX = 0;
for i1 = 1:3
    for i2 = 1:3
        iX = iX+1;
        pos = [1+2*(i2-1)+12*(i1-1):2+2*(i2-1)+12*(i1-1)...
            7+2*(i2-1)+12*(i1-1):8+2*(i2-1)+12*(i1-1)];
        subplot(2*4,2*3,pos)
        if i1 < 3; set(gca,'XTickLabel',[]); end
        if i2 > 1; set(gca,'YTickLabel',[]); end
        if i1 == 3 && i2 == 2; xlabel('$E^{\ell^2}$','interpreter','latex'); end
        if i1 == 2 && i2 == 1; ylabel('$\mathcal{E}^{\ell^2}$','interpreter','latex'); end
        text(bounds_x(1)/0.9, bounds_y(2)/1.1,"\xi = " + xis(iX));
        xlim(bounds_x); ylim(bounds_y)
        hold on
        for iN = min(store_N):max(store_N)
            idxs = find(store_N) == iN;
            scatter(reshape(E2dis(idxs,iX,:),1,[]),reshape(E2con(idxs,iX,:),1,[]),...
                'filled','MarkerFaceColor',color_list(iX),...
                'MarkerEdgeColor',color_list(iX),...
                'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0)
        end
        plot(0:5,0:5,'--k')

        
    end
end


for i = min(store_N):max(store_N)
    idxs = find(store_N) == i;
    err_mean_nodecount(i) = mean(diff_E2(idxs,:));
end
meandiff = mean(diff_E2,1);

subplot(2*4,2*3,[37:39 43:45])
hold on
for i = 1:9
    plot(xis(i),meandiff(i),'.','MarkerSize',20,'Color', color_list(i));
end
set(gca,'XScale','log','YScale','log')

xlabel('\xi')

ylabel('MAPE')
box off


subplot(2*4,2*3,[40:42 46:48])
hold on
for i = 1:9
    scatter(store_N,diff_E2(:,i),'filled','MarkerFaceColor', color_list(i))
end
xlim([min(store_N) max(store_N)])
xticks(unique(store_N))
xticklabels(unique(store_N))

ylabel('MAPE')
xlabel('N')
set(gca,'YAxisLocation','right')


function P = zipf(rank, expn, minP)


    H = sum(1./(rank.^expn));
    
    
    P = 1./rank.^expn./H;
    P = P/min(P)*minP;
    
end