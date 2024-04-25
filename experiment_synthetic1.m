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

%% EXPERIMENT TO GENERATE FIG. 1 OF THE PAPER

R = [.5;1.5];
N=2;
xi1 = 0:0.01:1;
xi2 = xi1;
N1 = [10000 100000 1000000]; N2 = N1;
C = zeros(2,2);
RG = zeros(length(xi1),length(xi2),length(N1),length(N2));
E2 = zeros(length(xi1),length(xi2),length(N1),length(N2));

for n1 = 1:length(N1)
    for n2 = 1:length(N2)
        for i1 = 1:length(xi1)
            for i2 = 1:length(xi2)

                C(1,1) = 1-xi1(i1);
                C(2,1) = xi1(i1);

                C(1,2) = xi2(i2);
                C(2,2) = 1-xi2(i2);

                RP = [N1(n1);N2(n2)];
                AP = C*RP;
                P = C.*RP'./AP;
                Z = P'*C;

                F = zeros(k*N,k*N);

                T = zeros(k*N,k*N);


                for ii = 1:k
                    F(1:N,1+N*(ii-1):N*ii) = sigma0*phi(ii)*repmat(R',N,1).*Z;
                    if ii < k
                        T(sub2ind(size(F),1+N*ii:N*(ii+1),1+N*(ii-1):N*ii))=sigma(ii);
                    end
                end

            L = F+T;
            NGM = F*inv(eye(k*N)-T);

            try 
                X = eigs(NGM,1,'largestreal'); 
            catch 
                X = NaN; 
            end

            try
                Y = max(svd(L));
            catch
                Y = NaN;
            end

        RG(i1,i2,n1,n2) = X;
        E2(i1,i2,n1,n2) = Y;

            end

        end
    end
end
%% GENERATE THE FIGURE
set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'defaultFigureRenderer', 'painters')
red = [0.8588235294117647 0.16862745098039217 0.2235294117647059];
white = [1 1 1];
blue = [0.1607843137254902 0.2 0.3607843137254902];

CMBWR = flipud([[linspace(red(1),white(1),256);linspace(red(2),white(2),256);...
    linspace(red(3),white(3),256)]' ;
    [linspace(white(1),blue(1),256);linspace(white(2),blue(2),256);...
    linspace(white(3),blue(3),256)]']);


[XX,YY] = meshgrid(xi1,xi2);
figure
for n1 = 1:length(N1)
    for n2 = 1:length(N2)
        subplot(length(N1),length(N2),3*(n1-1)+n2)
        RGv = squeeze(RG(:,:,n1,n2));
        imagesc(xi1,xi2,RGv')
        set(gca,'YDir','normal')
        hold on
        contourf(XX,YY,RGv'<=1,'LevelList',1,'FaceAlpha',0,'EdgeColor','k','linewidth',1)
        if n1==n2
            plot(0:1,flip(0:1),'k','linewidth',1)
        end
        colormap(CMBWR)
        axis equal
        colorbar
        clim([.5 1.5])
        if n2==1
            ylabel('\xi_2')
        end
        if n1==3
            xlabel('\xi_1')
        end
    end
end

%% EXPERIMENT TO GENERATE FIGURE 2 OF THE PAPER

R1 = 0:0.01:1.25; R2 = 0:0.01:1.25;

C = [.75 .25; .25 .75]; RP1 = [1e5; 1e5]; RP2 = [1e5; 5e5]; RP3 = [1e5; 10e5]; RP4 = [1e5; 50e5];
AP1 = C*RP1; AP2 = C*RP2; AP3 = C*RP3; AP4 = C*RP4;
P1 = C.*RP1'./AP1; P2 = C.*RP2'./AP2; P3 = C.*RP3'./AP3; P4 = C.*RP4'./AP4;
Z1 = P1'*C; Z2 = P2'*C;  Z3 = P3'*C; Z4 = P4'*C;

RG2 = zeros(length(R1),length(R2),4);
E11 = zeros(length(R1),length(R2),4);
E12 = zeros(length(R1),length(R2),4);
E22 = zeros(length(R1),length(R2),4);

X1 = zeros(2,2*21);
X2 = zeros(2,2*21);

for i = 1:2
    for j = 1:21
        X1(i,i+2*(j-1)) = 1;
        X21(i,i+2*(j-1)) = 1/RP1(i);
        X22(i,i+2*(j-1)) = 1/RP2(i);
        X23(i,i+2*(j-1)) = 1/RP3(i);
        X24(i,i+2*(j-1)) = 1/RP4(i);
    end
end

for r1 = 1:length(R1);
    for r2 = 1:length(R2);
        Rtemp = [R1(r1);R2(r2)];
        [L1,K1] = compute_leslie_matrix(phi,p_x,Rtemp,Z1,1);
        RG2(r1,r2,1) = compute_global_RN(K1);
        temp = compute_epidemicity(L1,X1,X21);
        E11(r1,r2,1) = temp.E1;
        E22(r1,r2,1) = temp.E2;
        temp = compute_epidemicity(L1,X21,X21);
        E12(r1,r2,1) = temp.E1;
        [L2,K2] = compute_leslie_matrix(phi,p_x,Rtemp,Z2,1);
        RG2(r1,r2,2) = compute_global_RN(K2);
        temp = compute_epidemicity(L2,X1,X22);
        E11(r1,r2,2) = temp.E1;
        E22(r1,r2,2) = temp.E2;
        temp = compute_epidemicity(L2,X22,X22);
        E12(r1,r2,2) = temp.E1;
        [L3,K3] = compute_leslie_matrix(phi,p_x,Rtemp,Z3,1);
        RG2(r1,r2,3) = compute_global_RN(K3);
        temp = compute_epidemicity(L3,X1,X23);
        E11(r1,r2,3) = temp.E1;
        E22(r1,r2,3) = temp.E2;
        temp = compute_epidemicity(L3,X23,X23);
        E12(r1,r2,3) = temp.E1;
        [L4,K4] = compute_leslie_matrix(phi,p_x,Rtemp,Z4,1);
        RG2(r1,r2,4) = compute_global_RN(K4);
        temp = compute_epidemicity(L4,X1,X24);
        E11(r1,r2,4) = temp.E1;
        E22(r1,r2,4) = temp.E2;
        temp = compute_epidemicity(L4,X24,X24);
        E12(r1,r2,4) = temp.E1;
    end
end

%% GENERATE THE FIGURE
[XX,YY] = meshgrid(R1,R2);
figure;
tiledlayout(2,2)
nexttile
hold on
contourf(XX,YY,squeeze(E22(:,:,1))'<=1,'LevelList',1,'FaceAlpha',1,'facecolor','#cce6f4','linestyle','none');
contourf(XX,YY,squeeze(RG2(:,:,1))','LevelList',1,'LineStyle',':','FaceColor','#db2b39');
contourf(XX,YY,squeeze(E11(:,:,1))'<=1,'LevelList',1,'LineStyle','none','FaceColor','#29335c');
contourf(XX,YY,squeeze(E12(:,:,1))'<=1,'LevelList',1,'FaceColor','none','EdgeColor','white','linewidth',1);

xlabel('$\mathcal{R}_1$','interpreter','latex')
ylabel('$\mathcal{R}_2$','interpreter','latex')
title('$n_1 = n_2 = 10^5$','Interpreter','latex')
axis equal
nexttile
hold on
contourf(XX,YY,squeeze(E22(:,:,2))'<=1,'LevelList',1,'FaceAlpha',1,'facecolor','#cce6f4','linestyle','none');
contourf(XX,YY,squeeze(RG2(:,:,2))','LevelList',1,'LineStyle',':','FaceColor','#db2b39');
contourf(XX,YY,squeeze(E11(:,:,2))'<=1,'LevelList',1,'LineStyle','none','FaceColor','#29335c');
contourf(XX,YY,squeeze(E12(:,:,2))'<=1,'LevelList',1,'FaceColor','none','EdgeColor','white','linewidth',1);
xlabel('$\mathcal{R}_1$','interpreter','latex')
title('$n_1 = 10^5, \ n_2 = 5\times 10^5$','Interpreter','latex')
axis equal
nexttile
hold on
contourf(XX,YY,squeeze(E22(:,:,3))'<=1,'LevelList',1,'FaceAlpha',1,'facecolor','#cce6f4','linestyle','none');
contourf(XX,YY,squeeze(RG2(:,:,3))','LevelList',1,'LineStyle',':','FaceColor','#db2b39');
contourf(XX,YY,squeeze(E11(:,:,3))'<=1,'LevelList',1,'LineStyle','none','FaceColor','#29335c');
contourf(XX,YY,squeeze(E12(:,:,3))'<=1,'LevelList',1,'FaceColor','none','EdgeColor','white','linewidth',1);
xlabel('$\mathcal{R}_1$','interpreter','latex')
title('$n_1 = 10^4, \ n_2 = 10^6$','Interpreter','latex')
axis equal
nexttile
hold on
contourf(XX,YY,squeeze(E22(:,:,4))'<=1,'LevelList',1,'FaceAlpha',1,'facecolor','#cce6f4','linestyle','none');
contourf(XX,YY,squeeze(RG2(:,:,4))','LevelList',1,'LineStyle',':','FaceColor','#db2b39');
contourf(XX,YY,squeeze(E11(:,:,4))'<=1,'LevelList',1,'LineStyle','none','FaceColor','#29335c');
contourf(XX,YY,squeeze(E12(:,:,4))'<=1,'LevelList',1,'FaceColor','none','EdgeColor','white','linewidth',1);
xlabel('$\mathcal{R}_1$','interpreter','latex')
title('$n_1 = 10^4, \ n_2 = 5 \times 10^6$','Interpreter','latex')
axis equal