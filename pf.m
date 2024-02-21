function [R,diagnostic,model_out] = pf(data,par,Q,ResPop,csi,alpha)
% PF Runs the particle filtering algorithm.
%   Inputs:
%       - data: the epidemiological data;
%       - par: the parameter set;
%       - Q: the initial contact matrix;
%       - ResPop: the resident population;
%       - csi: the fraction of outgoing mobility
%   Outputs:
%       - R: the estimated reproduction number;
%       - diagnostic: the diagnostics.
%       - model_out: the simulated cases according to the model's
%       predictions

    Nc = size(data,1);
    Np = par.Np;
    lik = par.lik; 

    alpha_min=par.alpha_min;  % minimum value of alpha to perform the resampling of parameters
    delta=par.delta;              % perform resampling when Neff<tau*N

    data(data==0) = NaN;
    
    % par.initialisation of vectors

    R.Q50 = zeros(Nc,size(data,2)); R.Q05 = zeros(Nc,size(data,2)); 
    R.Q25 = zeros(Nc,size(data,2)); R.Q75 = zeros(Nc,size(data,2));
    R.Q95 = zeros(Nc,size(data,2));

    diagnostic.loglike = zeros(1,size(data,2)); 
    diagnostic.ESS = zeros(1,size(data,2));
    
    diagnostic.sigma_r.Q50 = zeros(Nc,size(data,2));
    diagnostic.sigma_r.Q05 = zeros(Nc,size(data,2));
    diagnostic.sigma_r.Q95 = zeros(Nc,size(data,2));
    

    %% mean and std for parameter r (Rt, log-normally distributed)
    r_mu_new = 3*ones(Nc,1); 

    cv_r_0=par.cv_r_0;
    sigma_r_new=cv_r_0*r_mu_new;
    low_cv_r=par.low_cv_r;
    
    %% Initialize weights
    logW_old = -log(par.Np)*ones(1, par.Np);   % initial log weights (w=1/Np)

    %% sample initial candidates for parameters 
        
    cv_r=sigma_r_new./r_mu_new; % coefficient of variation of the lognormal distrib.

    si2_r=log(cv_r.^2+1); % sigma of the normal distribution associated to the lognormal
    mu_r=log(r_mu_new)-0.5*si2_r; % mean of the normal distribution associated to the lognormal

    r_cand = exp(normrnd(repmat(mu_r,1,Np),repmat(sqrt(si2_r),1,Np)));  % samples of the lognormal distribution for R

    % Implementation of the particle filtering
    for t = par.init:size(data,2)
   
        %disp(' ');
        %disp(['time: ',num2str(t),' ; ',num2str(t/size(data,2)*100),' %' ]);

        % Preparing contact matrix
        x=csi(:,t);
        C = diag(1-x)+Q*diag(x);
        ActPop=C*ResPop;

        % Computation of cases and weights
        mu=zeros(size(data,1), par.Np);
        w = zeros(1, par.Np);
        
        for i = 1:par.Np
            if strcmp(lik,'V1')
                mu(:,i) = (C'*((C*(r_cand(:,i).*alpha(:,t)))./ActPop)); % ENRICO/MARINO, PERCENTUALE
            elseif strcmp(lik, 'V2')
                mu(:,i) = (C'*(r_cand(:,i)./ActPop.*(C*alpha(:,t)))); % CRISTIANO, PERCENTUALE
            end
            w(i) = -log(sum((log(data(:,t)./ResPop./mu(:,i))).^2,'omitnan'));
            mu(:,i) = mu(:,i).*ResPop;
        end
       

        % Normalisation
        w = w + logW_old;
        w(isnan(w)) = -Inf;
        w = exp(w-max(max(w)));%.*w_old;
        
        w = w/sum(w);
        diagnostic.ESS(t) = 1/sum(w.^2);

        %disp(['Neff= ',num2str(diagnostic.ESS(t))]);

        if diagnostic.ESS(t) < par.Np*delta
            

            sel = systematic_resampling(w,par.Np);

            r_temp = zeros(size(r_cand));

            for nn = 1:Nc
                if alpha(:,t) > alpha_min
                    r_temp(nn,:) = r_cand(nn,sel);
                else
                    r_temp(nn,:) = r_cand(nn,:);
                end
            end

            r_mu_new= mean(r_temp,2);
            sigma_r_new = std(r_temp,0,2);
            cv_r=sigma_r_new./r_mu_new;
            cv_r(cv_r<low_cv_r)=low_cv_r;
    
            si2_r=log(cv_r.^2+1);         % variance of the normal distribution associated to the lognormal
            mu_r=log(r_mu_new)-0.5*si2_r; % mean of the normal distribution associated to the lognormal
    
    
            % samples of the lognormal distribution for R
            r_cand = exp(normrnd(repmat(mu_r,1,Np),...
                repmat(sqrt(si2_r),1,Np))); 

            % set weights to 1/Np after resampling
            logW_old = -log(par.Np)*ones(1, par.Np);

        end
    
        %To output
        % these statistics should be weighted statistics
        R.Q05(:,t) = prctile(r_cand,5,2);
        R.Q25(:,t) = prctile(r_cand,25,2);
        R.Q50(:,t) = median(r_cand,2);
        R.Q75(:,t) = prctile(r_cand,75,2);
        R.Q95(:,t) = prctile(r_cand,95,2);

        model_out.Q05(:,t) = prctile(mu,5,2);
        model_out.Q25(:,t) = prctile(mu,25,2);
        model_out.Q50(:,t) = median (mu,2);
        model_out.Q75(:,t) = prctile(mu,75,2);
        model_out.Q95(:,t) = prctile(mu,95,2);

        diagnostic.sigma_r.Q50(:,t) = median(sigma_r_new,2);
        diagnostic.sigma_r.Q05(:,t) = prctile(sigma_r_new,5,2);
        diagnostic.sigma_r.Q95(:,t) = prctile(sigma_r_new,95,2);


        %% median statistics
        if strcmp(lik, 'V2')
            mux = (C'*(R.Q50(:,t)./ActPop.*(C*alpha(:,t)))); 
        elseif strcmp(lik, 'V1')
            mux = (C'*((C*(R.Q50(:,t).*alpha(:,t)))./ActPop));
        end

        diagnostic.loglike(t) = -Nc/2*log(sum((log(data(:,t)./ResPop./mux)).^2));

    end

    R.Q50(R.Q50 == 0) = NaN;
end

%% systenatic resampling
function [indn]=systematic_resampling(w,NSample)
u=rand*1/NSample;
csum=0;
j=0;
for cont_sample=1:NSample
    while csum<u
        j=j+1;
        csum=csum+w(j);
        
    end
    indn(cont_sample)=j;
    u=u+1/NSample;
end
return
end