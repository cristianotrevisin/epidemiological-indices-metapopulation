function [L,K] = compute_leslie_matrix(phi,p,R,Z,sigma0)

    % Define nodes and length
    N = size(R,1);      %number of nodes
    t_len = size(R,2);      %time length
    q = size(phi,2); %number of phis
    sigma = p(2:q)./p(1:q-1);               % fraction of reaching next day

    % Build Leslie Matrix
    L = zeros(N*q,N*q,t_len);
    K = zeros(N*q,N*q,t_len);

    start = find(sum(R,1)>0, 1, 'first');
    if start > 1
        L(:,:,1:start-1) = NaN;
        K(:,:,1:start-1) = NaN;
    end

    for t = start:t_len
        T = zeros(N*q,N*q);
        S = zeros(N*q,N*q);

        for ii = 1:q
            T(1:N,1+N*(ii-1):N*ii) = sigma0*phi(ii)*repmat(R(:,t)',N,1).*squeeze(Z(:,:,t)); 
            if ii < q
                S(sub2ind(size(T),1+N*ii:N*(ii+1),1+N*(ii-1):N*ii))=sigma(ii);
            end
        end
        K(:,:,t) = squeeze(T)*inv(eye(q*N)-S);
        L(:,:,t) = T+S;
    end 
end