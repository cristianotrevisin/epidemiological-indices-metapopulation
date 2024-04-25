function R_global = compute_global_RN(K)
    
    T = size(K,3);
    start = find(K(1,1,:)>=0, 1, 'first');
    R_global = zeros(1,T);          % global eff. reprod. number
    if start > 1
            R_global(1:start-1)=NaN;
    end

    for t = start:T
        Kt = squeeze(K(:,:,t));
        R_global(t) = max(real(eig(Kt)));
    end

end

