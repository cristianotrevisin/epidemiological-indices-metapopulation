function alpha = compute_alpha(data,beta_x)

    if ndims(data)==2
       if size(data,1)==1 %there is just one node
            data = reshape(data,[1 size(data,1) size(data,2)]);
        else % there is just one variant
            data = reshape(data,[size(data,1) 1 size(data,2)]);
       end
    end

    alpha = zeros(size(data,1),size(beta_x,1),size(data,3));

    q = size(beta_x,2);
    
    for t = 2:size(data,3)
        keff = min(q,t-1);
        for v = 1:size(beta_x,1)
            alpha(:,v,t) = sum(flip(squeeze(data(:,v,t-keff:t-1)),2).*beta_x(v,1:keff),2);
        end
    end
end