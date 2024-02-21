function ES = compute_epidemic_subset(L,N,q,X1,X2)

ES.E1.value = zeros(N,size(L,3));
ES.E1.check = false(N,size(L,3));
ES.E2.value = zeros(N,size(L,3));
ES.E2.check = false(N,size(L,3));

for t = 1:size(L,3)
    Lt = squeeze(L(:,:,t));
    for n = 1:N
        for i = 1:q
            I = zeros(N*q,1);
            I(n+N*(i-1)) = 1;
            Y_old = X1*I;
            Y_new = X1*Lt*I;
            temp1(i) = norm(Y_new,1)/norm(Y_old,1);

            Y_old = X2*I;
            Y_new = X2*Lt*I;
            temp2(i) = norm(Y_new,2)/norm(Y_old,2);
        end
        ES.E1.value(n,t) = max(temp1);
        ES.E2.value(n,t) = max(temp2);
    end
end

ES.E1.check = ES.E1.value >=1;
ES.E2.check = ES.E2.value >=1;

end