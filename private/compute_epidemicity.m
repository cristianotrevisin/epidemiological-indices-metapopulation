function E = compute_epidemicity(L,X1, X2)
    
    if nargin == 2
        X2 = X1;
    end
    
    T = size(L,3);
    start = find(L(1,1,:)>=0, 1, 'first');
    E.E1 = zeros(1,T);    
    E.E2 = zeros(1,T);
    E.EI = zeros(1,T);
    if start > 1
            E.E1(1:start-1)=NaN;
            E.E2(1:start-1)=NaN;
            E.EI(1:start-1)=NaN;
    end

    for t = start:T
        Lt = squeeze(L(:,:,t));
        temp1 = zeros(size(Lt,1),1);
        temp2 = zeros(size(Lt,1),1);
        temp3 = zeros(size(Lt,1),1);
        for i = 1:size(Lt,1);
            I = zeros(size(L,1),1);
            I(i) = 1;
            Y_old = X1*I;
            Y_new = X1*Lt*I;
            temp1(i) = norm(Y_new,1)/norm(Y_old,1);

            Y_old = X2*I;
            Y_new = X2*Lt*I;
            temp2(i) = norm(Y_new,2)/norm(Y_old,2);
            temp3(i) = norm(Y_new,'inf')/norm(Y_old,'inf');

        end


        E.E1(t) = max(temp1);
        E.E2(t) = max(temp2);
        E.EI(t) = max(temp3);

    end

end

