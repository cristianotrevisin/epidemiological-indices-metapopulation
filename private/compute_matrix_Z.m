function Z = compute_matrix_Z(ResPop,Q,csi)
    Z = zeros(size(Q,1),size(Q,2),size(csi,2));
    for t = 1:size(csi,2)
        x=csi(:,t);
        C = diag(1-x)+Q*diag(x);
        ActPop=C*ResPop;
        P = C.*ResPop'./ActPop;
        Z(:,:,t) = P'*C;
    end
end