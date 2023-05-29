function [Dh,Dv]=TVOperatorGen(n)
    Dh=-eye(n^2)+diag(ones(1,n^2-1),1);
    Dh(n:n:n^2,:)=0;
    Dv=-eye(n^2)+diag(ones(1,n^2-n),n);
    Dv(n*(n-1)+1:n^2,:)=0;
end