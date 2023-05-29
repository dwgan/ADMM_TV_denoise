function xp=ADMM_TVdenoise(x,delta,lambda,iteratMax)
    [N,~]=size(x);
    x=reshape(x,N*N,1);
    [Dh,Dv]=TVOperatorGen(N);
    D=sparse([Dh',Dv']');
    d=D*x;
    p=ones(2*N*N,1)/delta;
    IdelDD=inv((eye(N*N)+delta*(D'*D)));
    for ii=1:iteratMax
        u=IdelDD*(x+delta*D'*(d-p));
        d=wthresh(D*u+p,'s',lambda/delta);
        p=p+D*u-d;
    end
    xp=reshape(u,N,N);
end