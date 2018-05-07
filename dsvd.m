function [U,S,V,dU,dS,dV]=dsvd(A,isNoJ)
    if nargin<2
        isNoJ=false;
    end
    [m,n]=size(A);
    assert(m==n)
    [U,S,V]=svd(A);
    if isNoJ
        return;
    end
    kVtUt=kron(V.',U.');
    dS=kVtUt((1:n)+(0:n-1)*m,:);
    S1=diag(S);
    s_eps=1e-5;
    if(abs(S1(1)-S1(2))<s_eps||abs(S1(2)-S1(3))<s_eps||abs(S1(3)-S1(1))<s_eps)
        warning('some sigular values are equal')
    end
    
    dk=repmat(S1,[1,m,m,m]);
    dl=repmat(S1.',[m,1,m,m]);
    uik=repmat(permute(U,[2,3,1,4]),[1,m,1,m]);
    uil=repmat(permute(U,[3,2,1,4]),[m,1,1,m]);
    vjk=repmat(permute(V,[2,3,4,1]),[1,m,m,1]);
    vjl=repmat(permute(V,[3,2,4,1]),[m,1,m,1]);
    
    invDD=1./(dl.^2-dk.^2);
    invDD(isinf(invDD))=0;

    dC=(dl.*uik.*vjl+dk.*uil.*vjk).*invDD;
    dD=(dk.*uik.*vjl+dl.*uil.*vjk).*invDD;
    dU=reshape(U*reshape(dC,m,[]),[],m*n);
    dV=reshape(V*reshape(dD,m,[]),[],m*n);

end