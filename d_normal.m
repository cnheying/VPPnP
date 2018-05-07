function [nr_out,dr]=d_normal(r,isNoJ)
    if nargin<2
        isNoJ=false;
    end
    [s1,s2,s3]=size(r);
    if s3>1
        r=permute(r,[1 3 2]);
    end
    r_norm=sqrt(sum(r.*r,1));
    nr_out=r./repmat(r_norm,[s1 1]);
    if isNoJ
        return;
    end
    nr=permute(nr_out,[1 3 2]);
    num=size(nr,3);
    r3=repmat(nr,[1 3 1]);
    rt3=repmat(permute(nr,[2 1 3]),[3 1 1]);
    I3=repmat(eye(3),[1 1 num]);
    r_norm3=repmat(reshape(r_norm,[1 1 num]),[3 3 1]);
    dr=(I3-r3.*rt3)./r_norm3;