function [R,T] = VPPnP(Pts,impts,initf)

if nargin<3
    initf=@initializer;
end
[R0,t0]=initf(Pts,impts);
npt=size(impts,2);
RTt0=[R0.'*t0];  %RTt0 is t in world frame
weights=ones(1,npt);
impts_ones=[impts ;ones(1,npt)];
v=normc(impts_ones);
data=[Pts; v;weights];

RTt = LM(@(x,isNoJ) opt_fun_pure_R_dist_orig(x,isNoJ,data),RTt0);
[y,~,Ropt] = opt_fun_pure_R_dist_orig(RTt,1,data) ;
T=Ropt.'*RTt(1:3);
R=Ropt.';
end
%%%%%%%%%%%%%%%%
function [y,J,Ropt] = opt_fun_pure_R_dist_orig(RTt,isNoJ,data) 
Pts=data(1:3,:);      v=data(4:6,:);    Ws=data(7,:);
npt=size(Pts,2);
% dR=angle2dcm(x(1),x(2),x(3));
Xc=Pts+repmat(RTt,1,npt);
% nc=normc(Xc);
if isNoJ
    [nc]=d_normal(Xc,isNoJ);
else
    [nc,Jnc_Xc]=d_normal(Xc);
end
Ws = Ws./sum(Ws);
vw=v.*repmat(Ws,[3 1]);
B = nc*vw.';          % Eq13b Humble Problems
if isNoJ
    [U, S, V] = dsvd(B,isNoJ);
else
    [U, S, V,JU_B,JS_B,JV_B] = dsvd(B);
end

Sig=[1 0 0;0 1 0;0 0 det(U*V.')];
Ropt = U*Sig*V.';

res=(nc-Ropt*v);
y=res(:);
if isNoJ
    J=[];
    return;
end
I3=eye(3);
JXc=repmat(I3,[1 1 npt]);
JB_nc=repmat(I3,[3 1 npt]).*repelem(permute(vw,[1 3 2]),3, 3);
Jnc=tmult(Jnc_Xc,JXc);
JB=squeeze(sum(tmult(JB_nc,Jnc),3));
JU=JU_B*JB;
JV=JV_B*JB;
JR_U=kron(V*Sig,I3);
JRt_V=kron(U*Sig,I3);
rearrange=vec(reshape(1:numel(Ropt),size(Ropt)).');
JR_V=JRt_V(rearrange,:);
JR=JR_U*JU+JR_V*JV;Jres=reshape(permute(Jnc,[1 3 2]),[],3)-kron(v.',I3)*JR;
J=Jres;
end
%%%%%%%%%%%%%%%%%%%%
function x = LM(FUN,x0)
x = x0;
epsx = 1e-5;
lambda=0.001;
dx=0;
isAccepted=true;
cnt=0;
max_iter=10;
while norm(dx)>epsx&&0<cnt&&cnt<max_iter ||cnt==0
    if isAccepted
        [r,J] = FUN(x,0);
        Jt=J.';
        JtJ=Jt*J;
        Jtr=Jt*r;
        D=eye(size(JtJ)).*(JtJ);
    end
    
    H=JtJ+lambda*D;
    dx=H\Jtr;
    
    r_new=FUN(x-dx,1);
    if norm(r_new)<norm(r)
        isAccepted=true;
        lambda=lambda/10;
        x=x-dx;
    else
        isAccepted=false;
        lambda=lambda*10;
    end
    cnt=cnt+1;
end
end