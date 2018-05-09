function [R,t]= DLT(XXw0,xx)
n= size(xx,2);
mXXw=mean(XXw0,2);
XXw=XXw0-repmat(mXXw,1,n);
% [U,S,V]=svd(XXw,'econ');
% if S(3,3)/S(2,2)<1e-3
%     S(3,3)=0;
%     XXw=U*S*V.'; % dont multiply U, do it at the end
% end
% D= zeros(n*2,12);
% for i= 1:n
%     xi= XXw(1,i); yi= XXw(2,i); zi= XXw(3,i);
%     ui= xx(1,i); vi= xx(2,i);
%     D_= [xi yi zi 0 0 0 -ui*xi -ui*yi -ui*zi 1 0 -ui;
%         0 0 0 xi yi zi -vi*xi -vi*yi -vi*zi 0 1 -vi];
%     D(i*2-1:i*2,:)= D_;
% end
xi= XXw(1,:).'; yi= XXw(2,:).'; zi= XXw(3,:).';
ui= xx(1,:).'; vi= xx(2,:).';
z0=zeros(n,1);
o1=ones(n,1);
D= [xi yi zi z0 z0 z0 -ui.*xi -ui.*yi -ui.*zi o1 z0 -ui;
    z0 z0 z0 xi yi zi -vi.*xi -vi.*yi -vi.*zi z0 o1 -vi];
% D=[D;null(D).'];
DD= D.'*D;
[V,~]= eig(DD);

v= V(:,1); 
v= v/sqrt(norm(v(1:3))*norm(v(4:6)));
% v= v*sign(v(12));
% v(7:9)=cross(v(1:3),v(4:6));

R= reshape(v(1:9),3,3).';
detR=det(R);
detR=sign(detR);%*abs(detR)^(-1/3);
R=R*detR;
t=detR* v(10:12);
t=t-R*mXXw;

XXc= R*XXw+repmat(t,1,size(XXw,2));
[R,t]= calcampose(XXc,XXw);
return

function [R2,t2] = calcampose(XXc,XXw)

n= length(XXc);

% K= eye(n)-ones(n,n)/n;
% mX=mean(X,2)
% ux= mean(X,2);
uy= mean(XXc,2);
% sigmx2= mean(sum((X*K).^2));
% SXY= Y*K*(X')/n;

% sigmx2= mean(sum((X).^2));
XXc=(XXc-repmat(uy,1,n));
SXY= XXc*XXw.';
[U, ~, V]= svd(SXY);
% S= eye(3);
% if det(U*V.') < 0
%     S(3,3)= -1;
% end

R2= U*[1 0 0;0 1 0;0 0 det(U*V.')]*V.';
% c2= trace(D*S)/sigmx2;
% t2= uy-c2*R2*ux;
t2=uy;

% X= R2(:,1);
% Y= R2(:,2);
% Z= R2(:,3);
% if norm(cross(X,Y)-Z) > 2e-2
%     R2(:,3)= -Z;
% end

return
