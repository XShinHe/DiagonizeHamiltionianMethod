clear;clf;
cd d:\\codefile\\mldatafile
%-------------------------参数------------------------------
L=20;
m=1;w=1;
nbase=50;
%-----------------------------------------------------------
x=-L/2:L/1000:L/2;
ps=@(x,n,L)sqrt(2/L)*sin(n*pi*(x./L+1/2));%hilbert空间基组
En=@(n,L,m)(n*pi/L)^2*1/(2*m);%hilbert空间基组能量
Vp=@(x)0.5*x.^2;%(视势阱底部为能量零点)
%Vp=@(x)0.5*(x.^2-2).^2;

hmt=zeros(nbase);
for i=1:nbase
    for j=i:nbase
        hmt(i,j)=En(i,L,m)*(i==j)+trapz(x,Vp(x).*ps(x,i,L).*ps(x,j,L));
        hmt(j,i)=hmt(i,j);
    end
end
[V,E]=eig(hmt)

cd d:\\codefile\\mldatafile
save hmt_HM.dat hmt -ASCII
save hmt_HM_V.dat V -ASCII
save hmt_HM_V.mat V -ASCII
save hmt_HM_E.dat E -ASCII

mx=20;
hold on
for p=1:mx
   f=0*x;
     for q=1:nbase
         f=f+ps(x,q,L).*V(q,p)
     end
     plot(x,f.^2+E(p,p),'r')%调整基线
     text(-L/2+L/50,E(p,p)+0.1,strcat('n=',num2str(p),' E=',num2str(E(p,p))))
end
title(strcat('谐振子势——前10个拟合波函数示意图（基线以能量为准）'))
plot(x,Vp(x),'b')
ylim([0,E(mx+1,mx+1)])
hold off

F=getframe(gcf)
imwrite(F.cdata,strcat('psof_LJ_tot','.png')) 

clf
for p=1:mx
   f=0*x;
     for q=1:nbase
         f=f+ps(x,q,L).*V(q,p)
     end
     hold on
     plot(x,f,'r')%调整基线
     title(strcat('n=',num2str(p),' E=',num2str(E(p,p))))
     hold off
     
     F=getframe(gcf)
     imwrite(F.cdata,strcat('psof_HarmOs_',num2str(p),'.png'))
     clf
end