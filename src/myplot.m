function []=myplot(mx,nbase,V,E,ps,Vp,L,x1,x2)
clf
x=x1:(x2-x1)/1000:x2;
for p=1:mx
     f=@(x)0;
     for q=1:nbase
         f=@(x)(f(x)+ps(x,q,L).*V(q,p))
     end
     hold on
     plot(x,f(x).^2,'r')%调整基线
     title(strcat('n=',num2str(p),' E=',num2str(E(p,p))))
     xlim([x1,x2])
     hold off
     %text(-0.9*s,E(p,p)+0.000001,strcat('n=',num2str(p),' E=',num2str(E(p,p))))
     F=getframe(gcf)
     imwrite(F.cdata,strcat('Rof_LJ_',num2str(p),'.png'))
     clf
end