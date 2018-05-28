cd d:\\codefile\\mldatafile
m=20.1797*1822.39/2%(ʹ��Լ��������ԭ�ӵ�λ�����Ե�������Ϊ1)
eg=35.6/315775%(0.00307eV=296J/mol,ʹ��ԭ�ӵ�λ��)
s=274.9/52.92%(lenard jones ��ƽ����룬5.195)
bs=2^(1/6)*s%(lenard jones ����׵ļ��Ϊ2^(1/6)s=1.122s)
%������LJ�Ƶ���ز���
%hatreeת��ΪJ/mol��*27.211*96485

b=0.3*s;
L=20*s;

nbase=300;
ps=@(x,n,L)sqrt(2/L)*sin(n*pi*(x./L-b/2));
VLJ=@(x)4*eg*((s./x).^12-(s./x).^6);%

hmt=zeros(nbase);
for i=1:nbase
    for j=i:nbase
        hmt(i,j)=(i*pi/L)^2/(2*m)*(i==j)+integral(@(x)VLJ(x).*ps(x,i,L).*ps(x,j,L),b,L+b);
        hmt(j,i)=hmt(i,j);
    end
end
[V,E]=eig(hmt)

cd d:\\codefile\\mldatafile
save hmt_LJ.dat hmt -ASCII
save hmt_LJ_V.dat V -ASCII
save hmt_LJ_V.mat V -ASCII
save hmt_LJ_E.dat E -ASCII

myallplot(6,nbase,V,E,ps,VLJ,L,b,5*bs)