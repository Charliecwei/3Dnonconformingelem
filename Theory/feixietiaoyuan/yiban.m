n=3;
alpha=Duoxiangs3(n);
beta=Duoxiangs2(n-1);
format rat

delta=Duoxiangs3(4);
E=shengchengMian(alpha,beta);
s=0;
k=5;
o0=n2+4-rank(E);
while s<4
    k=k+1;
    alpha1=[alpha;delta(k,:)];
    E=shengchengMian(alpha1,beta);
    o=n2+4-rank(E);
    if o<o0
       s=s+1;
       o0=o;
       alpha=alpha1;
    end
end
alpha(21:24,:);
% E= shengchengMian(delta,beta);
% rank(E)
lambda=[1/4,1/4,1/4,1/4];
G=[2,0,1,1;1,2,0,1;1,1,2,0;0,1,1,2];
alpha=[alpha(1:20,1:4);G];
E=shengchengMian(alpha,beta);
S=shengchengDian(alpha,lambda);
E=[E;S];
rank(E)





ES=shengchengMian(delta,beta);
S=ES(1:6,[12,8,6,14]);
X4=S(:,1)+S(:,2)-S(:,3)-S(:,4);

S=ES(7:12,[7,26,16,13]);
X3=S(:,1)+S(:,2)-S(:,3)-S(:,4);



S=ES(13:18,[19,8,18,9]);
X2=S(:,1)+S(:,2)-S(:,3)-S(:,4);

S=ES(19:24,[28,26,25,29]);
X1=S(:,1)+S(:,2)-S(:,3)-S(:,4);

X=[X1;X2;X3;X4;0];

E=[E,X];
rank(E)


