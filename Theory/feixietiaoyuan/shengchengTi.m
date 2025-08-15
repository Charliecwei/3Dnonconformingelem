function E=shengchengTi(alpha,beta)
%四面体体积分
n=size(alpha);
n2=n(1);
n=size(beta);
n1=n(1);
E=zeros(n1,n2);
for i=1:n1
    for j=1:n2
        a=alpha(j,:)+beta(i,:);
        E(i,j)=6*prod(factorial(a))/(factorial(sum(a)+3));
    end
end