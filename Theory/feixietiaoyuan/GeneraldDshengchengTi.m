function E=GeneraldDshengchengTi(alpha,beta)
%一般d维单纯性体积分
n=size(alpha);
n2=n(1);
d = n(2)-1;
n = size(beta);
n1=n(1);
E=zeros(n1,n2);
for i=1:n1
    for j=1:n2
        a=alpha(j,:)+beta(i,:);
        E(i,j)=factorial(d)*prod(factorial(a))/(factorial(sum(a)+d));
    end
end