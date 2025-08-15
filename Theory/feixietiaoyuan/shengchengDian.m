function E=shengchengDian(alpha,lambda)
%关于在lambda点的取值
n=size(alpha);
n1=n(1);
n=size(lambda);
n2=n(1);
%n3=n(2);
E=zeros(n2,n1);
for i=1:n2
    for j=1:n1
        a=alpha(j,:);
%         s=1;
%         for r=1:n3
%             s=s*power(lambda(i,r),a(r));
%         end
         s=prod(power(lambda(i,:),a));
        E(i,j)=s;
    end
end