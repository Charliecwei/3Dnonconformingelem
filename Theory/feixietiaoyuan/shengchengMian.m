function E=shengchengMian(alpha,beta)
%面积分
n2=size(alpha,1);
no=size(beta,1);
n1=4*no;
E=zeros(n1,n2);
as=zeros(1,3);
for i=1:4
    %io=4+1-i;
    io = i; %%2019.5.28 面1，2，3，4顺序
    for j=1:no
           g=(i-1)*no+j;
        for r=1:n2
               a=alpha(r,:);
               if DD(a(io))==0
                    E(g,r)=0;
               else
                   for k=1:4
                         s=pisju(io,k);
                      if s~=0
                         as(s)=a(k);
                      end
                   end
                   as=as+beta(j,:);
                   E(g,r)=2*prod(factorial(as))/(factorial(sum(as)+2));
               end
        end
    end
end