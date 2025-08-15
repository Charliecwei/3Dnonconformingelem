function E=GeneraldDshengchengMian(alpha,beta)
%面积分
n2=size(alpha,1);
no=size(beta,1);
d = size(alpha,2);%空间维数+1
n1=d*no;
E=zeros(n1,n2);
as=zeros(1,d-1);
for i=1:d  %%面1，2,...,d-1,d顺序
    for j=1:no
           g=(i-1)*no+j;
        for r=1:n2
               a=alpha(r,:);
               if a(i)     % 表示 lambda_i的系数不为0
                    E(g,r)=0;
               else
                   for k=1:d
                         s=Generalpisju(i,k,d);
                      if s~=0
                         as(s)=a(k);
                      end
                   end
                   as=as+beta(j,:);
                   E(g,r)=factorial(d-2)*prod(factorial(as))/(factorial(sum(as)+d-2));
               end
        end
    end
end



function s=Generalpisju(i,k,d)
%第k个全局坐标在第i个面的位置
if k<i
    s=d+k-i;
else
   s=k-i;
end