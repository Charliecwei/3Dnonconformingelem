D = zeros(4,4);
D2 = zeros(4,4);

for k = 1:4
    for j = 1:4
        D(k,j)=0.5*((-1)^(k+j))*cot(pi*(k-j)/4);
        D2(k,j) = -0.5*((-1)^(k+j))*(csc(pi*(k-j)/4))^2;
    end
end

D

D2


for i = 1:4
    D(i,i)=0;
end


D*D