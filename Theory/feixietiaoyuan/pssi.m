function s=pssi(i,j)
%��i����ĵ�j�����������Ӧ��ȫ��λ��
s=rem(i+j,4);
if s==0
    s=4;
end