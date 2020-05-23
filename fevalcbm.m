function out=fevalcbm(cbm,ut)
% ut is a column vector
[s1,s2]=size(cbm);
out1=cellfun(@(c) feval(c,ut),cbm);
out1=cell2mat(out1.blocks);
len=length(ut);
out=zeros(s1,s2,len);
for k=1:len
    out(:,:,k)=out1(k:len:(s1-1)*len+k,:);
end
end