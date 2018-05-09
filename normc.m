function B=normc(A)
B=A./repmat(sqrt(sum(A.*A,1)),size(A,1),1);