function[ListOfPerturbations] = find_perturbations_Hermitian(OperatorList,ZeroTol)

%This function returns a matrix whose columns are orthogonal perturbations
%of the operators!

%First we read off some information about the list of operators.
OpListDim=size(OperatorList);

%We assign this information to variables:
VectorSpaceDimension= OpListDim(2)^2;
NumberOperators=OpListDim(1);

%We then create the matrix that Ax=0 implies a linear dependency;
A=zeros(VectorSpaceDimension,NumberOperators);

for i=1:NumberOperators
    OperatorList(i,:,:)=1/2*(squeeze(OperatorList(i,:,:))+squeeze(OperatorList(i,:,:))');
    % Each column is an operator;
    A(:,i)=OptoVec(squeeze(OperatorList(i,:,:)));
end

%We then calculate the singular value decomposiiton of A
[U,S,V]=svd(A);

%We then find the number of non-zero elements of S (up to tolerance ZeroTol)
s=sum(sum(S>ZeroTol));
%And use this to get the null basis (by exploiting the descending order). Each of these is a perturbation in our
%Operators space! 
ListOfPerturbations = V(:,s+1:end);





