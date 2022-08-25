function[vector]=OptoVec(Operator)
dimL=size(Operator);
dim=dimL(1);
vector=zeros(dim^2,1);
vectorco=1;
for row=1:dim
for col=row:dim
if row==col
vector(vectorco)=Operator(row,col);
vectorco=vectorco+1;
else
vector(vectorco)=sqrt(2)*real(Operator(row,col));
vectorco=vectorco+1;
vector(vectorco)=-sqrt(2)*imag(Operator(row,col));
vectorco=vectorco+1;
end
end
end
