function[Operator]=VectoOp(Vector)
dim=length(Vector);
MDim=sqrt(dim);
Operator=zeros(MDim,MDim);
vectorco=1;
for row=1:MDim
for col=row:MDim
if row==col
Operator(row,col)=Vector(vectorco);
vectorco=vectorco+1;
else
Operator(row,col)=1/sqrt(2)*(Vector(vectorco)-1j*Vector(vectorco+1));
Operator(col,row)=1/sqrt(2)*(Vector(vectorco)+1j*Vector(vectorco+1));
vectorco=vectorco+2;
end
end
end