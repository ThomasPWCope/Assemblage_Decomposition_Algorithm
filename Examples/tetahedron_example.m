%This is example 2 in the paper.

tetroid=zeros(4,2,2,2);
QuantumState=1/2*[[1,0,0,1];[0,0,0,0];[0,0,0,0];[1,0,0,1]];
tetroid(1,1,:,:)=PartialTrace(kron([[1,1];[1,1]]/2,eye(2))*QuantumState,1);
tetroid(2,1,:,:)=PartialTrace(kron([[1,-1];[-1,1]]/2,eye(2))*QuantumState,1);
POVM1=1/4*(eye(2)+(sqrt(8/9))*Pauli(1)-(1/3)*Pauli(3));
tetroid(1,2,:,:)=PartialTrace(kron(POVM1,eye(2))*QuantumState,1);
POVM2=1/4*(eye(2)+(-sqrt(2/9))*Pauli(1)+(sqrt(2/3))*Pauli(2)-1/3*Pauli(3));
tetroid(2,2,:,:)=PartialTrace(kron(POVM2,eye(2))*QuantumState,1);
POVM3=1/4*(eye(2)+(-sqrt(2/9))*Pauli(1)-(sqrt(2/3))*Pauli(2)-1/3*Pauli(3));
tetroid(3,2,:,:)=PartialTrace(kron(POVM3,eye(2))*QuantumState,1);
POVM4=1/4*(eye(2)+1*Pauli(3));
tetroid(4,2,:,:)=PartialTrace(kron(POVM4,eye(2))*QuantumState,1);

OutListtet=zero_recursion_H([],tetroid,1,2,10^-6);
congregate_extremals(OutListtet,10^-4)


