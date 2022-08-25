function [OperatorList,OperatorAssociation] = find_operators(Assemblage,ZeroYN,ZeroTol);

%First we find the required parameters
AssemblageSize=size(Assemblage);
NumberOfInputs=AssemblageSize(2);
NumberOfOutputs=AssemblageSize(1);
Dimension=AssemblageSize(3);

%Let us define our array to store operators;
OperatorList=zeros(0,Dimension,Dimension);
NumberOperatorsFound=0;

%In this case we consider all operators, and remove the 1st non-zero
%diagonal one.
if ZeroYN == 0
error('not written yet')
    
    
    
    
%Otherwise, we only create the operators for 1 measurement.    
else


    for i=1:NumberOfOutputs
        %First we take the current Substate
        CurrentSubstate=squeeze(Assemblage(i,ZeroYN,:,:));
        
        %We calculate the eigenvectors and eigenvalues of the state.
        [EigenvectorMatrix,EigenvalueMatrix]=eigs(CurrentSubstate);
        
        
        
        %Then, we run through the possible operators in it's cokernel.
        for j=1:Dimension
            
            %We will use ZeroTol as the decider of whether something is in
            %the kernel or not.
            if EigenvalueMatrix(j,j)>ZeroTol
                
                %We now need to run through the other eigenvalues; to
                %create off diagonal operators.
                for k=1:j
                    
                    %If k==j we have a trace 1 operator.
                    if k==j 
                        OperatorList(NumberOperatorsFound+1,:,:)= EigenvectorMatrix(:,k)*EigenvectorMatrix(:,k)';
                        OperatorAssociation(NumberOperatorsFound+1,:)=[ZeroYN,i];
                        NumberOperatorsFound=NumberOperatorsFound+1;
                        
                    %Otherwise, we have to add two Hermitian trace 0 operators.    
                    elseif EigenvalueMatrix(k,k)>ZeroTol
                        OperatorList(NumberOperatorsFound+1,:,:)= ...
                      (1/sqrt(2))*(EigenvectorMatrix(:,j)*EigenvectorMatrix(:,k)' + EigenvectorMatrix(:,k)*EigenvectorMatrix(:,j)');
                        OperatorAssociation(NumberOperatorsFound+1,:)=[ZeroYN,i];
                        
                        
                        OperatorList(NumberOperatorsFound+2,:,:)= ...
                      (-1j/sqrt(2))*(EigenvectorMatrix(:,j)*EigenvectorMatrix(:,k)' - EigenvectorMatrix(:,k)*EigenvectorMatrix(:,j)');
                        OperatorAssociation(NumberOperatorsFound+2,:)=[ZeroYN,i];
                        
                        NumberOperatorsFound=NumberOperatorsFound+2;
                    end
                end
            end
        end
    end             
end
end
                        
                        
                        

