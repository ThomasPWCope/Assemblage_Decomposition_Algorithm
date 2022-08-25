function[ListOfPerturbations] = find_perturbations_marginal_Hermitian(OperatorLists,ZeroTol)

%This function returns a matrix whose columns are orthogonal perturbations
%of the operators. However, this function focusses on NON-Zero Marginal 
%Perturbations.

%First we read off some information about the number of inputs:
ListSize=size(OperatorLists);
NumberOfInputs=ListSize(2);


%For each input; we need the constraint that the traces add up to 0;
%and that pairwise the no-signalling rule applies.


%We will also keep track of the total number of operators: this will help
%us to properly construct our perturbation matrix.
CurrentNumberOfOperators=0;


% We now range over each Input, constructing our constraints.
for Input = 1:NumberOfInputs
    
    %First, we extract some useful information:
    %List of Operators, How many, and the dimension of the space.
   
    CurrentList=OperatorLists{Input};
    SizeOfList=size(CurrentList);
    NumberOfOperators=SizeOfList(1);
    VectorSpaceDimension=SizeOfList(2)^2;
    
    %We keep track of which operators have trace =1 (though we only use
    %this for our first input constraint.
    NewTraceMatrix=zeros(1,NumberOfOperators);
    
    %We also keep track of the matrix whose columns are operators.
    OperatorMatrix=zeros(VectorSpaceDimension,NumberOfOperators);
    for i=1:NumberOfOperators
            CurrentList(i,:,:)=1/2*(squeeze(CurrentList(i,:,:))+squeeze(CurrentList(i,:,:))');
            %Here we check the trace condition
            if abs(1-trace(squeeze(OperatorLists{Input}(i,:,:))))<=ZeroTol
                NewTraceMatrix(i)=1;
            end
            
            %Here we construct the operator-column matrix.
            OperatorMatrix(:,i)=OptoVec(squeeze(CurrentList(i,:,:)));
    end
        %In the special case of Input ==1, we do not add a no-signalling
        %requirement: but we set up that \sum input=1 - \sum input=Input>1 =0
        if Input==1
            TotalOperatorMatrix=repmat(OperatorMatrix,[NumberOfInputs-1,1]);
            TraceMatrix=NewTraceMatrix;
            NumberFirstOperators=NumberOfOperators;
        
        %For the other inputs, we set up the other half of the
        %no-signalling equality.
        else
            AdditionOperatorMatrix=repmat(0*OperatorMatrix,[NumberOfInputs-1,1]);
            AdditionOperatorMatrix((Input-2)*VectorSpaceDimension+1:(Input-1)*VectorSpaceDimension,1:NumberOfOperators)=-OperatorMatrix;
            TotalOperatorMatrix=[TotalOperatorMatrix,AdditionOperatorMatrix];
            %TotalTraceMatrix = [TotalTraceMatrix,NewTraceMatrix];
            
            
        
        
        end
    %We finally update the total number of operators considered.    
    CurrentNumberOfOperators=CurrentNumberOfOperators+NumberOfOperators;

end

%We need only constrain the trace for the first input: the no-signalling
%will ensure the trace condition for the other inputs.
TotalTraceMatrix=zeros(1,CurrentNumberOfOperators);
TotalTraceMatrix(1:NumberFirstOperators)=TraceMatrix;


%Our full constraint matrix then considers no-sig and trace conditions
%together.
FullPerturbationMatrix=[TotalOperatorMatrix;TotalTraceMatrix];    

%We then calculate the singular value decomposiiton of our Matrix.
[U,S,V]=svd(FullPerturbationMatrix);

%We then find the number of non-zero elements of S (up to tolerance ZeroTol)
s=sum(sum(S>ZeroTol));
%And use this to get the null basis (by exploiting the descending order). Each of these is a perturbation in our
%Operators space! 
ListOfPerturbations = V(:,s+1:end);





