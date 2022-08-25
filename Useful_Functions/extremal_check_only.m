function [] = extremal_check_only(Assemblage,CompletedInputs,ZeroTol)

%Function: This is our full recursive function: it will tell us whether 
%an assemblage is recursive.

%%%%%%%%%
%Inputs:


% assemblage:   A (#Outputs, #Inputs, Dim, Dim) array representing a
% steering assemblage;



%CompletedInputs:  A time-saving parameter. Any inputs ABOVE this value of
%considered already extremal. If 0, only non-zero perturbations are looked
%for. 

%ZeroTol : A tolerance, below which an eigenvalue or perturbation distance
% is to be considered ==0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First we find the required parameters: #Inputs, #Outputs , Dimension.
AssemblageSize=size(Assemblage);
NumberOfInputs=AssemblageSize(2);
NumberOfOutputs=AssemblageSize(1);
Dimension=AssemblageSize(3);


%If dealing with non-zero perturbations, we need full information!
if CompletedInputs==0 

%We will have to keep track of all the operators we need;
OperatorLists=cell(NumberOfInputs);

%And also which input/output they are from.
OperatorAssociations=cell(NumberOfInputs);

%We will also need to keep track of all the perturbations
PerturbationLists=cell(NumberOfInputs);

%For each input, we then calculate the following things:

%OperatorList: The operators which we need to check Linear Dependence

%OperatorTrack: an array stating which input and measurement each operator
%corresponds to

%PerturbationLists: the set of possible perturbations for that measurement.
for Input=1:NumberOfInputs
    [OperatorList,OperatorTrack]=find_operators(Assemblage,Input,ZeroTol);
    OperatorLists{Input}=OperatorList;
    OperatorAssociations{Input}=OperatorTrack;
end

%We now need to generate all the possible non-zero perturbations.
PerturbationList=find_perturbations_marginal(OperatorLists,ZeroTol);
MatrixSize=size(PerturbationList);
NumberofPerturbations=MatrixSize(2);

if NumberofPerturbations==0
    disp('Assemblage is extremal')
    return
else
    disp('Assemblage is not extremal')
   return 
end
                                     
                                     
    
    
else

%We start looking at the current input.
CurrentInput=CompletedInputs;

[OperatorList,OperatorTrack]=find_operators(Assemblage,CurrentInput,ZeroTol);
PerturbationMatrix =find_perturbations(OperatorList,ZeroTol);

%If there are no possible perturbations, then we need to move on to the
%next input.
MatrixSize=size(PerturbationMatrix);
NumberofPerturbations=MatrixSize(2);

if NumberofPerturbations==0
    %In this case we can move on to the next measurement.
    extremal_check_only(Assemblage,CompletedInputs-1,ZeroTol);
    
else
disp('Assemblage is not extremal')
return
end

end




    
    
    
    
    
    
    
    
    
    
    
    
