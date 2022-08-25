function [ListOfExtremalAssemblagesOut] = zero_recursion_H(ListOfExtremalAssemblages,Assemblage,Probability,CompletedInputs,ZeroTol)

%Function: This is our full recursive function: it will be completely
%decompose an assemblage.


%%%%%%%%%
%Inputs:

%List of Extremal Assemblages: Once as assemblage is completely extremal,
%it will be added to this list.

% assemblage:   A (#Outputs, #Inputs, Dim, Dim) array representing a
% steering assemblage;

%Probability: keeps track of the weight of this assemblage in the
%decomposition.

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

%We will Re-Hermitian the assemblage.
%It is also a good idea for stability to re-Hermitian both arrays.
for input=1:NumberOfInputs
    for output=1:NumberOfOutputs
    Assemblage(output,input,:,:)=1/2*(squeeze(Assemblage(output,input,:,:))+squeeze(Assemblage(output,input,:,:))');
   end
end

%For each input, we then calculate the following things:

%OperatorList: The operators which we need to check Linear Dependence

%OperatorTrack: an array stating which input and measurement each operator
%corresponds to

%PerturbationLists: the set of possible perturbations for that measurement.
for Input=1:NumberOfInputs
    [OperatorList,OperatorTrack]=find_operators_Hermitian(Assemblage,Input,ZeroTol);
    OperatorLists{Input}=OperatorList;
    OperatorAssociations{Input}=OperatorTrack;
end

%We now need to generate all the possible non-zero perturbations.
PerturbationList=find_perturbations_marginal_Hermitian(OperatorLists,ZeroTol);
MatrixSize=size(PerturbationList);
NumberofPerturbations=MatrixSize(2);


%I am not sure if we need the things below.
%To do this, we need to combine all our possible operators.
FullOperatorList=vertcat(OperatorLists{:});
FullOperatorAssociations=vertcat(OperatorAssociations{:});

if NumberofPerturbations==0
    %We will have already checked all zero-marginal perturbations, and so
    %we can add our current perturbation to the extremal list!
    ListOfExtremalAssemblages=[ListOfExtremalAssemblages;{Probability,Assemblage}];
    ListOfExtremalAssemblagesOut=ListOfExtremalAssemblages;
else
    %Otherwise, we have to do a perturbation. Like before, we will choose
    %the last perturbation of the matrix.
     
    
    RelevantPerturbation=PerturbationList(:,NumberofPerturbations);
    
    if norm(abs(imag(RelevantPerturbation)))>ZeroTol
        warning('Warning, imaginary perturbation components high')
    end
    
    [Sucess,OutputCell] = apply_perturbation_Hermitian(Probability,Assemblage,FullOperatorList...
                                         ,FullOperatorAssociations,real(RelevantPerturbation),0,ZeroTol);
    
                                     
    if Sucess~=1
    %We should never get a failed perturbation.
    error('Failed Perturbation')
    else
       %If we suceeded, we need to keep going with our perturbations! 
       
        NewAssemblage1=OutputCell{1,2};
        NewProb1=OutputCell{1,1};
        NewAssemblage2=OutputCell{2,2};
        NewProb2=OutputCell{2,1};
        
        ListOfExtremalAssemblagesOut = zero_recursion_H(ListOfExtremalAssemblages,NewAssemblage1,NewProb1,CompletedInputs,ZeroTol);
        ListOfExtremalAssemblages=ListOfExtremalAssemblagesOut;
        ListOfExtremalAssemblagesOut = zero_recursion_H(ListOfExtremalAssemblages,NewAssemblage2,NewProb2,CompletedInputs,ZeroTol);
        ListOfExtremalAssemblages=ListOfExtremalAssemblagesOut;                                
        
    end
end
                                     
                                     
    
    
    



else

%We start looking at the current input.
CurrentInput=CompletedInputs;

[OperatorList,OperatorTrack]=find_operators_Hermitian(Assemblage,CurrentInput,ZeroTol);
PerturbationMatrix =find_perturbations_Hermitian(OperatorList,ZeroTol);

%If there are no possible perturbations, then we need to move on to the
%next input.
MatrixSize=size(PerturbationMatrix);
NumberofPerturbations=MatrixSize(2);

if NumberofPerturbations==0
    %In this case we can move on to the next measurement.
    ListOfExtremalAssemblagesOut=zero_recursion_H(ListOfExtremalAssemblages,Assemblage,Probability,CompletedInputs-1,ZeroTol);
    ListOfExtremalAssemblages=ListOfExtremalAssemblagesOut;
else
    %If there is a valid perturbation, of course we have to apply it! 
    
    %We always will choose the last column as the perturbation we apply.
    RelevantPerturbation=PerturbationMatrix(:,NumberofPerturbations);
    
    if norm(abs(imag(RelevantPerturbation)))>ZeroTol
        warning('Warning, imaginary perturbation components high')
    end
    

    [Sucess,OutputCell] = apply_perturbation_Hermitian(Probability,Assemblage,OperatorList...
                                         ,OperatorTrack,real(RelevantPerturbation),CurrentInput,ZeroTol);
                                     
    if Sucess~=1
        %We should never get a failed perturbation.
        
        %We can first try to clean the assemblage: that is, remove small
        %negative eigenvalues, and make it again Hermitian.
        
        error('Failed Perturbation')
    else
       %If we suceeded, we need to keep going with our perturbations! 

        NewAssemblage1=OutputCell{1,2};
        NewProb1=OutputCell{1,1};
        NewAssemblage2=OutputCell{2,2};
        NewProb2=OutputCell{2,1};
        
        ListOfExtremalAssemblagesOut = zero_recursion_H(ListOfExtremalAssemblages,NewAssemblage1,NewProb1,CompletedInputs,ZeroTol);
        ListOfExtremalAssemblages=ListOfExtremalAssemblagesOut; 
        ListOfExtremalAssemblagesOut = zero_recursion_H(ListOfExtremalAssemblages,NewAssemblage2,NewProb2,CompletedInputs,ZeroTol);
        ListOfExtremalAssemblages=ListOfExtremalAssemblagesOut; 
                                            
        
    end
end

end




    
    
    
    
    
    
    
    
    
    
    
    
