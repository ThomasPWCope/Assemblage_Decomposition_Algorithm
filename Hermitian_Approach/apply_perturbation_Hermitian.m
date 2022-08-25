function[Sucess,OutputCell] = apply_perturbation_Hermitian(ProbIn,Assemblage,OperatorList,OperatorAssociation,Perturbation,ZeroYN,ZeroTol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:

%ProbIn: The probability of the current assemblage in the overall
%decompositon.

%Assemblage: The current assemblage in question.

%Operator List: The list of operators being considered under the current
%perturbation.

%Operator Association: Denotes which input/output each Operator belongs to.

%Peturbation: a column vector giving the coefficient of every operator in
%the perturbation.

%ZeroYN: is either non-zero and gives the input to consider,
%or is 0 and is a full assemblage perturbation.

%ZeroTol: Determines whether a perturbation is accepted as valid or not. 


%we need to use the dimensions of the assemblage:
AssemblageSize=size(Assemblage);
NumberOfInputs=AssemblageSize(2);
NumberOfOutputs=AssemblageSize(1);
Dimension=AssemblageSize(3);

%This is what we will add to our assemblage;
PerturbationAddition=0*Assemblage;

%We take the length of the Perturbation;
PerturbLength=length(Perturbation);

%And then create the Perturbation!
for i=1:PerturbLength
    %We need to know where to add the the perturbation on our assemblage.
    Substate=OperatorAssociation(i,:);
    Input=Substate(1);
    Output=Substate(2);
    %and then add it!
    PerturbationAddition(Output,Input,:,:)=squeeze(PerturbationAddition(Output,Input,:,:))+Perturbation(i)*squeeze(OperatorList(i,:,:));
end

%It is also a good idea for stability to re-Hermitian both arrays.
for input=1:NumberOfInputs
    for output=1:NumberOfOutputs
    Assemblage(output,input,:,:)=1/2*(squeeze(Assemblage(output,input,:,:))+squeeze(Assemblage(output,input,:,:))');
    PerturbationAddition(output,input,:,:)=1/2*(squeeze(PerturbationAddition(output,input,:,:))+squeeze(PerturbationAddition(output,input,:,:))');
   end
end
    

cvx_begin sdp quiet
variable d(1,1)
maximise d
subject to
if ZeroYN ~=0
for output=1:NumberOfOutputs
squeeze(Assemblage(output,ZeroYN,:,:) + d * PerturbationAddition(output,ZeroYN,:,:)) >=0;
end
else
for input=1:NumberOfInputs
for output=1:NumberOfOutputs
squeeze(Assemblage(output,input,:,:) + d * PerturbationAddition(output,input,:,:)) >=0;
end
end
end
cvx_end
maxd=cvx_optval;

cvx_begin sdp quiet
variable d(1,1)
minimize d
subject to
if ZeroYN ~=0
for output=1:NumberOfOutputs
squeeze(Assemblage(output,ZeroYN,:,:) + d * PerturbationAddition(output,ZeroYN,:,:)) >=0;
end
else
for input=1:NumberOfInputs
for output=1:NumberOfOutputs
squeeze(Assemblage(output,input,:,:) + d * PerturbationAddition(output,input,:,:)) >=0;
end
end
end
cvx_end
mind=cvx_optval;

%We only want to split if our perturbation is bigger than ZeroTol
if min(abs([maxd,mind])) > ZeroTol
Sucess=1;
Prob1=(-mind/(maxd-mind))*ProbIn;
Prob2=(maxd/(maxd-mind))*ProbIn;
Assemblage1=Assemblage + maxd *PerturbationAddition;
Assemblage2=Assemblage + mind *PerturbationAddition;
OutputCell={Prob1,Assemblage1;Prob2,Assemblage2};
else
Sucess=0;
OutputCell={ProbIn,Assemblage};
end
end
