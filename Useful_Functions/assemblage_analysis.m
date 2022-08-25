function [] = assemblage_analysis(OriginalAssemblage,Decomposition)

%This function tests whether a decomposition on an assemblage is a good
%one. 

%It performs 3 tests: reproduction of original, no-signalling, and
%positivity.

Dimensions=size(OriginalAssemblage);
NumberOutputs=Dimensions(1);
NumberInputs=Dimensions(2);
Dim=Dimensions(3);

DecompSize=size(Decomposition);
NumberOfExtremals=DecompSize(1);
ReSum=0*OriginalAssemblage;
for i=1:NumberOfExtremals
    ReSum=ReSum+Decomposition{i,1}*Decomposition{i,2};
end
disp('Distance from original')
disp('Total Sum Squares Difference')
OrigVec=reshape(OriginalAssemblage,[NumberOutputs*NumberInputs*Dim*Dim,1]);
ReSumVec=reshape(ReSum,[NumberOutputs*NumberInputs*Dim*Dim,1]);
disp((OrigVec-ReSumVec)'*(OrigVec-ReSumVec));
disp('largest absolute element distance')
disp(max(abs(OrigVec-ReSumVec)));

disp('Positivity')
LargestNeg=0;
AverageLNeg=0;
AverageNeg=0;

MaxNonH=0;
AverageMax=0;
SumSq=0;
WeightedSumSq=0;

for Ex=1:NumberOfExtremals
    CurrentEx=Decomposition{Ex,2};
    LargestLocalNeg=0;
    LargestLocalMaxNonH=0;
    for Input=1:NumberInputs
       
    for  Output=1:NumberOutputs
        [U,D,V]=svd(squeeze(CurrentEx(Output,Input,:,:)));
        PotNeg=min(min(D));
        
        if PotNeg< LargestNeg
            LargestNeg=PotNeg;
        end
        if PotNeg < LargestLocalNeg
            LargestLocalNeg = PotNeg;
        end
        AverageNeg=AverageNeg+PotNeg;
        
        HermDiff=squeeze(CurrentEx(Output,Input,:,:))-squeeze(CurrentEx(Output,Input,:,:))';
        MaxDiff=max(max(abs(HermDiff)));
        
        if MaxDiff>MaxNonH
            MaxNonH = MaxDiff;
        end
        if MaxDiff>LargestLocalMaxNonH
            LargestLocalMaxNonH = MaxDiff;
        end
        HermDiffVec=reshape(HermDiff,[Dim^2,1]);
     SumSqV=HermDiffVec'*HermDiffVec;
     SumSq=SumSq+SumSqV;
     WeightedSumSq=WeightedSumSq+Decomposition{Ex,1}*SumSqV;
    end
    end
    AverageLNeg=AverageLNeg+LargestLocalNeg;
    AverageMax=AverageMax+LargestLocalMaxNonH;
end
disp('Largest Negative Eigenvalue Overall');
disp(LargestNeg);
disp('Average Largest Negative Eigenvalue');
disp(AverageLNeg/NumberOfExtremals);
disp('Average Negative Eigenvalue');
disp(AverageNeg/(NumberOfExtremals*NumberInputs*NumberOutputs));

disp('Hermiticity')
disp('Largest Elementwise Non-Hermicity')
disp(MaxNonH)
disp('Average Largest Element Non-Hermicity')
disp(AverageMax/NumberOfExtremals);
disp('Non-Weighted total non-Hermiticty (Sum of Squares)')
disp(SumSq)
disp('Weighted total non-Hermiticty (Sum of Squares)')
disp(WeightedSumSq)



disp('No-Signalling: values displayed in sum of squares')
LargestGap=0;
AverageLGap=0;
AverageGap=0;
if NumberInputs < 2 
disp('No no-signalling constraints!')
else
    for Ex=1:NumberOfExtremals
        CurrentEx=Decomposition{Ex,2};
        LargestLocalGap=0;
            for Input=2:NumberInputs
                for Input2=1:Input-1
                    MarginalDifference=sum(CurrentEx(:,Input,:,:),1)-sum(CurrentEx(:,Input2,:,:),1);
                    MarginalVec=reshape(MarginalDifference,[Dim^2,1]);
                    NSDiff=norm(MarginalVec)^2;
                    if NSDiff>LargestGap
                        LargestGap=NSDiff;
                    end
                    if NSDiff>LargestLocalGap
                        LargestLocalGap=NSDiff;
                    end
                    AverageGap=AverageGap+NSDiff;
                    
                end
            end
            AverageLGap=AverageLGap+LargestLocalGap;
    end
disp('Largest No-Sig violation');
disp(LargestGap);
disp('Average Largest no-Sig violation');
disp(AverageLGap/NumberOfExtremals);
disp('Average no-Sig violation');
disp(AverageGap/(NumberOfExtremals*NumberInputs*NumberOutputs));
end


        

