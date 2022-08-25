function [CleanAssemblage] = assemblage_clean(Assemblage,ZeroTol)

%%%%%%%%%%%%%%%%%%%%
%This function "cleans up" the assemblage by making it again Hermitian, and
%getting rid of small negatives.

%Inputs: 

%Assemblage: an #Out * #In * Dim * Dim array

%ZeroTol: If something is unbalanced / too negative, we will spit out a
%warning.

Dimensions=size(Assemblage);
NumberOfOutputs=Dimensions(1);
NumberOfInputs=Dimensions(2);
Dimension=Dimensions(3);

%We will output a definitely Hermitian, definitely positive set of
%substates.
CleanAssemblage=Assemblage;

for Input=1:NumberOfInputs
    for Output=1:NumberOfOutputs
        
        %First: we make each substate Hermitian
        CurrentState=squeeze(Assemblage(Output,Input,:,:));
        SymState=1/2*(CurrentState+CurrentState');
        
        %If the symmetry processs changes the state too much, we spit out
        %an error.
        if max(max(abs(SymState-CurrentState))) > ZeroTol
            warning('Assemblage became asymmetric over tolerance');
        end
        
        
        %Next, we get rid of any small eigenvalues, or negative
        %eigenvalues.
        %It is ok to get rid of such small eigenvalues, since they would be
        %ignored by the perturbation generator as well.
        
        [U,D,V]=svd(SymState);
        MinD=min(diag(D));
        
        %Once again, we issue a warning if too large a negative value
        %appears.
        if MinD< -ZeroTol
            warning('Assemblage became negative over tolerance');
            
        end
        PosD=D>=ZeroTol;
        ImprovedState=U*D*PosD*V';
        
        %We even re-symmetrise, just in case!
        CleanAssemblage(Output,Input,:,:)=1/2*(ImprovedState+ImprovedState');
        
    end
end

        
        
        
        
       
