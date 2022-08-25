function [OutputListOfExtremals]=congregate_extremals(InputListOfExtremals,EqualityTol)

%%%%%%%%%%%%%%%%%%%%
%This function takes in a list of extremal assemblages, and sorts them into
%a smaller list in which all identical assemblages are combined. 

%%%%%%%%%%%%%%
%Inputs:

%InputListOfExtremals: The original list of extremal assemblages
% EqualityTol: The difference (in 1-norm) between two assemblages, for them
%to be classed as the same assemblage.

%Outputs: 

%OutputListOfExtremals: a smaller list of extremal inequalities. 

%%%%%%%%%%%%%%%%%%%%%%%%%

NumberExtremals=length(InputListOfExtremals);
NumberNewExtremals=0;

OutputListOfExtremals=cell(0,2);

for i=1:NumberExtremals
    
CurrentRow=InputListOfExtremals(i,:);
CurrentProb=CurrentRow{1};
CurrentAssem=CurrentRow{2};
Match=0;
CurrentNewExtremal=1;

while Match==0 && CurrentNewExtremal<=NumberNewExtremals
CurrentNewAssem=OutputListOfExtremals{CurrentNewExtremal,2};
DistanceVal=sum(abs(CurrentNewAssem-CurrentAssem),'all');

if DistanceVal < EqualityTol
Match=1;
OutputListOfExtremals{CurrentNewExtremal,1}=OutputListOfExtremals{CurrentNewExtremal,1}+CurrentProb;
else
CurrentNewExtremal=CurrentNewExtremal+1;
end
end
if Match==0
OutputListOfExtremals{CurrentNewExtremal,1}=CurrentProb;
OutputListOfExtremals{CurrentNewExtremal,2}=CurrentAssem;
NumberNewExtremals=NumberNewExtremals+1;
end
end

