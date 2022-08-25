function [idassem] = identity_assemblage_creator(InputNo,OutputNo,Dim);

%This function outputs an (#outs,#ins,Dim,Dim) array which corresponds to
%an assemblage where every substate = identity(Dim)/Dim.
idassem=zeros(OutputNo,InputNo,Dim,Dim);
for Input=1:InputNo
    for Output=1:OutputNo
        idassem(Output,Input,:,:)=eye(Dim)/(OutputNo*Dim);
    end
end
end
        