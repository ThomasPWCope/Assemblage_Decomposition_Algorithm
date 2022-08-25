%This is example 1 in the main paper.

PentagonExtremal=zeros(5,1,2,2);
for i=1:5
BlochVector=[cos((i-1)*(2*pi/5)),sin((i-1)*(2*pi/5)),0];
PentagonExtremal(i,1,:,:)=(1/10)*(eye(2) + BlochVector(1)*Pauli(1)+BlochVector(2)*Pauli(2)+BlochVector(3)*Pauli(3));
end
PentDecomp=zero_recursion_H([],PentagonExtremal,1,1,10^-6);


PentDecomp=congregate_extremals(PentDecomp,10^-5);

% We want to summarise our data.
Weights=zeros(1,length(PentDecomp));
BlochVectors=zeros(3,length(PentDecomp));
for i = 1:length(PentDecomp)
Weights(i)=PentDecomp{i,1};
Traces=[trace(squeeze(PentDecomp{i,2}(1,1,:,:))), trace(squeeze(PentDecomp{i,2}(2,1,:,:))) ,...
    trace(squeeze(PentDecomp{i,2}(3,1,:,:))), trace(squeeze(PentDecomp{i,2}(4,1,:,:))), ...
    trace(squeeze(PentDecomp{i,2}(5,1,:,:)))];
[MaxT,PMax]=max(Traces);
[EVec,EVal]=eigs(squeeze(PentDecomp{i,2}(PMax,1,:,:)));
Rho=EVec(:,1)*ctranspose(EVec(:,1));
BlochVectors(1,i)=2*real(Rho(1,2));
BlochVectors(2,i)=2*imag(Rho(1,2));
BlochVectors(3,i)=2*real(Rho(1,1))-1;
end

%save('weights.mat','weights');
%save('blochvectors.mat','blochvectors');
     