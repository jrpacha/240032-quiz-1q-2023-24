clearvars
close all

nodA = 1; %Node A
nodB = 2; %Node B
nodC = 3; %Node C

plt = 1;              %plt = 1 --> draw figures, 
                      %plt = 0 --> no figures (only computations)                         
secArea=2550.0;       %section Area (in mm^2);
youngM=175.0;         %Young modulus, in kN/mm^2 (1GP=1.0 kN/mm^2)

barLengthBL = 1.8e3;    %mm
normForce = 250.0; %kN

%Geometry
theta = linspace(0, 2*pi-pi/3, 6);
z = barLengthBL*exp(1i*theta);
nodes = [real(z), 0; imag(z), 0]';
elem = [
    7, 1;
    7, 2;
    7, 3;
    7, 4;
    7, 5;
    7, 6;
    1, 2;
    2, 3;
    3, 4;
    4, 5;
    6, 1
    ];
      
numNod=size(nodes,1);
numElem=size(elem,1);
dim=size(nodes,2);

if plt
    numbering=1;
    %plotElements(nodes, elem, numbering);
    plotElementsOld(nodes, elem, numbering);
end

%Real constants: Materials and sections area
A=secArea*ones(1,numElem);
E=youngM*ones(1,numElem);

%Assembly
u=zeros(dim*numNod,1);
Q=zeros(dim*numNod,1);
K=zeros(dim*numNod);

for e=1:numElem
    Ke=planeLinkStiffMatrix(nodes,elem,e,E,A);
    rows=[2*elem(e,1)-1;2*elem(e,1);2*elem(e,2)-1;2*elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke; %Assembly
end

%Loads
%node 1:
nod=1;
Q(dim*nod-1)=normForce;                 %kN
Q(dim*nod)=0;                           %kN

%node 2:
nod=2;
Q(dim*nod-1)=normForce*cos(theta(nod)); %kN
Q(dim*nod)=normForce*sin(theta(nod));   %kN

%node 3:
nod=3;
Q(dim*nod-1)=normForce*cos(theta(nod)); %kN
Q(dim*nod)=normForce*sin(theta(nod));   %kN

%Boundary Conditions
fixedNods=[];
%node 4:
nod=4;
fixedNods=[fixedNods,dim*nod-1];    %(u4_x=0);
u(dim*nod-1)=0.0;

%node 5:
nod=5;
fixedNods=[fixedNods,dim*nod-1];    %(u5_x=0); 
fixedNods=[fixedNods,dim*nod];      %(u5_y=0);  

%node 6;
nod=6;
fixedNods=[fixedNods,dim*nod-1];    %(u6_x=0); 
fixedNods=[fixedNods,dim*nod];      %(u6_y=0);  

%Reduced system
freeNods=setdiff(1:dim*numNod,fixedNods);
Qm=Q(freeNods,1)-K(freeNods,fixedNods)*u(fixedNods);
Km=K(freeNods,freeNods);

%%
%Solve the reduced system
um=Km\Qm;
u(freeNods)=um;

%Postp-rocess
%Reaction forces
R = K*u-Q;

%Print out position, displacements, applied forces and reaction forces
displ = [u(1:2:end),u(2:2:end)];
forces = [Q(1:2:end), Q(2:2:end)];
reactF = [R(1:2:end), R(2:2:end)];
fprintf('\n%6s%10s%14s%14s%14s%14s%14s\n','NOD.','UX(mm)','UY(mm)','QX(kN)','QY(kN)','RX(kN)','RY(kN)')
fprintf('%4d%14.5e%14.5e%14.5e%14.5e%14.5e%14.5e\n',[(1:numNod)',displ,forces,reactF]')

%Post-process
%Show the original structure and the deformed one
%figure()
if plt
    esc=100; %scale factor to magnify displacements
    plotDeformedTruss(nodes, elem, u, esc);
end

%
%Part (a)
%
%Final distance (after the forces are applied) between node A (node 1) and C (node 3)
%
finalDistAC = norm(nodes(3,:)-nodes(1,:)+displ(3,:)-displ(1,:));
fprintf('\n(a) Final distance between nodes A and C, dist = %.4e\n',finalDistAC)
fprintf('    Hint. The horizontal displacement of node B is UBX = %.4e\n',u(2*nodB-1))

%
%Part (b)
%
% Elements'' final lenghts'' sum
finalBarLength = zeros(numElem,1);
for e = 1:numElem
    nod1 = elem(e,1); nod2 = elem(e,2);
    finalBarLength(e) = norm(nodes(nod2,:)-nodes(nod1,:)+displ(nod2,:)-displ(nod1,:));
end
fprintf('(b) Elements'' final lengths'' sum, sum = %.4e\n', sum(finalBarLength))

%
%Part (c)
%Max stress of elements
%
sigma = zeros(numElem,1);
for e = 1:numElem
    nod1 = elem(e,1);
    nod2 = elem(e,2);
    vr = nodes(nod2,:)-nodes(nod1,:);
    normVRSquared = vr*vr';
    sigma(e) = E(e)*(displ(nod2,:)-displ(nod1,:))*vr'/normVRSquared;
end
fprintf('(c) Maximum of absolute value''s elements'' stress, max(|sigma|) = %.4e\n',norm(sigma,inf))