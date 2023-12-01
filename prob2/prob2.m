clearvars
close all

dudn = -3.8; nodQ = 7; ElemMeanTemp = 3; nodHint = 7;

temp = 30.0;
p = [0.25, 0.5];

 
plt = 1;     %plt = 1 --> draw figures, 
             %plt = 0 --> no figures (only computations)

kc = 1.0;
a11 = 1.0;
a12 = 0.0;
a21 = a12;
a22 = a11;
a00 = 0.0;
f = 0;

nodes = [ 
    0.0, 0.0;
    0.0, 1.0;
    1.0, 0.0;
    1.0, 1.0;
    2.0, 1.0;
    0.5, 0.5;
    1.5, 0.5;
    ];

elem = [
    1, 6, 2;
    1, 3, 6;
    3, 4, 6;
    4, 2, 6;
    3, 7, 4;
    5, 4, 7
    ];

numNodes = size(nodes,1);
numElem = size(elem,1);

if plt
    numbering = 1;
    %plotElements(nodes, elem, numbering)
    plotElementsOld(nodes, elem, numbering);
end

K = zeros(numNodes);
F = zeros(numNodes,1);
Q = zeros(numNodes,1);

coeff = [a11, a12, a21, a22, a00, f];

for e = 1:numElem
    rows = elem(e,:);
    cols = rows;
    [Ke,Fe]=linearTriangElement(coeff,nodes,elem,e);
    K(rows, cols) = K(rows,cols) + Ke;
    if (coeff(6) ~= 0)
        F(rows) = F(rows) + Fe;
    end
end

%BC 
fixedNodes = [1,2];
freeNodes = [3,4,5,6,7];

%Natural BC
q0 = kc*dudn;
indBC = [3,7,5];
[Q]=applyConstantNaturalBC(nodes,elem,indBC,q0,Q);
indBC = [5,4];
[Q]=applyConstantNaturalBC(nodes,elem,indBC,q0,Q);

indBC = [3,4,5,7]; %All nodes with natural BC
fprintf('Computation of Q''s using function ''applyConstantNaturalBC'':\n')
printNaturalBC(indBC,Q);

%We can compute Q's by hand
Q(3) = q0 * norm(nodes(3,:)-nodes(7,:))/2;
Q(4) = q0 * norm(nodes(4,:)-nodes(5,:))/2;
Q(5) = q0 * norm(nodes(4,:)-nodes(5,:))/2 + q0 * norm(nodes(5,:)-nodes(7,:))/2;
Q(7) = q0 * norm(nodes(3,:)-nodes(7,:))/2 + q0 * norm(nodes(7,:)-nodes(5,:))/2;

fprintf('Computation of Q''s ''by hand'':\n')
printNaturalBC(indBC,Q);

%Essential BC
u = zeros(numNodes,1);
u(fixedNodes) = temp;

%Reduced System;
Qm = Q(freeNodes) + F(freeNodes) - K(freeNodes,fixedNodes)*u(fixedNodes);
Km = K(freeNodes,freeNodes);

um = Km\Qm;
u(freeNodes) = um;

%Post Process
%
%Compute all the Q's
Q = K*u - F;

solutionTable = [(1:numNodes)',u,Q];
asts = repelem('*',28);
fprintf('\nPost-process: print table and plot temperatures (if plt = 1):\n')
% fprintf(['%28s\n',...
%          'Post Process:\n',...
%     'Print table results, and\n',...
%     'plot the temperatures\n',...
%     '%28s\n'],asts,asts)
fprintf('%4s%8s%11s\n','nod.','T','Q')
fprintf('%4d%12.4e%12.4e\n',solutionTable')

%Plot solution
if plt 
    title = 'Distribution of temperatures';
    colorPalette = 'jet';
    plotContourSolution(nodes,elem,u,title,colorPalette)
end

%Print out soluitons
%
meanTemp = zeros(numElem,1);

for e = 1:numElem
    nods = elem(e,:);
    vertexs = nodes(nods,:);
    [alphas, isInside] = baryCoord(vertexs, p);
    meanTemp(e) = sum(u(nods))/length(nods);
    if isInside >= 1
        interpTemp = alphas * u(nods);
        %break;
    end
end

%Solutions
fprintf('\nSolutions\n')
fprintf('(a) Q(%d) = %.4e\n', nodQ, Q(nodQ))
fprintf('(b) interpTemp  = %.4e\n', interpTemp)
fprintf('    Hint. Temperature at node %d: T = %.4e\n',nodHint,u(nodHint))
fprintf('(c) <T> at element %d: %.4e\n',ElemMeanTemp,meanTemp(ElemMeanTemp))


%%
% To print the natural BC
%
%%
function printNaturalBC(nodsBC, Q)
    for nod = nodsBC
        if nod == nodsBC(end)
            fprintf('Q(%d) = %.4e\n',nod, Q(nod))
        else
            fprintf('Q(%d) = %.4e, ',nod, Q(nod))
        end
    end
end