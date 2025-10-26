function BFTO_Buckling(lelx, lely, vol_con)
%% PARAMETERS DEFINITION
scale = 1; h = 1;
lelx = 240*scale; lely =120*scale;
BDY = [0, 0; lelx, lely] ;
vol_con = 0.35;
xmin = 1.e-9; dbeta = 0.5;
E0 = 1; nu = 0.3;
D = E0/(1-nu^2)*[1  nu  0;  nu  1  0;  0  0  (1-nu)/2] ;
rmin =  2; % FILTER RADIUS
Meshdiff = 5*h; d1 = 0.5*h; d2 = 1.0*h;
%% INITIAL MESHING
[xn, yn] = meshgrid(0: h :lelx, 0: h: lely) ;
[p, t, ~, ~, Ve, pmid, K, S, B1, BB, dGdu] = GenerateMesh(xn, yn, h, BDY, D, []) ;
% INITIAL DESIGN VARIABLES
xupp = ones(length(t), 1); xlow = ones(length(t), 1)*xmin;
xphys = max(xlow,min(xupp,ones(length(t),1)*vol_con));
x = xphys; xfilt = xphys;
% PREPARE FILTER
[ H, Hs] = FilterIndex(pmid, rmin, Ve);  % PREPARE FILTER
%% MAIN LOOP
loop = 0 ; loopbeta = 0;  tol = 0;
beta = 3;
while (tol == 0)  && (loop < 1000)
    xold = x; xpold = xphys; loop = loop + 1 ; loopbeta = loopbeta + 1;   beta1 = 1.2 * beta;
    % FEA AND SENSITIVITY ANALYSIS
    [dCe, C, str1, phi, lambda, mu, dmu, mu1, dmu1, muKS, dmuKS, muKS2, dmuKS2, BSE] = FEA(t, p, Ve, BDY, x, K, S, B1, BB, dGdu) ;
    BSE1 = reshape(log10(BSE(:, 1)/max(BSE(:, 1))), length(t), 1);
    figure(2); clf;
    colormap(brewermap([],'-Spectral')); clim([min(BSE1), max(BSE1)]);
    patch('Faces', t, 'Vertices', p, 'EdgeColor', 'none', 'FaceVertexCData', BSE1, 'FaceColor', 'flat') ; % PLOT OPTIMIZED STRUCTURE
    colorbar;  axis off equal tight ; drawnow ;
    dVe = Ve./sum(Ve); % SENSITIVITY OF VOLUME
    vol = sum( Ve.* xphys)/( lelx * lely) ;
    % SENSITIVITY TRANSFORMATION
    dc = dmuKS.*(1-tanh(beta1*(xfilt-0.5)).^2)*beta1/2/tanh(beta1*0.5);
    dv = dVe.*(1-tanh(beta1*(xfilt-0.5)).^2)*beta1/2/tanh(beta1*0.5);
    % SENSITIVITY AVERAGE
    if loopbeta > 1; dc = (dc+olddc)/2.; end
    olddc = dc;
    % GRAYNESS DEGREE
    WetE = Ve .* min(abs(xphys - 0), abs(xphys - 1));
    delta = sum(WetE) / sum(Ve);
    % OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES
    l1 = 0 ; l2 = 100 ;  move = 0.05 ;
    while l2-l1 > 1.0e-12
        lmid = 0.5*(l2+l1) ;
        x = max(xlow, min( xupp, xold.*(1+(-dc./(lmid.*dv)-1)/5)));
        xfilt = max(xlow,min(xupp,H*x(:)./Hs));
        xphys = (tanh(beta1*0.5)+tanh(beta1*(xfilt-0.5)))/2/tanh(beta1*0.5);
        xphys = max(xlow, max(xpold-move, min( xupp, min(xpold+move, xphys)))) ;
        if sum(xphys(:).*Ve(:)) > vol_con*lelx*lely, l1 = lmid; else l2 = lmid ; end
    end
    % IMPLICIT FLOATING PROJECTION (FP) CONSTRAINT FOR DESIGN VARIABLES
    l1 = 0; l2 = 1; sum_xnew = sum(x(:).*Ve(:));
    while (l2-l1) > 1.0e-12
        th = (l1+l2)/2.0;
        x = max(xlow,min(xupp,max(xold-move,min(xold+move,(tanh(beta*th)...
            +tanh(beta*(xfilt-th)))/(tanh(beta*th)+tanh(beta*(1.-th)))))));
        if sum(x(:).*Ve(:)) > sum_xnew; l1 = th; else; l2 = th; end
    end
    % RESULT DISPLAY
    figure(1); clf; colormap summer ;
    patch( 'Faces', t, 'Vertices', p, 'FaceVertexCData', xphys, 'FaceColor', 'flat') ;
    colorbar ; hold on ; axis off equal tight ; drawnow ;
    change = max(abs(xold(:)-x(:))) ;
    fprintf( 'It.:%5i Obj.:%8.4f lambda.:%8.4f mu1.:%8.4f muKS.:%8.4f Vol.:%7.3f ch.:%7.3f 0/1:%7.3f\n', loop, C, lambda(1), mu1, muKS, vol, change, delta) ;
    plotL1(loop) = lambda(1); plotL2(loop) = lambda(2); strL='\lambda_1(-)';
    plotL3(loop) = lambda(3); plotL4(loop) = lambda(4);
    figure(3); plot(1:loop, plotL1, 1:loop, plotL2, 1:loop, plotL3, 1:loop, plotL4); ylabel(strL); pause(1e-6);
    % INCREASE PROJECTION PARAMETER
    if ((loopbeta >= 5 && change < 0.01) || loopbeta >= 100)
        if delta < 0.01 && change < 0.02
            tol = 1;
        end
        if delta >= 0.01
            beta = beta + dbeta; loopbeta = 0;
            fprintf('Parameter beta increased to %g.\n',beta) ;
        end
    end
end
%% FINAL RESULT PLOTTTING
xBF = griddata(pmid(:,1), pmid(:,2), xphys , xn, yn, 'linear') ;
xBD = griddata(pmid(:,1), pmid(:,2), xphys , xn, yn, 'nearest') ;
xBF(isnan(xBF)) = xBD(isnan(xBF));
figure(5); clf;
[c] = ContourPoints(contour(xn, yn, xBF, [0.5 0.5]), d1, d2) ;
[p,~,t1,t2] = GenerateMesh(xn, yn, h, BDY, D, c, xBF, Meshdiff) ;
figure(5); clf;
patch('Faces', t1, 'Vertices', p, 'EdgeColor', 'k', 'FaceColor',[0 127 102]/255) ;  axis off equal tight
patch('Faces', t2, 'Vertices', p, 'EdgeColor', 'k', 'FaceColor',[255 255 102]/255) ; set(gca,'YDir','reverse') ;
end
%% FINITE ELEMENT ANALYSIS
function [dCe, C, str1, phi, lambda, mu, dmu, mu1, dmu1, muKS, dmuKS, muKS2, dmuKS2, BSE] = FEA(t, p, Ve, BDY, x, K, S, B1, BB, dGdu)
NT = length(t); NDOF = 2*length(p); nEig = 12; pnorm = 160;
elemDof = zeros(NT,6); elemDof(:,[1 3 5]) = 2*t-1; elemDof(:,[2 4 6]) = 2*t;
iK = reshape( kron( elemDof, ones(6,1))', 36*NT,1) ;
jK = reshape( kron( elemDof, ones(1,6))', 36*NT,1) ;
xi = repelem( x, 36, 1) ;
sK = xi.*reshape( K, 36*NT, 1) ;
NK = fsparse( iK, jK, sK, [NDOF, NDOF]) ;
fixedNodes = find(p(:,1)==BDY(1,1)) ;
forceNodes = find(p(:,1)==BDY(2,1) & (p(:,2)>=(BDY(1,2)+BDY(2,2))/2-4) & (p(:,2)<=(BDY(1,2)+BDY(2,2))/2+4));
fixedDof = [2*fixedNodes-1; 2*fixedNodes] ;
AllDof = 1:NDOF ;
freeDofs = setdiff( AllDof, fixedDof) ;
FM = 0.1/(length(forceNodes)-1);
F = sparse( 2*forceNodes-1, 1, -FM, NDOF, 1) ;
[F(2*forceNodes(2)-1),F(2*forceNodes(end)-1)] = deal(F(2*forceNodes(2)-1)/2,F(2*forceNodes(end)-1)/2);
U = zeros(NDOF, 1) ;
KT = decomposition(NK(freeDofs,freeDofs),'chol','lower');
U(freeDofs) = KT \ F(freeDofs) ;
C = F' * U / 2;
%% Buckling Analysis
phi = zeros(NDOF,nEig);  BSE = zeros(NT,nEig); adjL = phi; adjV = phi;
intp = x(:).^(1/3);
[dCe,str1,vms,Ge] = deal(zeros(NT,1) , zeros(NT,3) , zeros(NT,1), zeros(6,6*NT));
for i = 1 : NT
    dCe(i) = - U(elemDof(i,:))'*K(:,6*i-5:6*i)*U(elemDof(i,:)) ;
    str1(i,:) = intp(i).* S(:,6*i-5:6*i) * U(elemDof(i,:));
    vms(i) = sqrt(str1(i,1)^2+str1(i,2)^2-str1(i,1)*str1(i,2)+3*str1(i,3)^2);
    str2 = [str1(i,1) str1(i,3) 0 0;
        str1(i,3) str1(i,2) 0 0;
        0 0 str1(i,1) str1(i,3);
        0 0 str1(i,3) str1(i,2)];
    Ge(:,6*i-5:6*i) =  x(i) .* Ve(i) .* B1(:,6*i-5:6*i)' * str2 * B1(:,6*i-5:6*i);
end
sG = reshape( Ge, 36*NT, 1);
NG = fsparse( iK, jK, sG, [NDOF, NDOF]);
matFun = @(x) KT\(NG(freeDofs,freeDofs)*x);
[eivecs,D] = eigs(matFun,length(freeDofs),nEig,'sa','FailureTreatment','drop');
[mu,ii] = sort(diag(-D),'descend');
eivSort = eivecs(:,ii(1:nEig));
phi(freeDofs,:) = eivSort./sqrt(diag(eivSort'*NK(freeDofs,freeDofs)*eivSort)');
%% SENNSITIVITY ANALYSIS
[phiDKphi,phiDGphi,adj] = deal(zeros(NT,nEig));
for j = 1 : nEig
    adji = zeros(NDOF,1);
    phij = phi(:,j);  % Buckling vector j
    for i = 1 : NT
        BSE(i, j) = x(i).*phij(elemDof(i,:))'*K(:,6*i-5:6*i)*phij(elemDof(i,:));  % BUCKLING STRAIN ENERGY
        phiDKphi(i,j) = phij(elemDof(i,:))'*K(:,6*i-5:6*i)*phij(elemDof(i,:));
        phiDGphi(i,j) = phij(elemDof(i,:))'*Ge(:,6*i-5:6*i)*phij(elemDof(i,:));
        for k = 1 : 6
          phiDGDuphi = phij(elemDof(i,:))'*(intp(i).*x(i).*dGdu(:,36*i+6*k-41:36*i+6*k-36))*phij(elemDof(i,:));
          adji(elemDof(i,k)) = adji(elemDof(i,k))+ phiDGDuphi;
        end
    end
    adjL(:,j) = adji;
end
adjV(freeDofs,:) = KT \ adjL(freeDofs,:);
for j = 1 : nEig
    V = adjV(:,j);
    for i = 1 : NT
        adj(i,j) = U(elemDof(i,:))'*K(:,6*i-5:6*i)*V(elemDof(i,:));
    end
end
lambda = 1./mu;
dmu = -(phiDGphi+mu(1:nEig )'.*phiDKphi-adj);
mu1 = mu(1);
dmu1 = reshape(dmu(:,1), NT, 1);
fKS = @(p,v)max(v)+log(sum(exp(p*(v-max(v)))))/p;  % KS AGGREGATION
dKS = @(p,v,dv)sum(exp(p.*(v-max(v)))'.*dv,2)./sum(exp(p.*(v-max(v))));
muKS = fKS(pnorm,mu(1:nEig));
dmuKS = reshape(dKS(pnorm,mu(1:nEig),dmu(:, 1:nEig)), NT, 1);
muKS2 = fKS(pnorm,mu(2:nEig));
dmuKS2 = reshape(dKS(pnorm,mu(2:nEig),dmu(:, 2:nEig)), NT, 1);
end
%% FIND POINTS ON THE CONTOUR
function [c] = ContourPoints(c,d1,d2)
num = 1; col = 1;
while col < size(c,2)
    idx = col+1:col+c(2,col);
    s(num).ck = c(:,idx);
    s(num).sed = norm(s(num).ck(:,1)-s(num).ck(:,end));
    s(num).isopen = s(num).sed>1e-12;
    num = num+1; col = col+c(2,col)+1;
end
c = [];
for k = 1:num-1
    ck = s(k).ck;
    if length(ck)>4
        ndist = vecnorm(diff(ck,1,2),2,1);
        ck(:,find(ndist < d1) +1) = 0.5*(ck(:,ndist < d1)+ck(:,find(ndist < d1) +1));  % REMOVE TOO-CLOSE NODES
        ck(:,ndist < d1) = [];
        if s(k).sed < d1
            ck(:, end) = [];
        end
        if ~isempty(ck)
            if ~s(k).isopen && ~isempty(ck)
                ck = [ck ck(:,1)];
            end
            ndist = vecnorm(diff(ck,1,2),2,1);
            ct = ck; ck = ct(:,1);
            for i = 1:length(ndist)
                if  ndist(i) > d2
                    ck = [ck 0.5*(ct(:,i)+ct(:,(i+1))) ct(:,i+1)]; % INCLUDE ADDITIONAL NODES
                else
                    ck = [ck ct(:,i+1)];
                end
            end
        end
        c = [c; ck(:,1:end-(~s(k).isopen))'];
    end
end
end
%% BODY FITTED MESH GENERATOR
function [p,t,t1,t2,Ve,pmid,K,S,B1,BB,dGdu] = GenerateMesh(xn, yn, h, BDY, D, c, dN, Meshdiff)
if isempty(c) == 1
    [x,y] = meshgrid(BDY(1,1):h: BDY(2,1), BDY(1,2):h: BDY(2,2)) ;
    x(2:2:end,:)=x(2:2:end,:)+0.5*h; pi=[x(:),y(:)];
    pi(:,1) = min(BDY(2,1),max(BDY(1,1),pi(:,1))); pi(:,2) = min(BDY(2,2),max(BDY(1,2),pi(:,2)));
    p = [ BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1),(BDY(1,2)+BDY(2,2))/2; pi];
    p = unique( p, 'rows', 'stable');
    t = delaunay(p);
    t1 = [] ; t2 = [] ;
else
    [x,y] = meshgrid(BDY(1,1): h: BDY(2,1), BDY(1,2): sqrt(3)/2*h : BDY(2,2)) ;
    x(2:2:end,:)=x(2:2:end,:)+0.5*h; pi=[x(:),y(:)];
    dist = min(pdist2(pi, c), [], 2);
    d = rescale(dist, h, Meshdiff);
    r0 = 1./d.^2;
    pfix=[c; BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1),(BDY(1,2)+BDY(2,2))/2];
    pfix = unique(pfix,'rows','stable');
    load('randommatrixmc.mat','ranmatrix');   % LOAD SAVED RANDOM MATRIX
    p=[pfix; pi(ranmatrix(1:size(pi,1))<r0./max(r0),:)];  % DELETE NODES USING REJECTION METHOD
    p(:,1) = min(BDY(2,1),max(BDY(1,1),p(:,1))); p(:,2) = min(BDY(2,2),max(BDY(1,2),p(:,2)));
    L1xy = reshape(sqrt(min((xn(:) - c(:,1)').^2 + (yn(:) - c(:,2)').^2, [], 2)),size(xn));
    L1xy = rescale(L1xy, h, Meshdiff); % DESIRED EDGE LENGTH L1
    p1 = 0;
    % NOVE-MOVING LOOPS
    for i = 1:200
        if max(sum((p-p1).^2,2))>0.1*h
            p = unique(p,'rows','stable');
            t = delaunay(p);  % DELAUNAY TRIANGULATION
            edges = unique(sort([t(:,[1,2]);t(:,[1,3]);t(:,[2,3])],2),'rows') ;
            p1 = p;
        end
        L = sqrt(sum((p(edges(:,1),:)-p(edges(:,2),:)).^2,2));
        pm = (p(edges(:,1),:)+p(edges(:,2),:))/2;
        L1 = interp2(xn,yn,L1xy,pm(:,1),pm(:,2),'cubic');  % INTERPOLATE DENSITY INTO EDGE CENTROIDS
        L0 = 1.2*L1*sqrt(sum(L.^2)/sum(L1.^2));
        Fb = max(L0-L,0)./L *[1,1].*(p(edges(:,1),:)-p(edges(:,2),:));
        Fp = full(sparse(edges(:,[1,1,2,2]),ones(size(L1))*[1,2,1,2],[Fb,-Fb],size(p,1),2));
        Fp(1:size(pfix,1),:) = 0 ;
        p = p + 0.2*Fp ;  % MOVE NODE ACCORDING TO VIRTUAL FORCE
        p(:,1) = min(BDY(2,1),max(BDY(1,1),p(:,1))) ;
        p(:,2) = min(BDY(2,2),max(BDY(1,2),p(:,2))) ;
        if max(sum(0.2*Fp,2))<0.05*h, break; end
    end
    Ve = zeros(length(t),1) ;
    for i = 1:length(t)
        J = [p(t(i,1),1)-p(t(i,3),1) p(t(i,1),2)-p(t(i,3),2) ; p(t(i,2),1)-p(t(i,3),1) p(t(i,2),2)-p(t(i,3),2)];
        Ve(i) = 0.5.*det(J);  % ELEMENT VOLUME
    end
    t(Ve == 0,:) = [];
    pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
    dE = interp2(xn,yn,dN,pmid(:,1),pmid(:,2),'cubic');  % INTERPOLATE DENSITY INTO ELEMENT CENTROIDS
    t1=[t(dE<0.5,:)]; t2=t(dE>=0.5,:); t=[t1;t2];
end
pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
Ve = zeros(length(t),1) ;
ps = p;
for i = 1:length(t)
    J = [ps(t(i,1),1)-ps(t(i,3),1) ps(t(i,1),2)-ps(t(i,3),2) ; ps(t(i,2),1)-ps(t(i,3),1) ps(t(i,2),2)-ps(t(i,3),2)];
    Ve(i) = 0.5.*det(J);  % ELEMENT VOLUME
    Be = 1/det(J)*[J(2,2) 0 -J(1,2) 0 -J(2,2)+J(1,2) 0; 0 -J(2,1) 0 J(1,1) 0 J(2,1)-J(1,1); -J(2,1) J(2,2) J(1,1) -J(1,2) J(2,1)-J(1,1) -J(2,2)+J(1,2)];
    B1(:,6*i-5:6*i) = 1/det(J)*[J(2,2) 0 -J(1,2) 0 -J(2,2)+J(1,2) 0; -J(2,1) 0 J(1,1) 0 J(2,1)-J(1,1) 0; 0 J(2,2) 0 -J(1,2) 0 -J(2,2)+J(1,2); 0 -J(2,1) 0 J(1,1) 0 J(2,1)-J(1,1)];
    K(:,6*i-5:6*i) = det(J)/2*Be'*D*Be;
    S(:,6*i-5:6*i) = D * Be;  % STRESS MATRIX
    BB(:,6*i-5:6*i) = [J(2,2)^2  J(2,2)*(-J(1,2))  J(2,2)*(-J(2,2)+J(1,2))  (-J(2,1))^2  (-J(2,2)+J(1,2))*(-J(2,1))  (-J(2,2)+J(1,2))^2;
        (-J(2,1))^2  (-J(2,1))*J(1,1)  (-J(2,1))*(J(2,1)-J(1,1))  J(1,1)^2  (J(2,1)-J(1,1))*J(1,1)  ((J(2,1)-J(1,1)))^2;
        2*(-J(2,1))*J(2,2)  -J(2,1)*(-J(1,2))+J(2,2)*J(1,1)  -J(2,1)*(-J(2,2)+J(1,2))+J(2,2)*(J(2,1)-J(1,1)) ...
        2*J(1,1)*(-J(1,2))  -J(1,1)*(-J(2,2)+J(1,2))+(-J(1,2))*(J(2,1)-J(1,1))  2*(J(2,1)-J(1,1))*(-J(2,2)+J(1,2))];
end
%% Buckling analysis for dG/du
dGdu = zeros(6,36*length(t));
for i = 1:length(t)
    for ii = 1 : 6
        Uvec = zeros(6,1); Uvec(ii,1) = 1;
        str1 = S(:,6*i-5:6*i)*Uvec;
        str2 = [str1(1) str1(3) 0 0;
            str1(3) str1(2) 0 0;
            0 0 str1(1) str1(3);
            0 0 str1(3) str1(2)];
        dGduii = Ve(i) .*B1(:,6*i-5:6*i)' * str2 * B1(:,6*i-5:6*i);
        dGdu(:,36*i+6*ii-41:36*i+6*ii-36) = dGduii;
    end
end
end
%% PREPARE FILTER
function [H, Hs] = FilterIndex(pmid,rmin,Ve)
maxN =  ceil (100 * rmin);    NT = length (pmid);
iH = ones(NT*maxN, 1) ;
jH = ones(NT*maxN, 1) ;
sH = zeros(NT*maxN, 1) ;
Start = 1;
for i = 1 : NT
    pni = find((pmid(:,1) - pmid(i,1)).^2 + (pmid(:,2) - pmid(i,2)).^2 <= rmin^2);
    Weight = max((rmin - pdist2(pmid(i,:),pmid(pni,:))),0) .* (Ve(pni))';
    pEnd = Start + length(Weight);
    iH(Start:pEnd-1) = i;
    jH(Start:pEnd-1) = pni;
    sH(Start:pEnd-1) = Weight;
    Start = pEnd;
end
iH = iH(1:Start - 1);
jH = jH(1:Start - 1);
sH = sH(1:Start - 1);
H = sparse( iH, jH, sH) ;
Hs = sum( H, 2) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by ZC. Zhuang and Y. M. Xie                %
% Centre for Innovative Structures and Materials, RMIT University         %
% Please send your comments to: zhuanginhongkong@outlook.com              %
%                                                                         %
% The program is proposed for educational purposes and is introduced in   %
% the paper - An Integrated Framework for Body-Fitted Topology            %
% Optimization Considering Buckling Load Factor, Stiffness,               %
% and von Mises Stress, CMAME, 2025                                       %
%                                                                         %
% Disclaimer:                                                             %
% The authors reserve all rights but do not guarantee that the code is    %
% free from errors. Furthermore, we shall not be liable in any event.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%