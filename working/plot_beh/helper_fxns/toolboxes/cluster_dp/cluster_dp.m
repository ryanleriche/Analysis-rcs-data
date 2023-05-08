function [dec_fig, cl,halo] = cluster_dp(cfg, datcluster)

% do clustering
% get distance matrix
D                = pdist(datcluster, @naneucdist);
distmat          = squareform(D);
distMatrices     = squareform(distmat,'tovector');

% get row indices
rows             = repmat(1:size(distmat,1),size(distmat,2),1)';
idx              = logical(eye(size(rows)));
rows(idx)        = 0;
rowsColmn        = squareform(rows,'tovector');

% get column idices
colmns           = repmat(1:size(distmat,1),size(distmat,2),1);
idx              = logical(eye(size(colmns)));
colmns(idx)      = 0;
colsColmn        = squareform(colmns,'tovector');
% save data for rodriges
distanceMat      = [];
distanceMat(:,1) = rowsColmn;
distanceMat(:,2) = colsColmn;
distanceMat(:,3) = distMatrices;




%%
xx = distanceMat;
% xx=load('example_distances.dat');
ND = max(xx(:,2));
NL = max(xx(:,1));
if (NL>ND)
  ND=NL;
end
N=size(xx,1);

for i=1:ND
  for j=1:ND
    dist(i,j) = 0;
  end
end

for i=1:N
  ii = xx(i,1);
  jj = xx(i,2);
  dist(ii,jj) = xx(i,3);
  dist(jj,ii) = xx(i,3);
end

percent= 1.25;
fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);

position = round(N*percent/100);
sda=sort(xx(:,3));
dc = sda(position);

fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);


for i=1:ND
  rho(i)=0.;
end
%
% Gaussian kernel
%
for i=1 : ND-1
  for j=i+1:ND
     rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
     rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
  end
end
%
% "Cut off" kernel

% for i=1:ND-1
%  for j=i+1:ND
%    if (dist(i,j)<dc)
%       rho(i)=rho(i)+1.;
%       rho(j)=rho(j)+1.;
%    end
%  end
% end

maxd=max(max(dist));

[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;

for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);
     end
   end
end
delta(ordrho(1))=max(delta(:));



for i=1:ND
  ind(i)   = i;
  gamma(i) = rho(i)*delta(i);
end

dec_fig = figure('Units', 'Inches', 'Position', [0, 0, 7, 7]);

subplot(211);

plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');

title ('Decision Graph','FontSize',15.0)
xlabel ('Local Density, \rho'); 
ylabel ('min(Dist to higher \rho point), \delta'); ylim([0, max(delta(:))*1.2]);
set(gca, 'FontSize', 10)

%% Draw rectangle of interest on density versus distance plot
if strcmp(cfg.CBDP_method, 'manual')

    disp('Select a rectangle enclosing cluster centers')
    
    rect =  drawrectangle;
    
    rho_min    = rect.Position(1);
    delta_min  = rect.Position(2);

elseif strcmp(cfg.CBDP_method, 'top_two')

    [~, i_sort] = sort(gamma, "descend");

    rho_min   = min(rho(i_sort(1:2)));

    %[~, i_sort] = sort(delta, "descend");
    delta_min = min(delta(i_sort(1:2)));

    

end
%%
t = '';

NCLUST=0;

for i=1:ND
  cl(i)=-1;
end

for i=1:ND
  if ((rho(i)>=rho_min) && (delta(i)>=delta_min))
     NCLUST=NCLUST+1;
     cl(i)=NCLUST;
     icl(NCLUST)=i;
  end
end

fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);
disp('Performing assignation')

%assignation
for i=1:ND
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(nneigh(ordrho(i)));
  end
end
%halo
for i=1:ND
  halo(i)=cl(i);
end
if (NCLUST>1)
  for i=1:NCLUST
    bord_rho(i)=0.;
  end
  for i=1:ND-1
    for j=i+1:ND
      if ((cl(i)~=cl(j))&& (dist(i,j)<=dc))
        rho_aver=(rho(i)+rho(j))/2.;
        if (rho_aver>bord_rho(cl(i))) 
          bord_rho(cl(i))=rho_aver;
        end
        if (rho_aver>bord_rho(cl(j))) 
          bord_rho(cl(j))=rho_aver;
        end
      end
    end
  end
  for i=1:ND
    if (rho(i)<bord_rho(cl(i)))
      halo(i)=0;
    end
  end
end
for i=1:NCLUST
  nc=0;
  nh=0;
  for j=1:ND
    if (cl(j)==i) 
      nc=nc+1;
    end
    if (halo(j)==i) 
      nh=nh+1;
    end
  end
    t = [t, newline, 'CLUSTER: ', num2str(i), ' CENTER : ', num2str(icl(i)),...
      ' ELEMENTS: ', num2str(nc), ' CORE: ', num2str(nh), ' HALO: ', num2str(nc - nh)];

end

TextLocation(t);
cmap = colormap(brewermap([],"Dark2"));

for i = 1 : NCLUST
   ic = int8((i*64.)/(NCLUST*1.));
   hold on
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:));
end


subplot(212)

ranked_gamma = sort(gamma, "descend") / max(gamma);

for i = 1: length(ranked_gamma)

    if (ranked_gamma(i) >= rho_min * delta_min / max(gamma))
    
        plot(i, ranked_gamma(i), 'o','Markersize', 10, 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', cmap(i,:))
        hold on
    
    else

    plot(i, ranked_gamma(i), 'o', 'Markersize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); 

    end
end


title ('Cluster Center Density Distribution','FontSize',15.0)

xlabel ('Every point (possible cluster centers)')
ylabel ('\rho * \delta / max(\rho * \delta)')
set(gca, 'FontSize', 10, 'TickLength', [0, 0])

% perform nonclassical multidimensional scaling before visualization
%{
subplot(2,1,2)


disp('Performing 2D nonclassical multidimensional scaling')

Y1 = mdscale(dist, 2, 'criterion','metricsstress');

plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');

title ('2D Nonclassical multidimensional scaling','FontSize',15.0)
xlabel ('X')
ylabel ('Y')
for i=1:ND
 A(i,1)=0.;
 A(i,2)=0.;
end
for i=1:NCLUST
  nn=0;
  ic=int8((i*64.)/(NCLUST*1.));
  for j=1:ND
    if (halo(j)==i)
      nn=nn+1;
      A(nn,1)=Y1(j,1);
      A(nn,2)=Y1(j,2);
    end
  end
  hold on
  plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end

%for i=1:ND
%   if (halo(i)>0)
%      ic=int8((halo(i)*64.)/(NCLUST*1.));
%      hold on
%      plot(Y1(i,1),Y1(i,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
%   end
%end
faa = fopen('CLUSTER_ASSIGNATION', 'w');
disp('Generated file:CLUSTER_ASSIGNATION')
disp('column 1:element id')
disp('column 2:cluster assignation without halo control')
disp('column 3:cluster assignation with halo control')
for i=1:ND
   fprintf(faa, '%i %i %i\n',i,cl(i),halo(i));
end
%}

function D2 = naneucdist(XI,XJ)  
%NANEUCDIST Euclidean distance ignoring coordinates with NaNs
n = size(XI,2);
sqdx = (XI-XJ).^2;
nstar = sum(~isnan(sqdx),2); % Number of pairs that do not contain NaNs
nstar(nstar == 0) = NaN; % To return NaN if all pairs include NaNs
D2squared = sum(sqdx,2,'omitnan').*n./nstar; % Correction for missing coordinates
D2 = sqrt(D2squared);
end
end