%% setup

[r,c]=size(ics);

% calculate appropriate bin number
nbin=ceil(range(ics)./(2*iqr(ics)*r.^(-1/3)));
nbin=repmat(median(nbin),1,c);

% assign values to bins
z=bsxfun(@minus,ics,min(ics));
z=ceil(bsxfun(@times,z,nbin./max(z)));

% ensure values correspond to appropriate integers
for n=1:c
    z(z(:,n)>nbin(n),n)=nbin(n);
end
z(~z)=1;

%% pairwise pdfs

o=triu(ones(c),1);
[r1 c1]=ind2sub(size(o),find(o==1));
K=cell(size(r1));
nr1=length(r1);
parfor k=1:nr1
  K{k}=sparse(z(:,r1(k)),z(:,c1(k)),1).*(1/r);  
end
[idx]=sub2ind(size(o),r1,c1);
% P=cell(c);
% for k=1:length(idx)
% P(r1(k),c1(k))=K(k);
% end

i_nan    = any(isnan([REDcap.RCS02.unpleasantVAS, REDcap.RCS02.painVAS]),2);
vas      = REDcap.RCS02.painVAS(~i_nan);
unpl_vas = REDcap.RCS02.unpleasantVAS(~i_nan);
%% find minimum of pairwise MI

% K == (X-by-1) sparse 2-dim jpdfs
MI = zeros(nr1,1);
for k=1:nr1
    p_mn=full(K{k}); %populate the jpdf with zeros
    [p_m,p_n]=meshgrid(sum(p_mn,2),sum(p_mn,1));
    tMI = p_mn.*log2(p_mn*(1./p_m)*(1./p_n));
    % NOTE: when p_m == 0 | p_n == 0, p_mn = 0
    % 0/0 :: NaN
    MI(k)=nansum(tMI(:));
end
[~,bestMI]=min(MI);

%%
% %
% % cp=cumprod(nbin);
% % % sort(z)
% % mfactor=[1 cp(1:c-1)];
% % nz=sum(bsxfun(@times,z,mfactor),2);
% % [count,uz]=hist(nz,unique(nz));
% %
% % v=zeros(numel(uz),c);
% % for n=c:-1:2
% %     v(:,n)=floor(uz./cp(n-1))+1;
% %     uz = mod(uz,cp(n-1));
% % end
% % v(:,1)=uz+1;

%% 2 tier sort

sc1=1;%1st priorty sorting component #
sc2=2;

%1st tier sort
[~,i]=sort(z(:,sc1));
y=z(i,:);
%find first occur of 1st tier
binind=1:r;
jumps=[1;diff(y(:,sc1))];
binoc=[binind(jumps>0),r+1];
%2nd tier sort
tmp=y(:,sc2);
for n=1:length(binoc)-1
    b1=binoc(n);
    b2=binoc(n+1);
    [~,i]=sort(tmp(b1:b2-1));
    ttmp=y(b1:b2-1,:);
    y(b1:b2-1,:)=ttmp(i,:);
end
%compile all bin change occur
jumps=[1;diff(y(:,sc2))];
binoc=sort(unique([binind(jumps>0),binoc]));

%% n dimensional pdf
fprintf('start ndim\n');
tic
w=cell(1,numel(binoc));
parfor h=1:numel(binoc)-1
    b1=binoc(h);
    b2=binoc(h+1);
    tz=y(b1:b2-1,:);
    if b2-b1>1
        wh=zeros(c+1,b2-b1+1);
        co=1;
        for k=1:b2-b1
            zv=tz(k,1:c);
            v=all(bsxfun(@eq,wh(1:c,1:co),zv'));
            %v=sum(bsxfun(@eq,wh(1:c,1:co),zv'),1)==c;
            if any(v)
                wh(c+1,v)=wh(c+1,v)+1;
            else
                co=co+1;
                wh(:,co)=[zv 1];
            end
        end
%         wh(:,[1 co+1:b2-b1+1])=[];
%         w{h}=wh;
        w{h}=wh(:,2:co);
    else
        w{h}=[tz 1]';
    end
    %fprintf('%d',h);
end
fprintf('\n');
w=cell2mat(w);
w(c+1,:)=w(c+1,:).*(1/r);
toc;

%% 2 tier sort

fprintf('start ndim-1\n');
sc1=1;%1st priorty sorting component #
sc2=2;
b=1;
y=z(:,setdiff(1:c,b));
c_1=c-1;
%1st tier sort
[~,i]=sort(y(:,sc1));
y=y(i,:);
%find first occur of 1st tier
binind=1:r;
jumps=[1;diff(y(:,sc1))];
binoc=[binind(jumps>0),r+1];
%2nd tier sort
tmp=y(:,sc2);
for n=1:length(binoc)-1
    b1=binoc(n);
    b2=binoc(n+1);
    [~,i]=sort(tmp(b1:b2-1));
    ttmp=y(b1:b2-1,:);
    y(b1:b2-1,:)=ttmp(i,:);
end
%compile all bin change occur
jumps=[1;diff(y(:,sc2))];
binoc=sort(unique([binind(jumps>0),binoc]));

%% n-1 dimensional pdf
% w=zeros(d+1,r+1);
nw=cell(1,numel(unique(y(:,1))));
parfor h=1:numel(unique(y(:,1)))
    b1=binoc(h);
    b2=binoc(h+1);
    tz=y(b1:b2-1,:);
    if b2-b1>1
        wh=zeros(c_1+1,b2-b1+1);
        co=1;
        for k=1:b2-b1
            zv=tz(k,1:c_1);
            v=all(bsxfun(@eq,wh(1:c_1,1:co),zv'));
            %v=sum(bsxfun(@eq,wh(1:c_1,1:co),zv'),1)==c_1;
            if any(v)
                wh(c_1+1,v)=wh(c_1+1,v)+1;
            else
                co=co+1;
                wh(:,co)=[zv 1];
            end
        end
        %wh(:,[1 co+1:b2-b1+1])=[];
        %nw{h}=wh;
        nw{h}=wh(:,2:co);
    else
        nw{h}=[tz 1]';
    end
end
nw=cell2mat(nw);
nw(c_1+1,:)=nw(c_1+1,:)/r;
%% n dimensional pdf ?
% % d=c;
% % 
% % w=zeros(d+1,r+1);
% % co=1;
% % for k=1:r
% %     zv=z(k,1:d);
% %     %         v=ismember(w(:,1:d),zv,'rows');
% %     %         v=sum(w(1:co,1:d)==repmat(zv,co,1),2)==d;
% %     v=sum(bsxfun(@eq,w(1:d,1:co),zv'),1)==d;
% %     if any(v)
% %         w(d+1,v)=w(d+1,v)+1;
% %     else
% %         co=co+1;
% %         w(:,co)=[zv 1];
% %     end
% % %     fprintf('.')
% % end
% % % fprintf('\n');
% % 
% % w(:,[1 co+1:r+1])=[];
% % save('mi_out.mat','w','-v7.3')
%% n-1 dimensional pdf
% % b=1;
% % ci=setdiff(1:c,b);
% % 
% % nw=zeros(d,size(w,2)+1);
% % co=1;
% % for k=1:size(w,2)
% %     wv=w(ci,k);
% %     v=sum(bsxfun(@eq,nw(1:d-1,1:co),wv),1)==d-1;
% %     if any(v)
% %         nw(d,v)=nw(d,v)+w(d+1,k);
% %     else
% %         co=co+1;
% %         nw(:,co)=[wv w(d+1,k)];
% %     end
% %     
% % end
% % nw(:,[1 co+1:size(w,2)+1])=[];
% % save('mi_out.mat','nw','-append')
%% 1 dimensional pdf
fprintf('start 1-dim\n');
pb = hist(ics(:,b),nbin(b))/r;

%% MI
fprintf('start MI\n');

ci=setdiff(1:c,b);
[~,i]=sort(w(ci(1),:));
w=w(:,i);
[~,i]=sort(nw(1,:));
nw=nw(:,i);

%find first occur
iw=[1 diff(w(ci(1),:))]';
binocw=1:size(w,2);
binocw=[binocw(iw>0),size(w,2)+1];    

inw=[1 diff(nw(1,:))]';
binocnw=1:size(nw,2);
binocnw=[binocnw(inw>0),size(nw,2)+1];   

% MI=zeros(1,size(w,2));
MI=cell(1,numel(unique(w(ci(1),:))));
parfor h=1:numel(unique(w(ci(1),:)))
    bw1=binocw(h);
    bw2=binocw(h+1);
    bnw1=binocnw(h);
    bnw2=binocnw(h+1);
    
    tMI=zeros(1,bw2-bw1);
    tw=w(:,bw1:bw2-1);
    tnw=nw(:,bnw1:bnw2-1);
    for k = 1:bw2-bw1
        k
        p_b = pb(tw(b,k));
        p_n = tw(c+1,k);
        p_n_1 = sum(tnw(c,sum(bsxfun(@eq,tnw(1:c-1,:),tw(ci,k)),1)==c));
        tMI(k) = p_n * log2(p_n/(p_b*p_n_1));
    end
    MI{h}=tMI;
end
MI=sum(cell2mat(MI));
%% end
fprintf('done\n');