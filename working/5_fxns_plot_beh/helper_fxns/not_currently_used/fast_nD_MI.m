function [MI,pMI,nbin]=fast_nD_MI(x,nbin)
% [MI,pMI,nbin]=fast_nD_MI(x,nbin)
%
% x    : input timeseries       (time-by-chans 2-D matrix)
% nbin : # of bins for pdfs     (1-by-1 scalar or 1-by-chans vector, if 
%                                   none passed then Freedman-Diaconis rule
%                                   (1981) used to estimate nbin)
%
% MI   : mutual information in bits     (1-by-chans vector)
% pMI  : normalized mutual information  (1-by-chans vector)
% nbin : # of bins for pdfs             (especially useful if no nbin 
%                                           argument passed in initial 
%                                           function call)
%
%
% Thomas Wozny & Zach Jessen, last updated 11/02/22



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% orient data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r,c]=size(x);
if c>r
    x=x';
    [r,c]=size(x);
end
invr=1/r;

%%%%%%%%%%%%%%%%%%%%%%%%%%% establish bin number %%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('nbin','var') || isempty(nbin)%Freedman-Diaconis rule (1981)
    nbin=ceil(range(x)./(2*iqr(x)*r.^(-1/3)));
    nbin=repmat(ceil(mean(nbin)),1,c);
elseif length(nbin)~=c
    nbin=repmat(nbin(1),1,c);
end

%%%%%%%%%%%%%%%%%%%%%%% generate histogram bin IDs %%%%%%%%%%%%%%%%%%%%%%%%
% %normalize raw values to assign bin IDs
% x=bsxfun(@minus,x,min(x));
% x=ceil(bsxfun(@times,x,nbin./max(x)));
%ensure bin IDs are not greater than max bin number
for n=1:c 
    x(x(:,n)>nbin(n),n)=nbin(n); 
end
%ensure bin IDs are not less than min bin number (1)
x(~x)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N dim pdf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,ii,J]=unique(x,'rows','stable');
%updated with histcounts, formerly: W=hist(J,1:max(J))'.*invr;
W=histcounts(J,max(J))'.*invr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N-1 dim pdfs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NW=zeros(size(W,1),c);
parfor n=1:c
    [~,~,tNj]=unique(x(:,setdiff(1:c,n)),'rows','stable');
    %updated with histcounts, formerly: h=hist(tNj,1:max(tNj))'.*invr;
    h=histcounts(tNj,max(tNj))'.*invr;
    h=h(tNj);
    NW(:,n)=h(ii);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 dim pdfs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MW=zeros(size(W,1),c);
H=zeros(1,c);
parfor n=1:c
    %updated with histcounts, formerly: tMW=hist(x(:,n),nbin(n)).*(invr);
    tMW=histcounts(x(:,n),nbin(n)).*(invr);
    MW(:,n)=tMW(x(ii,n));
    H(n)=-nansum(tMW.*log2(tMW));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Mutual information %%%%%%%%%%%%%%%%%%%%%%%%%%%%
MI=sum(bsxfun(@times,log2(bsxfun(@times,1./(MW.*NW),W)),W));
pMI=MI./H;

end