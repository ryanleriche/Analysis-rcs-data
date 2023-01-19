function [entropy,inflation,nbin,pdf]=entropy1D(data,nbin)
% [entropy,inflation,nbin]=entropy1D(data,nbin)
%     calculates entropy of channels in data
%  
% data: [time,channels]
% nbin: # of bins to use in the histogram/PDF (if no arg passed then nbin 
%       is automatically estimated using Freedman-Diaconis rule (1981)
% 
% entropy: entropy (bits) of each channel
% inflation: estimated error (bits) of how inflated entropy estimation is 
%       given the # of samples
% nbin: # of bins actually used for histogram (especially useful if no nbin
%       argument passed in initial function call)
%     
% 
% TAW_062415 (edits 011618)


%orient data
[N_rows,N_cols]=size(data);
if N_cols>N_rows 
    data=data';
    [N_rows,N_cols]=size(data);
end

%get nbin
if ~exist('nbin','var') || isempty(nbin)%Freedman-Diaconis rule (1981)
    nbin=ceil( range(data) ./ (2.*iqr(data).*N_rows.^(-1/3)) );
    nbin=ceil(mean(nbin));
elseif length(nbin)~=N_cols
    nbin=nbin(1);
end

%make PDF
invr=1./N_rows;
pdf=zeros(nbin,N_cols);


parfor i=1:N_cols
    pdf(:,i)=histcounts(data(:,i),nbin)*invr; % histcounts uses bin edges
end

% entropy in bits from PDF (eps is added to avoid log2(0) = -Inf)
entropy=-sum(pdf.*log2(pdf + eps), 'omitnan');

%estimate inflation of entropy estimate given N sample points
inflation=(nbin-1)/(2*N_rows*log(2));

end