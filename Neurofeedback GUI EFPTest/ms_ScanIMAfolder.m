function [series] = ms_ScanIMAfolder(P)
%% scans a folder and returns the IMAs info in a 'series' structure; 
% series is empty if no IMAs were found

start = pwd;
if nargin<1
    P = spm_select(1, 'dir', 'Select IMA dir');
end
cd(P)
series=[];
%% this approach is not so good.. takes too long
% test=spm_dicom_headers(ls('*.IMA'));

%% first.. look if there are IMAs at all
tmp=ls('*.IMA');

if ~isempty(tmp)
    %% check for "series"
    for ix=1:size(tmp,1)
        loc=strfind(tmp(ix,:),'.');
        seriesStr{ix}=tmp(ix,loc(3)+1:loc(4)-1); % series number as string
    end
    allSeries=unique(seriesStr);
    seriesName=tmp(1,1:loc(3));
    
    %% create a structure with all infos
    for ix=1:length(allSeries)
        series(ix).seriesStr = allSeries{ix};
        series(ix).seriesNumber = str2double(series(ix).seriesStr);
        series(ix).IMAs = cellstr(ls([seriesName series(ix).seriesStr '*.IMA']));
        hdr = spm_dicom_headers(series(ix).IMAs{1});
        series(ix).description = hdr{1}.SeriesDescription;
    end
end

cd(start)