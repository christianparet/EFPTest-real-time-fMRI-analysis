function [outr, res]=SNiP_tbxvol_extract_fast(srcimgs, roispec, avg)

% functionm [outr, res]=SNiP_tbxvol_extract_fast(srcimages, roispec, avg)
%
% Reads raw data from a series of images (nifti only) and extracts from 
% a given ROI (mask, sphere, single coordinate). 
%
% Example code:
% F = spm_select(Inf,'image');
% 
% % roi (mask)
% roispec{1}.srcimg = spm_select(1, 'image');
% 
% % sphere around center
% roispec{2}.roisphere.roicent = [0 0 0]';
% roispec{2}.roisphere.roirad = 6;
% 
% % coordinate
% roispec{3}.roilist = [12 12 12]';
% 
% [outr, res] = SNiP_tbxvol_extract_fast(F, roispec, 'none',0);
%
% scrcimages    char array of input files
% roispec       roi structure (see example above)
% avg           none or vox
%
% outr          the output format depends on option avg
%       none    outr is a #rois cell array of matrices (tpt x voxels in
%               ROI)
%       vox     outr is a tpt x ROI array with meaned time series
%
% This code adapted from the Volumes toolbox for SPM.
%
% Axel Schäfer, axel.schaefer@zi-mannheim.de, March 16 2013
% Systems Neuroscience in Psychiatry (SNiP), CIMH, Mannheim


interp=0;
res = struct('raw',[], 'adj',[], 'posmm',[], 'Vspace',[]);

V = spm_vol(srcimgs);
for l = 1:numel(roispec)
    fprintf(1,'roi # %d\n',l);
    res(l).Vspace = rmfield(V(1),{'fname','private'});
    [res(l).posmm res(l).posvx] = get_pos(roispec{l},V(1));
    for k=1:numel(V)
        raw = spm_sample_vol(V(k), res(l).posvx(1,:), ...
            res(l).posvx(2,:), res(l).posvx(3,:), ...
            interp);
        switch avg
            case 'none'
                outr{l}(k,:) = raw;
            case 'vox'
                outr(k,l) = mean(raw);
        end;
    end;
end;


function [posmm, posvx] = get_pos(roispec, V)

if isfield(roispec, 'srcimg')
    % resample mask VM in space of current image V
    VM = spm_vol(roispec.srcimg);
    x = []; y = []; z = [];
    [x1 y1] = ndgrid(1:V.dim(1),1:V.dim(2));
    for p = 1:V.dim(3)
        B = spm_matrix([0 0 -p 0 0 0 1 1 1]);
        M = VM.mat\(V.mat/B);
        msk = find(spm_slice_vol(VM,M,V.dim(1:2),0));
        if ~isempty(msk)
            z1 = p*ones(size(msk(:)));
            x = [x; x1(msk(:))];
            y = [y; y1(msk(:))];
            z = [z; z1];
        end;
    end;
    posvx = [x'; y'; z'];
    xyzmm = V.mat*[posvx;ones(1,size(posvx,2))];
    posmm = xyzmm(1:3,:);
elseif isfield(roispec, 'roisphere')
    cent = round(V.mat\[roispec.roisphere.roicent; 1]);
    tmp = spm_imatrix(V.mat);
    vdim = tmp(7:9);
    vxrad = ceil((roispec.roisphere.roirad*ones(1,3))./ ...
        vdim)';
    [x y z] = ndgrid(-vxrad(1):sign(vdim(1)):vxrad(1), ...
        -vxrad(2):sign(vdim(2)):vxrad(2), ...
        -vxrad(3):sign(vdim(3)):vxrad(3));
    sel = (x./vxrad(1)).^2 + (y./vxrad(2)).^2 + (z./vxrad(3)).^2 <= 1;
    x = cent(1)+x(sel(:));
    y = cent(2)+y(sel(:));
    z = cent(3)+z(sel(:));
    posvx = [x y z]';
    xyzmm = V.mat*[posvx;ones(1,size(posvx,2))];
    posmm = xyzmm(1:3,:);
elseif isfield(roispec, 'roilist')
    posmm = roispec.roilist;
    posvx = V.mat\[posmm; ones(1,size(posmm,2))];
    posvx = posvx(1:3,:);
end;
