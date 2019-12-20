function update_orthview_new(fname, roispec)
spmfig=spm_figure('GetWin', 'Graphics');
spm_clf(spmfig);
H = spm_orthviews('Image',fname);

for i=1:numel(roispec)
    col = zeros(1,3);
    col(i) = 1;
    spm_orthviews('AddColouredImage', H, roispec{i}.srcimg, col);
end
spm_orthviews('Redraw');

end