function ken_test_open

where_is_matlab_utilities

nd2_file_string = 'demo_nikon.nd2'
nd2_file = bfopen(nd2_file_string)
meta_data = nd2_file{1,2}; % this is a hashtable a java variant of a structure
im = nd2_file{1,1} % these are the images. they accompanied with planes and channels. Nicole's should have 3 probably.

sp = layout_subplots(panels_high = 1, ...
                panels_wide = 2)

colormap('gray')
for i = 1 : size(im,2)

    temp_imv = [];
    temp_im = im{i,1};
    
    subplot(sp(i))
    center_image_with_preserved_aspect_ratio(temp_im);


end

exportgraphics(gcf,'ken_test_open.png')
end