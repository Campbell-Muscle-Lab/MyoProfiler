function test_binary



I = imread("myomesin.png");

Img = imclose(I,strel('line',50,45));

B = fibermetric(I,2,ObjectPolarity='bright');

B = imbinarize(B,'adaptive');

thresh_mask = imbinarize(Img);
thresh_mask = imerode(thresh_mask,strel('disk',5));

B(~thresh_mask) = 0;
% B = bwareaopen(B,100);
B = imfill(B,'holes');
imshowpair(thresh_mask,B,'montage')



end