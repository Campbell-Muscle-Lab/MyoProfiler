function test_spline

profile = load('profile.mat','-mat');

profile = profile.im_profile';

profile = rescale(profile)

[~,peaks] = findpeaks(profile,'MinPeakDistance',25);
[~,dips] = findpeaks(-profile+1,'MinPeakDistance',50);

se = strel('disk', 200);
im_open = (imopen(profile,se));

% number_of_elements = numel(x);
% fraction = app.FractionSpinner.Value/100;
% fraction_ix = ceil(fraction*number_of_elements);
% app.c_x = [];
% app.c_y = [];
% ix_1 = [];
% ix_2 = [];
% ix_1 = 1:fraction_ix;
% ix_2 = number_of_elements:-1:(number_of_elements - fraction_ix + 1);
% app.c_x = [ix_1 ix_2];
% app.c_y = x(app.c_x);
% s = app.SmoothingEditField.Value;
% app.gel_data.background(box_no).x_back = csaps(app.c_x,app.c_y,s,0:numel(x)-1)'

figure(1)
clf
subplot(2,1,1)
hold on
plot(profile,'k')
plot(peaks,profile(peaks),'ro','MarkerSize',10)
plot(dips,profile(dips),'ms','MarkerSize',10)

for i = 1 : numel(dips) - 1

cx = [dips(i) dips(i+1)];
cy = profile(cx);

spl = csaps(cx,cy,0.5,cx(1):cx(end))'

plot(cx(1):cx(end),spl)

end






subplot(2,1,2)
hold on
plot(profile,'k')
plot(im_open,'r')





end