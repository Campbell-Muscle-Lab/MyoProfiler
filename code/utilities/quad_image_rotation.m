function [im_profile_rot,rot_prof_x,rot_prof_y] = quad_image_rotation(im,roi_width,prof_x,prof_y,x_profile)

f = waitbar(0,'1','Name','ROI Expansion');
for i = 1 : (numel(prof_x) - 1)
    ang(i) = atan2((prof_y(i+1) - prof_y(i)), (prof_x(i+1) - prof_x(i)));
end

for i = 1 : (numel(prof_x)-1)

    sx = 1;
    sy = 1;
    [H, W] = size(im);
    phi = ang(i);
    x0 = prof_x(i);
    y0 = prof_y(i);
    d = round(roi_width/2);

    rot_prof_x(i,1) =  prof_x(i) + roi_width*cos(phi)*sin(phi);
    rot_prof_y(i,1) = prof_y(i) + roi_width * sin(phi)*sin(phi) - roi_width;
    rot_prof_x(i,2) =  prof_x(i) - roi_width*cos(phi)*sin(phi);
    rot_prof_y(i,2) = prof_y(i) - roi_width * sin(phi)*sin(phi) + roi_width;

    T = [sx*cos(phi), -sx*sin(phi), 0
        sy*sin(phi),  sy*cos(phi), 0
        (W+1)/2-((sx*x0*cos(phi))+(sy*y0*sin(phi))), (H+1)/2+((sx*x0*sin(phi))-(sy*y0*cos(phi))), 1];

    tform = affine2d(T);

    im_rot = imwarp(im, tform, 'OutputView', imref2d([H, W]), 'Interp', 'cubic');

    cx = W/2;
    cy = H/2;
    cy_t = cy+d;
    cy_b = cy-d;

    % im_profile_rot(x_profile(i)) = mean(im_rot(cx, (cy_b):(cy_t)));
    im_profile_rot(x_profile(i)) = mean(im_rot((cy_b):(cy_t),cx));

    % figure(99999)
    % % plot(im_profile_rot)


    if (mod(i,10) == 1)

        mes = sprintf('Extracting profile: %i%% Completed',round(100*i/(numel(x_profile)-1)));
        waitbar(i/(numel(x_profile)-1),f,mes);
        % figure(99999)
        % plot(im_profile_rot)

    end

end
delete(f)

end