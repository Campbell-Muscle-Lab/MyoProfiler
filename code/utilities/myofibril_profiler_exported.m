classdef myofibril_profiler_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        MyoProfilerUIFigure            matlab.ui.Figure
        Menu                           matlab.ui.container.Menu
        LoadImageMenu                  matlab.ui.container.Menu
        NikonMenu                      matlab.ui.container.Menu
        StandardFormatsMenu            matlab.ui.container.Menu
        LoadAnalysisMenu               matlab.ui.container.Menu
        ExportAnalysisMenu             matlab.ui.container.Menu
        ZLineStainCheckBox             matlab.ui.control.CheckBox
        SarcomerePanel                 matlab.ui.container.Panel
        SummaryTable                   matlab.ui.control.Table
        Sarcomeres                     matlab.ui.control.UIAxes
        SarcomereMean                  matlab.ui.control.UIAxes
        ControlsPanel                  matlab.ui.container.Panel
        ProfilePanel                   matlab.ui.container.Panel
        ZLineProminenceEditField       matlab.ui.control.NumericEditField
        ZLineProminenceEditFieldLabel  matlab.ui.control.Label
        PeakDistancepxEditField        matlab.ui.control.NumericEditField
        PeakDistancepxEditFieldLabel   matlab.ui.control.Label
        CalibrationumpxEditField       matlab.ui.control.NumericEditField
        CalibrationumpxEditFieldLabel  matlab.ui.control.Label
        ABandProminenceEditField       matlab.ui.control.NumericEditField
        ABandProminenceEditFieldLabel  matlab.ui.control.Label
        ROIPanel                       matlab.ui.container.Panel
        ROIWidthpxEditField            matlab.ui.control.NumericEditField
        ROIWidthpxEditFieldLabel       matlab.ui.control.Label
        SelectPointsButton             matlab.ui.control.Button
        ChannelColormapPanel           matlab.ui.container.Panel
        ChannelColorDropDown           matlab.ui.control.DropDown
        ProfilerPanel                  matlab.ui.container.Panel
        ProfileIntensityYCoord         matlab.ui.control.UIAxes
        ProfileIntensityXCoord         matlab.ui.control.UIAxes
        ProfileIntensity               matlab.ui.control.UIAxes
        ImageDisplayPanel              matlab.ui.container.Panel
        TabGroup                       matlab.ui.container.TabGroup
        ImageTab                       matlab.ui.container.Tab
        ImageAxes                      matlab.ui.control.UIAxes
    end


    properties (Access = public)
        myofibril_data = []
        shaded_stats
    end

    properties (Access = private)
        Tabs
        ChannelAxes
        ChannelColors
        SplineLine
        LoadedAnalysis % Description
        Patches = [] % Description
    end

    methods (Access = public)

        function ExtractProfiles(app,channel_no)

            im = app.myofibril_data.image{channel_no};
            xs = app.myofibril_data.profile(channel_no).xs;
            ys = app.myofibril_data.profile(channel_no).ys;
            prominence = app.ZLineProminenceEditField.Value;
            peak_distance = app.PeakDistancepxEditField.Value;
            roi_width = app.ROIWidthpxEditField.Value;
            col = app.ChannelColors(channel_no,:);

            [prof_x, prof_y, im_profile] = improfile(im, xs, ys);

            x_profile = 1:numel(im_profile);


            if roi_width > 1
                [im_profile_rot,rot_prof_x,rot_prof_y] = quad_image_rotation(im,roi_width,prof_x,prof_y,x_profile);
                x_profile(end) = [];
                prof_x(end) = [];
                prof_y(end) = [];
                im_profile = [];
                im_profile = im_profile_rot;
                
                
                app.GeneratePatches(rot_prof_x,rot_prof_y,channel_no);

            end

            [pks_z_line, locs_z_line] = findpeaks(-rescale(im_profile), ...
                'MinPeakProminence', prominence * range(rescale(im_profile)), ...
                'MinPeakDistance',peak_distance);

            hold(app.ProfileIntensity,'on')
            hold(app.ProfileIntensityXCoord,'on')
            hold(app.ProfileIntensityYCoord,'on')


            plot3(app.ProfileIntensity,prof_x,prof_y,rescale(im_profile),'color',col,'LineWidth',1.7);
            ang = atan2((prof_y(end) - prof_y(1)), (prof_x(end) - prof_x(1)));
            view(app.ProfileIntensity,[rad2deg(ang),30])
            plot(app.ProfileIntensityXCoord,prof_x,rescale(im_profile),'Color',col,'LineWidth',1.7);
            plot(app.ProfileIntensityYCoord,prof_y,rescale(im_profile),'Color',col,'LineWidth',1.7);
            
            x_limits_horizontal = [prof_x(1) prof_x(end)];
            x_limits_horizontal = sort(x_limits_horizontal,'ascend');
            x_limits_vertical = [prof_y(1) prof_y(end)];
            x_limits_vertical = sort(x_limits_vertical,'ascend');
            xlim(app.ProfileIntensityXCoord,x_limits_horizontal)
            xlim(app.ProfileIntensityYCoord,x_limits_vertical)

            col(4) = 0.7;
            for i = 1 : numel(locs_z_line)
                plot(app.ProfileIntensityXCoord,prof_x(locs_z_line((i)))*ones(1,10),linspace(0,1,10),'LineStyle','--','color',col,'LineWidth',1.7)
                plot(app.ProfileIntensityYCoord,prof_y(locs_z_line((i)))*ones(1,10),linspace(0,1,10),'LineStyle','--','color',col,'LineWidth',1.7)
            end
            col(4) = [];

            app.ExtractSarcomeres(channel_no,im_profile,locs_z_line,prof_x,prof_y,col);

        end

        function ExtractSarcomeres(app,channel_no,im_profile,locs_z_line,prof_x,prof_y,col)

            prominence = app.ABandProminenceEditField.Value;
            peak_distance = app.PeakDistancepxEditField.Value;
            calibration = app.CalibrationumpxEditField.Value;
            sarcs_to_remove = [];

            no_of_sarcomeres = numel(locs_z_line) - 1;

            for i = 1 : no_of_sarcomeres
                distance = arclength(prof_x(locs_z_line(i):locs_z_line(i+1)),prof_y(locs_z_line(i):locs_z_line(i+1)),'spline');
                sarc_len(i) = calibration*distance;              
            end

            % for i = 1 : no_of_sarcomeres
            %     profile_indices = locs_z_line(i) : locs_z_line(i+1);
            %     x_coord = prof_x(profile_indices);
            %     y_coord = prof_y(profile_indices);
            %     sarc_profile = im_profile(profile_indices);
            %     sarc_profile = rescale(sarc_profile)
            %     locs_all = [];
            %     locs_m_line = [];
            % 
            %     % u_fwhm = [];
            %     try
            %         [~, locs_all] = findpeaks((sarc_profile), ...
            %             'MinPeakDistance',peak_distance, ...
            %             'MinPeakProminence',prominence);
            %         [~,locs_m_line] = max(-sarc_profile(locs_all(1):locs_all(2)));
            %         locs_m_line = locs_m_line + locs_all(1) - 1;
            %     catch
            %         sarcs_to_remove = [sarcs_to_remove;i];
            %         u_fwhm(i) = NaN;
            %         continue
            %     end
            % 
            %     u_fwhm_ix(i,1) = find(sarc_profile >= 0.5*sarc_profile(locs_all(1)),1,'first');
            %     u_fwhm_ix(i,2) = find(sarc_profile >= 0.5*sarc_profile(locs_all(2)),1,'last');
            %     u_fwhm(i) = calibration*arclength(x_coord(u_fwhm_ix(i,1):u_fwhm_ix(i,2)),y_coord(u_fwhm_ix(i,1):u_fwhm_ix(i,2)),'spline')
            % 
            % end

            sarc_col = return_color_scheme(col,no_of_sarcomeres);

            for i = 1 : no_of_sarcomeres

                mean_sarc_len = mean(sarc_len);
                x_sarc_profile(i,:) = linspace(0,sarc_len(i), 1000);
                profile_indices = locs_z_line(i) : locs_z_line(i+1);
                x_temp = sarc_len(i)*normalize(profile_indices, 'range');
                x_sarc_profile(i,end) = x_temp(end);
                y_sarc_profile(i, :) = interp1(x_temp, im_profile(profile_indices), ...
                    x_sarc_profile(i,:));
                y_sarc_profile(i,:) = normalize(y_sarc_profile(i,:),'range');

                locs_all = [];
                locs_m_line = [];
                % fwhm = [];
                try
                    [~, locs_all] = findpeaks(y_sarc_profile(i,:), ...
                        'MinPeakDistance',peak_distance, ...
                        'MinPeakProminence',prominence);
                    [~,locs_m_line] = max(-y_sarc_profile(i,locs_all(1):locs_all(2)));
                    locs_m_line = locs_m_line + locs_all(1) - 1;
                catch
                    sarcs_to_remove = [sarcs_to_remove;i];
                    fwhm(i) = NaN;
                    continue
                end

                x_sarc_profile(i,:) = x_sarc_profile(i,:)  - x_sarc_profile(i,locs_m_line);

                fwhm_ix(i,1) = find(y_sarc_profile(i,:) >= 0.5*(y_sarc_profile(i,locs_all(1)) + y_sarc_profile(i,1)),1,'first');
                fwhm_ix(i,2) = find(y_sarc_profile(i,:) >= 0.5*(y_sarc_profile(i,locs_all(2)) + y_sarc_profile(i,end)),1,'last');
                fwhm(i) = x_sarc_profile(i,fwhm_ix(i,2)) - x_sarc_profile(i,fwhm_ix(i,1))

                hold(app.Sarcomeres,'on')

                plot(app.Sarcomeres,x_sarc_profile(i,:),y_sarc_profile(i,:),'Color',sarc_col(i,:),'LineWidth',1.7)
                scatter(app.Sarcomeres,x_sarc_profile(i,locs_m_line),y_sarc_profile(i,locs_m_line),80,'d','MarkerEdgeColor',[0 0 0],...
                    'MarkerFaceColor',[1 0 1],...
                    'MarkerFaceAlpha',0.5)
                scatter(app.Sarcomeres,x_sarc_profile(i,locs_all),y_sarc_profile(i,locs_all),80,'d','MarkerEdgeColor',[0 0 0],...
                    'MarkerFaceColor',[1 0 0],...
                    'MarkerFaceAlpha',0.5)

            end
            if ~isempty(sarcs_to_remove)
                y_sarc_profile(sarcs_to_remove,:) = [];
                x_sarc_profile(sarcs_to_remove,:) = [];
                sarc_len(sarcs_to_remove) = [];
                fwhm(sarcs_to_remove) = [];
                % u_fwhm(sarcs_to_remove) = [];
            end
            
            if no_of_sarcomeres > 0
            mean_sarc_profile = mean(y_sarc_profile,1);

            x_sarc_profile_mean = mean(x_sarc_profile,1);
            hold(app.SarcomereMean,"on")
            if size(y_sarc_profile,1) == 1
                plot(app.SarcomereMean, x_sarc_profile_mean,y_sarc_profile,'-','LineWidth',1.7,'Color',col)
            else
                shadedErrorBar2(x_sarc_profile_mean, y_sarc_profile, {@mean, @std},{'-','LineWidth',1.7,'Color',col},1,app.SarcomereMean);
            end

            xlim(app.Sarcomeres,[-0.51*mean_sarc_len 0.51*mean_sarc_len])
            xlim(app.SarcomereMean,[-0.51*mean_sarc_len 0.51*mean_sarc_len])


            app.myofibril_data.profile(channel_no).intensity = im_profile;
            app.myofibril_data.profile(channel_no).sarcomere_intensities = y_sarc_profile;
            app.myofibril_data.profile(channel_no).sarcomere_location = x_sarc_profile;
            app.myofibril_data.profile(channel_no).mean_sarcomere_location = mean(x_sarc_profile,1);
            app.myofibril_data.profile(channel_no).mean_sarcomere_profile = mean_sarc_profile;
            app.myofibril_data.profile(channel_no).std_sarcomere_profile = std(y_sarc_profile,0,1);
            app.myofibril_data.profile(channel_no).sem_sarcomere_profile = std(y_sarc_profile,0,1)./sqrt(length(y_sarc_profile));
            app.myofibril_data.profile(channel_no).sarcomere_lengths = sarc_len;
            app.myofibril_data.profile(channel_no).fwhm = fwhm;
            app.myofibril_data.calibration = calibration;
            
            no_of_sarcomeres = numel(sarc_len);
            if no_of_sarcomeres > 0
                app.UpdateSummaryTable(channel_no,no_of_sarcomeres,sarc_col);
            end
            end

        end

        function RefreshDisplay(app)

            ax = {'ProfileIntensity','Sarcomeres','SarcomereMean','ProfileIntensityXCoord','ProfileIntensityYCoord'};

            for i = 1:numel(ax)
                cla(app.(ax{i}))
            end

            app.SummaryTable.Data = [];
        end

        function ColorScheme(app)

            app.ChannelColors = [];
            
            selected_scheme = app.ChannelColorDropDown.Value;

            num_of_channels = size(app.Tabs,2)-1;

            switch selected_scheme
                case 'Em. Wavelength'
                    for i = 1 : num_of_channels
                        app.ChannelColors(i,:) = wavelength2color(app.myofibril_data.em_wavelengths(i));
                    end
                case 'Parula'
                    app.ChannelColors = parula(num_of_channels);
            end

        end

        function PseudoColoring(app)

            raw_images = app.myofibril_data.image;

            for i = 1 : numel(raw_images)

                r_channel = app.ChannelColors(i,1) * raw_images{i,1};
                g_channel = app.ChannelColors(i,2) * raw_images{i,1};
                b_channel = app.ChannelColors(i,3) * raw_images{i,1};
                pseudo_color_images{i,1} = cat(3,r_channel,g_channel,b_channel);
                pseudo_color_images{i,1} = rescale(pseudo_color_images{i,1});

                if i == 1
                    fused_image = pseudo_color_images{i,1};
                else
                    fused_image = imfuse(fused_image,pseudo_color_images{i,1},'blend', 'Scaling', 'joint');
                end
            end
            app.myofibril_data.pseudo_color_images = pseudo_color_images;
            app.myofibril_data.fused_image = fused_image;
            center_image_with_preserved_aspect_ratio(app.myofibril_data.fused_image,app.ChannelAxes{end})


        end


        function UpdateSummaryTable(app,channel_no,no_of_sarcomeres,sarc_col)

            for i = 1 : no_of_sarcomeres
                st.channel(i,:) = channel_no;
                st.band_no(i,:) = i;
                st.color{i,:} = '';
            end
            st.sarcomere_lengths = app.myofibril_data.profile(channel_no).sarcomere_lengths';
            st.fwhm = app.myofibril_data.profile(channel_no).fwhm';

            if numel(st.sarcomere_lengths) ~= numel(st.fwhm)
                st.fwhm(end+1:numel(st.sarcomere_lengths),1) = NaN;
            end


            app.SummaryTable.Data = [app.SummaryTable.Data; struct2table(st)];

            if channel_no > 1
                starting_ix = size(app.SummaryTable.Data,1) - no_of_sarcomeres;
            else
                starting_ix = 0;
            end

            for i = 1 : no_of_sarcomeres
                s = uistyle("BackgroundColor",sarc_col(i,:));
                addStyle(app.SummaryTable,s,"cell",[starting_ix+i 3])
            end

        end
        
        
        function GeneratePatches(app,rot_prof_x,rot_prof_y,channel_no)
            app.Patches{channel_no} = patch(app.ChannelAxes{channel_no},[rot_prof_x(:,1); flip(rot_prof_x(:,2))], ...
                [rot_prof_y(:,1); flip(rot_prof_y(:,2))], ...
                [1 0 1]*0.8, 'EdgeColor','none', ...
                'FaceAlpha',0.25);
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            addpath(genpath('utilities'))
            movegui(app.MyoProfilerUIFigure,'center')
            colormap(app.ImageAxes, 'gray');
        end

        % Menu selected function: NikonMenu
        function NikonImageSelected(app, event)
            f = figure('Renderer', 'painters', 'Position', [-100 -100 0 0]);

            [file_string,path_string]=uigetfile2( ...
                {'*.nd2','ND2'}, ...
                'Select Image File');
            delete(f)

            if (path_string~=0)
                app.myofibril_data = [];
                app.myofibril_data.image_file_string = fullfile(path_string,file_string);
                im_ax = app.ImageAxes;
                app.Tabs = [];
                app.ChannelAxes = [];

                if contains(file_string,'.nd2')
                    app.myofibril_data.image_file = bfopen(app.myofibril_data.image_file_string);
                    app.myofibril_data.meta_file = app.myofibril_data.image_file{1,2};
                    h = app.myofibril_data.meta_file;
                    app.CalibrationumpxEditField.Value = h.get('Global dCalibration');
                    im = app.myofibril_data.image_file;
                    im = im{1,1}
                    im(:,2) = [];
                    app.myofibril_data.image = im;
                    sz = size(im);
                    delete(app.TabGroup.Children)
                    tab_no = 1;
                    brightfield_im = 0;
                    for i = 1:sz(1)
                        try
                        w_text = sprintf('Global CSU-W1, FilterChanger(EM) #%i',i);
                        w_text = h.get(w_text);
                        em_wavelength = str2double(extractBetween(w_text,"(","/"));
                        app.myofibril_data.em_wavelengths(tab_no) = em_wavelength;
                        app.Tabs{tab_no} = uitab(app.TabGroup,'Title',['Channel ' num2str(tab_no)]);
                        app.ChannelAxes{tab_no} = copyobj(im_ax,app.Tabs{tab_no});
                        app.ChannelAxes{tab_no}.Visible = 'on';
                        center_image_with_preserved_aspect_ratio(app.myofibril_data.image{i,1},app.ChannelAxes{tab_no})
                        % app.myofibril_data.pseudo_color_images{tab_no,1} = app.myofibril_data.image{i,1};
                        tab_no = tab_no + 1;
                        catch
                            brightfield_im = i;
                            continue
                        end

                    end
                    if brightfield_im
                        app.myofibril_data.image(i) = [];
                    end
                    app.Tabs{end+1} = uitab(app.TabGroup,'Title','Merged ');
                    app.ChannelAxes{end+1} = copyobj(im_ax,app.Tabs{end});
                    app.ChannelAxes{end}.Visible = 'on';
                    app.ColorScheme;
                    app.PseudoColoring
                    app.RefreshDisplay
                else
                    im = imread(app.myofibril_data.image_file_string);
                    if (ndims(im)==3)
                        im = rgb2gray(im);
                    end
                    app.myofibril_data.image = im;
                    app.ChannelAxes{1} = copyobj(im_ax,app.ImageTab);
                    app.ChannelAxes{1}.Visible = 'on';
                end

                % center_image_with_preserved_aspect_ratio(app.myofibril_data.image,app.ImageAxes)

            end
        end

        % Menu selected function: StandardFormatsMenu
        function StandardFormatsImageSelected(app, event)
            prompt = {'Total Number of Channels:','Emission Wavelengths (nm):'};
            dlgtitle = 'Standard Formats Channel Input';
            fieldsize = [1 65; 1 65];
            definput = {'',''};
            user_input = inputdlg(prompt,dlgtitle,fieldsize,definput);
            
            if ~isempty(user_input)
                app.myofibril_data = [];
                total_number_of_channels = str2double(user_input{1});
                em_wavelengths = user_input{2};
                em_wavelengths = split(em_wavelengths,",");
                em_wavelengths = str2double(em_wavelengths)';
                if isnan(em_wavelengths)
                    app.ChannelColorDropDown.Value = 'Parula';
                else
                    app.myofibril_data.em_wavelengths = em_wavelengths;
                end

                for loaded_files = 1 : total_number_of_channels
                    f = figure('Renderer', 'painters', 'Position', [-100 -100 0 0]);
                    dialog_text = sprintf('Select Image File for Channel %i',loaded_files);
                    [file_string,path_string]=uigetfile2( ...
                        {'*.tif','TIF';'*.tiff','TIFF';'*.png','PNG'}, ...
                        dialog_text);
                    delete(f)
                    if path_string~=0
                       app.myofibril_data.image_file_string{loaded_files} = fullfile(path_string,file_string);
                       app.myofibril_data.image{loaded_files,1} = imread(app.myofibril_data.image_file_string{loaded_files});
                       if (ndims(app.myofibril_data.image{loaded_files,1})==3)
                           app.myofibril_data.image{loaded_files,1} = rgb2gray(app.myofibril_data.image{loaded_files,1});
                       end
                    end
                end

                if isfield(app.myofibril_data,'image_file_string')

                im_ax = app.ImageAxes;
                app.Tabs = [];
                app.ChannelAxes = [];
                delete(app.TabGroup.Children)
                for i = 1:numel(app.myofibril_data.image)
                    app.Tabs{i} = uitab(app.TabGroup,'Title',['Channel ' num2str(i)]);
                    app.ChannelAxes{i} = copyobj(im_ax,app.Tabs{i});
                    app.ChannelAxes{i}.Visible = 'on';
                    center_image_with_preserved_aspect_ratio(app.myofibril_data.image{i},app.ChannelAxes{i})
                    % app.myofibril_data.pseudo_color_images{i,1} = app.myofibril_data.image{i};
                end
                app.Tabs{end+1} = uitab(app.TabGroup,'Title','Merged ');
                app.ChannelAxes{end+1} = copyobj(im_ax,app.Tabs{end});
                app.ChannelAxes{end}.Visible = 'on';
                app.ColorScheme;
                app.PseudoColoring;
                app.RefreshDisplay;
                end


                % if (path_string~=0)
                %     app.myofibril_data = [];
                %     app.myofibril_data.image_file_string{i} = fullfile(path_string,file_string);
                %     if contains(file_string,'.nd2')
                %         app.myofibril_data.image_file = bfopen(app.myofibril_data.image_file_string);
                %         app.myofibril_data.meta_file = app.myofibril_data.image_file{1,2};
                %         h = app.myofibril_data.meta_file;
                %         app.CalibrationumpxEditField.Value = h.get('Global dCalibration');
                %         im = app.myofibril_data.image_file;
                %         im = im{1,1};
                %         im(:,2) = [];
                %         app.myofibril_data.image = im;
                %         sz = size(im);
                %         delete(app.TabGroup.Children)
                %         for i = 1:sz(1)
                %             app.Tabs{i} = uitab(app.TabGroup,'Title',['Channel ' num2str(i)]);
                %             app.ChannelAxes{i} = copyobj(im_ax,app.Tabs{i});
                %             app.ChannelAxes{i}.Visible = 'on';
                %             center_image_with_preserved_aspect_ratio(app.myofibril_data.image{i,1},app.ChannelAxes{i})
                %             app.myofibril_data.pseudo_color_images{i,1} = app.myofibril_data.image{i,1};
                %             w_text = sprintf('Global CSU-W1, FilterChanger(EM) #%i',i);
                %             w_text = h.get(w_text);
                %             em_wavelength = str2double(extractBetween(w_text,"(","/"));
                %             app.myofibril_data.em_wavelengths(i) = em_wavelength;
                %         end
                %         app.Tabs{end+1} = uitab(app.TabGroup,'Title','Merged ');
                %         app.ChannelAxes{end+1} = copyobj(im_ax,app.Tabs{end});
                %         app.ChannelAxes{end}.Visible = 'on';
                %         app.ColorScheme;
                %         app.PseudoColoring
                %         app.RefreshDisplay
                %     else
                %         im = imread(app.myofibril_data.image_file_string);
                %         app.myofibril_data.image = im;
                %         app.ChannelAxes{1} = copyobj(im_ax,app.ImageTab);
                %         app.ChannelAxes{1}.Visible = 'on';
                %     end
                % 
                %     % center_image_with_preserved_aspect_ratio(app.myofibril_data.image,app.ImageAxes)
                % 
                % end
            end
        end

        % Button pushed function: SelectPointsButton
        function SelectPointsButtonPushed(app, event)

            selected_tab = app.TabGroup.SelectedTab.Title;
            patches = app.Patches;
            try
                image_axis_index = str2double(regexp(selected_tab,'\d*','Match'));
                im_axis = app.ChannelAxes{image_axis_index};
                tab_number_excluding_merged = (numel(app.TabGroup.Children) - 1);

                tabs_to_copy = 1 : tab_number_excluding_merged+1;

            catch
                return
            end

            if  isfield(app.myofibril_data,'profile')
                app.myofibril_data.profile = [];

                for i = 1 : numel(tabs_to_copy)
                    cla(app.ChannelAxes{i})
                    if i ~= tabs_to_copy(end)
                        center_image_with_preserved_aspect_ratio(app.myofibril_data.image{i,1},app.ChannelAxes{i})
                    else
                        center_image_with_preserved_aspect_ratio(app.myofibril_data.fused_image,app.ChannelAxes{i})
                    end
                end
                patches = [];
                app.RefreshDisplay
                app.Patches = patches;
                if ~isempty(app.Patches)
                    for i = 1 : numel(app.Patches)
                        app.Patches{i}.FaceAlpha = 0;
                    end
                end
            end

            tabs_to_copy(image_axis_index) = [];



            app.myofibril_data.roi{image_axis_index} = drawpolyline(im_axis,'LineWidth',1E-32,'MarkerSize',6);


            x = app.myofibril_data.roi{image_axis_index}.Position(:,1);
            y = app.myofibril_data.roi{image_axis_index}.Position(:,2);
            cs = spline(x, y);
            app.myofibril_data.profile(image_axis_index).xs = linspace(x(1), x(end), 1000);
            app.myofibril_data.profile(image_axis_index).ys = ppval(cs, app.myofibril_data.profile(image_axis_index).xs);

            hold(im_axis,'on')
            app.SplineLine = cell(1,numel(tabs_to_copy));
            app.SplineLine{image_axis_index} = plot(im_axis,app.myofibril_data.profile(image_axis_index).xs,app.myofibril_data.profile(image_axis_index).ys,"Color",'m','LineWidth',2);

            addlistener(app.myofibril_data.roi{image_axis_index},"ROIMoved",@(src,evt) update_profile(evt,app));

            for i = 1 : numel(tabs_to_copy)
                app.myofibril_data.roi{tabs_to_copy(i)} = copyobj(app.myofibril_data.roi{image_axis_index}, app.ChannelAxes{tabs_to_copy(i)});
                addlistener(app.myofibril_data.roi{tabs_to_copy(i)},"ROIMoved",@(src,evt) update_profile(evt,app));
                app.SplineLine{tabs_to_copy(i)} = copyobj(app.SplineLine{image_axis_index}, app.ChannelAxes{tabs_to_copy(i)});


            end

            for channel_no = 1 : tab_number_excluding_merged
                app.myofibril_data.profile(channel_no).xs = app.myofibril_data.profile(image_axis_index).xs;
                app.myofibril_data.profile(channel_no).ys = app.myofibril_data.profile(image_axis_index).ys;
                app.ExtractProfiles(channel_no)
            end


            function update_profile(evt,app)
                x_upt = evt.CurrentPosition(:,1);
                y_upt = evt.CurrentPosition(:,2);
                cs_upt = spline(x_upt, y_upt);
                xs_upt = linspace(x_upt(1), x_upt(end), 1000);
                ys_upt = ppval(cs_upt, xs_upt);
                app.myofibril_data.profile(1).xs = linspace(x_upt(1), x_upt(end), 1000);
                app.myofibril_data.profile(1).ys = ppval(cs_upt, app.myofibril_data.profile(1).xs);

                for spline_count = 1 : size(app.SplineLine,2)
                    app.SplineLine{spline_count}.XData = app.myofibril_data.profile(1).xs;
                    app.SplineLine{spline_count}.YData = app.myofibril_data.profile(1).ys;
                    app.myofibril_data.roi{spline_count}.Position(:,1) = x_upt;
                    app.myofibril_data.roi{spline_count}.Position(:,2) = y_upt;
                end

                app.RefreshDisplay;
                if ~isempty(app.Patches)
                    for patch_no = 1 : size(app.Tabs,2)
                        app.Patches{patch_no}.FaceAlpha = 0;
                    end
                end
                for ch_no = 1 : size(app.Tabs,2)-1
                    app.myofibril_data.profile(ch_no).xs = app.myofibril_data.profile(1).xs;
                    app.myofibril_data.profile(ch_no).ys = app.myofibril_data.profile(1).ys;
                    app.ExtractProfiles(ch_no)
                end
            end

        end

        % Value changed function: ZLineProminenceEditField
        function ZLineProminenceEditFieldValueChanged(app, event)
            app.RefreshDisplay;
            if isfield(app.myofibril_data,'profile')
                for channel_no = 1 : size(app.Tabs,2)-1
                    app.ExtractProfiles(channel_no)
                end
            end
        end

        % Value changed function: ABandProminenceEditField
        function ABandProminenceEditFieldValueChanged(app, event)
            app.RefreshDisplay;
            if isfield(app.myofibril_data,'profile')
                for channel_no = 1 : size(app.Tabs,2)-1
                    app.ExtractProfiles(channel_no)
                end
            end
        end

        % Value changed function: PeakDistancepxEditField
        function PeakDistancepxEditFieldValueChanged(app, event)
            app.RefreshDisplay;
            if isfield(app.myofibril_data,'profile')
                for channel_no = 1 : size(app.Tabs,2)-1
                    app.ExtractProfiles(channel_no)
                end
            end
        end

        % Value changed function: CalibrationumpxEditField
        function CalibrationumpxEditFieldValueChanged(app, event)
            app.RefreshDisplay;
            if isfield(app.myofibril_data,'profile')
                for channel_no = 1 : size(app.Tabs,2)-1
                    app.ExtractProfiles(channel_no)
                end
            end
        end

        % Value changed function: ChannelColorDropDown
        function ChannelColorDropDownValueChanged(app, event)
            app.ColorScheme;
            app.PseudoColoring;
            app.RefreshDisplay;
            if isfield(app.myofibril_data,'roi')
                app.myofibril_data.roi{size(app.Tabs,2)} = copyobj(app.myofibril_data.roi{1}, app.ChannelAxes{size(app.Tabs,2)});
                addlistener(app.myofibril_data.roi{size(app.Tabs,2)},"ROIMoved",@(src,evt) update_profile_3(evt,app));
                app.SplineLine{size(app.Tabs,2)} = copyobj(app.SplineLine{1}, app.ChannelAxes{size(app.Tabs,2)});
                if isfield(app.myofibril_data,'profile')
                    for channel_no = 1 : size(app.Tabs,2)-1
                        app.ExtractProfiles(channel_no)
                    end
                end
            end
            function update_profile_3(evt,app)
                x_upt = evt.CurrentPosition(:,1);
                y_upt = evt.CurrentPosition(:,2);
                cs_upt = spline(x_upt, y_upt);
                xs_upt = linspace(x_upt(1), x_upt(end), 1000);
                ys_upt = ppval(cs_upt, xs_upt);
                app.myofibril_data.profile(1).xs = linspace(x_upt(1), x_upt(end), 1000);
                app.myofibril_data.profile(1).ys = ppval(cs_upt, app.myofibril_data.profile(1).xs);

                for spline_count = 1 : size(app.SplineLine,2)
                    app.SplineLine{spline_count}.XData = app.myofibril_data.profile(1).xs;
                    app.SplineLine{spline_count}.YData = app.myofibril_data.profile(1).ys;
                    app.myofibril_data.roi{spline_count}.Position(:,1) = x_upt;
                    app.myofibril_data.roi{spline_count}.Position(:,2) = y_upt;
                end

                app.RefreshDisplay;
                if ~isempty(app.Patches)
                    for patch_no = 1 : size(app.Tabs,2)
                        app.Patches{patch_no}.FaceAlpha = 0;
                    end
                end
                for ch_no = 1 : size(app.Tabs,2)-1
                    app.myofibril_data.profile(ch_no).xs = app.myofibril_data.profile(1).xs;
                    app.myofibril_data.profile(ch_no).ys = app.myofibril_data.profile(1).ys;
                    app.ExtractProfiles(ch_no)
                end
            end
        end

        % Menu selected function: ExportAnalysisMenu
        function ExportAnalysisMenuSelected(app, event)
            [file_string,path_string] = uiputfile2( ...
                {'*.xlsx','Excel file'},'Enter Excel File Name For Analysis Results');

            if (path_string~=0)
                output_file_string = fullfile(path_string,file_string);

                try
                    delete(output_file_string);
                end

                no_of_channels = size(app.myofibril_data.profile,2);

                summary_fields = {'image_file','px_to_um_calibration',...
                    'channel_no',...
                    'mean_sarcomere_length_um','std_sarcomere_length_um',...
                    'sem_sarcomere_length_um','mean_fwhm_um',...
                    'std_fwhm_um','sem_fwhm_um'};

                sarcomere_summary_fields = {'channel_no','sarcomere_index',...
                    'sarcomere_length_um','fwhm_um'};

                for i = 1 : numel(summary_fields)
                    sum_out.(summary_fields{i}) = [];
                end

                for i = 1 : numel(sarcomere_summary_fields)
                    sarcomere_out.(sarcomere_summary_fields{i}) = [];
                end


                dat_type = app.myofibril_data.image_file_string;

                if ischar(dat_type)
                    sum_out.image_file = app.myofibril_data.image_file_string;
                else
                    sum_out.image_file{1} = app.myofibril_data.image_file_string{1};
                end
                sum_out.px_to_um_calibration = num2str(app.myofibril_data.calibration);

                for i = 1 : no_of_channels

                    sum_out.channel_no(i,1) = i;
                    sum_out.mean_sarcomere_length_um(i,1) = mean(app.myofibril_data.profile(i).sarcomere_lengths);
                    sum_out.std_sarcomere_length_um(i,1) = std(app.myofibril_data.profile(i).sarcomere_lengths);
                    sum_out.sem_sarcomere_length_um(i,1) = std(app.myofibril_data.profile(i).sarcomere_lengths)/sqrt(numel(app.myofibril_data.profile(i).sarcomere_lengths));

                    sum_out.mean_fwhm_um(i,1) = mean(app.myofibril_data.profile(i).fwhm);
                    sum_out.std_fwhm_um(i,1) = std(app.myofibril_data.profile(i).fwhm);
                    sum_out.sem_fwhm_um(i,1) = std(app.myofibril_data.profile(i).fwhm)/sqrt(numel(app.myofibril_data.profile(i).fwhm));

                end




                for i = 1 : no_of_channels
                    sarcomere_out.channel_no = [sarcomere_out.channel_no;...
                        i*ones(numel(app.myofibril_data.profile(i).sarcomere_lengths),1)];
                    sarcomere_out.sarcomere_index = [sarcomere_out.sarcomere_index;...
                        (1:numel(app.myofibril_data.profile(i).sarcomere_lengths))'];
                    sarcomere_out.sarcomere_length_um = [sarcomere_out.sarcomere_length_um;...
                        (app.myofibril_data.profile(i).sarcomere_lengths)'];
                    sarcomere_out.fwhm_um = [sarcomere_out.fwhm_um;...
                        (app.myofibril_data.profile(i).fwhm)'];

                end

                if numel(app.myofibril_data.image_file_string) > 1 && ~ischar(dat_type)
                    for i = 2 : no_of_channels
                        sum_out.image_file{i,1} = app.myofibril_data.image_file_string{i};
                    end
                else
                    sum_out.image_file(2:no_of_channels,1) = NaN;
                end
                sum_out.px_to_um_calibration(2:no_of_channels,1) = sum_out.px_to_um_calibration(1);
                writetable(struct2table(sum_out),output_file_string,'Sheet','Analysis Summary')
                writetable(struct2table(sarcomere_out),output_file_string,'Sheet','Sarcomere Summary')


                for i = 1 : no_of_channels
                    channel_name = sprintf('channel_%i',i);
                    sum_sheet_name = sprintf('Channel %i Summary Profiles',i);
                    sheet_name = sprintf('Channel %i Sarcomere Profiles',i);

                    sarcomere.sum_profiles.(channel_name).location_um = app.myofibril_data.profile(i).mean_sarcomere_location';
                    sarcomere.sum_profiles.(channel_name).mean_sarcomere_profile = app.myofibril_data.profile(i).mean_sarcomere_profile';
                    sarcomere.sum_profiles.(channel_name).std_sarcomere_profile = app.myofibril_data.profile(i).std_sarcomere_profile';
                    sarcomere.sum_profiles.(channel_name).sem_sarcomere_profile = app.myofibril_data.profile(i).sem_sarcomere_profile';

                    writetable(struct2table(sarcomere.sum_profiles.(channel_name)),output_file_string,'Sheet',sum_sheet_name)

                    
                    if ~isempty(app.myofibril_data.profile(i).sarcomere_lengths)
                        for j = 1 : numel(app.myofibril_data.profile(i).sarcomere_lengths)
                            var_name = sprintf('sarcomere_intensity_%i',j);
                            location_name = sprintf('location_um_%i',j);
                            sarcomere.(channel_name).(location_name) = app.myofibril_data.profile(i).sarcomere_location(j,:)';
                            sarcomere.(channel_name).(var_name) = app.myofibril_data.profile(i).sarcomere_intensities(j,:)';
                        end

                        writetable(struct2table(sarcomere.(channel_name)),output_file_string,'Sheet',sheet_name)
                    end

                end

                output_file_string = replace(output_file_string,'.xlsx','.myoprof');
                analysis_session = app.myofibril_data;
                analysis_session.roi_width = app.ROIWidthpxEditField.Value;
                analysis_session.z_prominence = app.ZLineProminenceEditField.Value;
                analysis_session.a_prominence = app.ABandProminenceEditField.Value;
                analysis_session.peak_distance = app.PeakDistancepxEditField.Value;

                save(output_file_string,'analysis_session')

            end
        end

        % Value changed function: ROIWidthpxEditField
        function ROIWidthpxEditFieldValueChanged(app, event)
            roi_width = app.ROIWidthpxEditField.Value;
            if  isfield(app.myofibril_data,'profile')
                app.RefreshDisplay;
                if ~isempty(app.Patches)
                    for i = 1 : numel(app.Patches)
                    app.Patches{i}.FaceAlpha = 0;
                    end
                end
                for channel_no = 1 : size(app.Tabs,2)-1
                    app.ExtractProfiles(channel_no)
                end
                if roi_width > 1
                    app.Patches{numel(app.Tabs)} = copyobj(app.Patches{1}, app.ChannelAxes{numel(app.Tabs)});
                end
            end
        end

        % Menu selected function: LoadAnalysisMenu
        function LoadAnalysisMenuSelected(app, event)
            f = figure('Renderer', 'painters', 'Position', [-100 -100 0 0]);

            [file_string,path_string] = uigetfile2( ...
                {'*.myoprof','MyofibrilProfiler file'},'Select MyoProf File To Load Analysis');
            delete(f);

            if (path_string~=0)
                app.SplineLine = [];
                app.Patches = [];

                temp = load(fullfile(path_string,file_string),'-mat','analysis_session');
                analysis_session = temp.analysis_session;
                app.myofibril_data = [];
                app.myofibril_data.image_file_string = analysis_session.image_file_string;
                im_ax = app.ImageAxes;
                app.Tabs = [];
                app.ChannelAxes = [];
                app.RefreshDisplay;
                app.ROIWidthpxEditField.Value = analysis_session.roi_width;
                app.ZLineProminenceEditField.Value = analysis_session.z_prominence;
                app.ABandProminenceEditField.Value = analysis_session.a_prominence ;
                app.PeakDistancepxEditField.Value = analysis_session.peak_distance;
                app.CalibrationumpxEditField.Value = analysis_session.calibration;

                if contains(analysis_session.image_file_string,'.nd2')
                    app.myofibril_data.image_file = analysis_session.image_file;
                    app.myofibril_data.meta_file = app.myofibril_data.image_file{1,2};
                    app.myofibril_data.em_wavelengths = analysis_session.em_wavelengths;
                    h = app.myofibril_data.meta_file;
                    im = analysis_session.image;
                    app.myofibril_data.image = im;
                else
                    app.myofibril_data = deal(analysis_session);
                end

                sz = numel(app.myofibril_data.image);
                delete(app.TabGroup.Children)
                for i = 1:sz
                    app.Tabs{i} = uitab(app.TabGroup,'Title',['Channel ' num2str(i)]);
                    app.ChannelAxes{i} = copyobj(im_ax,app.Tabs{i});
                    app.ChannelAxes{i}.Visible = 'on';
                    center_image_with_preserved_aspect_ratio(app.myofibril_data.image{i,1},app.ChannelAxes{i})
                    % app.myofibril_data.pseudo_color_images{i,1} = app.myofibril_data.image{i,1};
                end
                app.Tabs{end+1} = uitab(app.TabGroup,'Title','Merged ');
                app.ChannelAxes{end+1} = copyobj(im_ax,app.Tabs{end});
                app.ChannelAxes{end}.Visible = 'on';
                app.ColorScheme;
                app.PseudoColoring;
                

                number_of_rois = numel(analysis_session.roi);

                pos = analysis_session.roi{1}.Position;
                app.myofibril_data.roi{1} = drawpolyline(app.ChannelAxes{1},'Position',pos,'LineWidth',1E-32,'MarkerSize',6);

                x = app.myofibril_data.roi{1}.Position(:,1);
                y = app.myofibril_data.roi{1}.Position(:,2);
                cs = spline(x, y);
                app.myofibril_data.profile(1).xs = linspace(x(1), x(end), 1000);
                app.myofibril_data.profile(1).ys = ppval(cs, app.myofibril_data.profile(1).xs);

                hold(app.ChannelAxes{1},'on')
                app.SplineLine = cell(1,numel(app.ChannelAxes));
                app.SplineLine{1} = plot(app.ChannelAxes{1},app.myofibril_data.profile(1).xs,app.myofibril_data.profile(1).ys,"Color",'m','LineWidth',2);

                addlistener(app.myofibril_data.roi{1},"ROIMoved",@(src,evt) update_profile_2(evt,app));

                for i = 2 : numel(app.ChannelAxes)
                    app.myofibril_data.roi{i} = copyobj(app.myofibril_data.roi{1}, app.ChannelAxes{i});
                    addlistener(app.myofibril_data.roi{i},"ROIMoved",@(src,evt) update_profile_2(evt,app));
                    app.SplineLine{i} = copyobj(app.SplineLine{1}, app.ChannelAxes{i});
                end
                
                tab_number_excluding_merged = numel(app.ChannelAxes)-1;

                for channel_no = 1 : tab_number_excluding_merged
                    app.myofibril_data.profile(channel_no).xs = app.myofibril_data.profile(1).xs;
                    app.myofibril_data.profile(channel_no).ys = app.myofibril_data.profile(1).ys;
                    app.ExtractProfiles(channel_no)
                end
            end

            function update_profile_2(evt,app)
                x_upt = evt.CurrentPosition(:,1);
                y_upt = evt.CurrentPosition(:,2);
                cs_upt = spline(x_upt, y_upt);
                xs_upt = linspace(x_upt(1), x_upt(end), 1000);
                ys_upt = ppval(cs_upt, xs_upt);
                app.myofibril_data.profile(1).xs = linspace(x_upt(1), x_upt(end), 1000);
                app.myofibril_data.profile(1).ys = ppval(cs_upt, app.myofibril_data.profile(1).xs);

                for spline_count = 1 : size(app.SplineLine,2)
                    app.SplineLine{spline_count}.XData = app.myofibril_data.profile(1).xs;
                    app.SplineLine{spline_count}.YData = app.myofibril_data.profile(1).ys;
                    app.myofibril_data.roi{spline_count}.Position(:,1) = x_upt;
                    app.myofibril_data.roi{spline_count}.Position(:,2) = y_upt;
                end

                app.RefreshDisplay;
                if ~isempty(app.Patches)
                    for patch_no = 1 : size(app.Tabs,2)
                        app.Patches{patch_no}.FaceAlpha = 0;
                    end
                end
                for ch_no = 1 : size(app.Tabs,2)-1
                    app.myofibril_data.profile(ch_no).xs = app.myofibril_data.profile(1).xs;
                    app.myofibril_data.profile(ch_no).ys = app.myofibril_data.profile(1).ys;
                    app.ExtractProfiles(ch_no)
                end
            end

            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create MyoProfilerUIFigure and hide until all components are created
            app.MyoProfilerUIFigure = uifigure('Visible', 'off');
            app.MyoProfilerUIFigure.Position = [92 92 1336 601];
            app.MyoProfilerUIFigure.Name = 'MyoProfiler';

            % Create Menu
            app.Menu = uimenu(app.MyoProfilerUIFigure);
            app.Menu.Text = 'Menu';

            % Create LoadImageMenu
            app.LoadImageMenu = uimenu(app.Menu);
            app.LoadImageMenu.Text = 'Load Image';

            % Create NikonMenu
            app.NikonMenu = uimenu(app.LoadImageMenu);
            app.NikonMenu.MenuSelectedFcn = createCallbackFcn(app, @NikonImageSelected, true);
            app.NikonMenu.Text = 'Nikon';

            % Create StandardFormatsMenu
            app.StandardFormatsMenu = uimenu(app.LoadImageMenu);
            app.StandardFormatsMenu.MenuSelectedFcn = createCallbackFcn(app, @StandardFormatsImageSelected, true);
            app.StandardFormatsMenu.Text = 'Standard Formats';

            % Create LoadAnalysisMenu
            app.LoadAnalysisMenu = uimenu(app.Menu);
            app.LoadAnalysisMenu.MenuSelectedFcn = createCallbackFcn(app, @LoadAnalysisMenuSelected, true);
            app.LoadAnalysisMenu.Text = 'Load Analysis';

            % Create ExportAnalysisMenu
            app.ExportAnalysisMenu = uimenu(app.Menu);
            app.ExportAnalysisMenu.MenuSelectedFcn = createCallbackFcn(app, @ExportAnalysisMenuSelected, true);
            app.ExportAnalysisMenu.Text = 'Export Analysis';

            % Create ImageDisplayPanel
            app.ImageDisplayPanel = uipanel(app.MyoProfilerUIFigure);
            app.ImageDisplayPanel.Title = 'Image Display';
            app.ImageDisplayPanel.Position = [10 15 432 388];

            % Create ImageAxes
            app.ImageAxes = uiaxes(app.ImageDisplayPanel);
            app.ImageAxes.XTick = [];
            app.ImageAxes.YTick = [];
            app.ImageAxes.Box = 'on';
            app.ImageAxes.Visible = 'off';
            app.ImageAxes.Position = [17 10 399 314];

            % Create TabGroup
            app.TabGroup = uitabgroup(app.ImageDisplayPanel);
            app.TabGroup.Position = [4 6 424 354];

            % Create ImageTab
            app.ImageTab = uitab(app.TabGroup);
            app.ImageTab.Title = 'Image';

            % Create ProfilerPanel
            app.ProfilerPanel = uipanel(app.MyoProfilerUIFigure);
            app.ProfilerPanel.Title = 'Profiler Panel';
            app.ProfilerPanel.Position = [451 16 432 579];

            % Create ProfileIntensity
            app.ProfileIntensity = uiaxes(app.ProfilerPanel);
            title(app.ProfileIntensity, 'Extracted Profile')
            xlabel(app.ProfileIntensity, 'Horizontal Location (px)')
            ylabel(app.ProfileIntensity, 'Vertical Location (px)')
            zlabel(app.ProfileIntensity, {'Normalized Intensity'; 'Over Full Range'})
            app.ProfileIntensity.Box = 'on';
            app.ProfileIntensity.Position = [8 366 415 182];

            % Create ProfileIntensityXCoord
            app.ProfileIntensityXCoord = uiaxes(app.ProfilerPanel);
            title(app.ProfileIntensityXCoord, 'Profile Intensity')
            xlabel(app.ProfileIntensityXCoord, 'Horizontal Location (px)')
            ylabel(app.ProfileIntensityXCoord, {'Normalized Intensity'; 'Over Full Range'})
            zlabel(app.ProfileIntensityXCoord, 'Z')
            app.ProfileIntensityXCoord.Box = 'on';
            app.ProfileIntensityXCoord.Position = [10 190 415 160];

            % Create ProfileIntensityYCoord
            app.ProfileIntensityYCoord = uiaxes(app.ProfilerPanel);
            title(app.ProfileIntensityYCoord, 'Profile Intensity')
            xlabel(app.ProfileIntensityYCoord, 'Vertical Location (px)')
            ylabel(app.ProfileIntensityYCoord, {'Normalized Intensity'; 'Over Full Range'})
            zlabel(app.ProfileIntensityYCoord, 'Z')
            app.ProfileIntensityYCoord.Box = 'on';
            app.ProfileIntensityYCoord.Position = [10 13 415 160];

            % Create ControlsPanel
            app.ControlsPanel = uipanel(app.MyoProfilerUIFigure);
            app.ControlsPanel.Title = 'Controls';
            app.ControlsPanel.Position = [11 409 432 187];

            % Create ChannelColormapPanel
            app.ChannelColormapPanel = uipanel(app.ControlsPanel);
            app.ChannelColormapPanel.Title = 'Channel Colormap';
            app.ChannelColormapPanel.Position = [6 105 164 57];

            % Create ChannelColorDropDown
            app.ChannelColorDropDown = uidropdown(app.ChannelColormapPanel);
            app.ChannelColorDropDown.Items = {'Em. Wavelength', 'Parula'};
            app.ChannelColorDropDown.ValueChangedFcn = createCallbackFcn(app, @ChannelColorDropDownValueChanged, true);
            app.ChannelColorDropDown.Position = [5 6 121 22];
            app.ChannelColorDropDown.Value = 'Em. Wavelength';

            % Create ROIPanel
            app.ROIPanel = uipanel(app.ControlsPanel);
            app.ROIPanel.Title = 'ROI';
            app.ROIPanel.Position = [174 105 253 57];

            % Create SelectPointsButton
            app.SelectPointsButton = uibutton(app.ROIPanel, 'push');
            app.SelectPointsButton.ButtonPushedFcn = createCallbackFcn(app, @SelectPointsButtonPushed, true);
            app.SelectPointsButton.Position = [5 7 85 23];
            app.SelectPointsButton.Text = 'Select Points';

            % Create ROIWidthpxEditFieldLabel
            app.ROIWidthpxEditFieldLabel = uilabel(app.ROIPanel);
            app.ROIWidthpxEditFieldLabel.HorizontalAlignment = 'center';
            app.ROIWidthpxEditFieldLabel.Position = [95 6 114 22];
            app.ROIWidthpxEditFieldLabel.Text = 'ROI Width (px)';

            % Create ROIWidthpxEditField
            app.ROIWidthpxEditField = uieditfield(app.ROIPanel, 'numeric');
            app.ROIWidthpxEditField.ValueChangedFcn = createCallbackFcn(app, @ROIWidthpxEditFieldValueChanged, true);
            app.ROIWidthpxEditField.Position = [200 6 32 22];
            app.ROIWidthpxEditField.Value = 1;

            % Create ProfilePanel
            app.ProfilePanel = uipanel(app.ControlsPanel);
            app.ProfilePanel.Title = 'Profile';
            app.ProfilePanel.Position = [6 6 421 96];

            % Create ABandProminenceEditFieldLabel
            app.ABandProminenceEditFieldLabel = uilabel(app.ProfilePanel);
            app.ABandProminenceEditFieldLabel.HorizontalAlignment = 'right';
            app.ABandProminenceEditFieldLabel.Position = [205 40 111 22];
            app.ABandProminenceEditFieldLabel.Text = 'A Band Prominence';

            % Create ABandProminenceEditField
            app.ABandProminenceEditField = uieditfield(app.ProfilePanel, 'numeric');
            app.ABandProminenceEditField.Limits = [0 1];
            app.ABandProminenceEditField.ValueChangedFcn = createCallbackFcn(app, @ABandProminenceEditFieldValueChanged, true);
            app.ABandProminenceEditField.Position = [368 40 39 22];
            app.ABandProminenceEditField.Value = 0.05;

            % Create CalibrationumpxEditFieldLabel
            app.CalibrationumpxEditFieldLabel = uilabel(app.ProfilePanel);
            app.CalibrationumpxEditFieldLabel.HorizontalAlignment = 'right';
            app.CalibrationumpxEditFieldLabel.Position = [205 9 106 22];
            app.CalibrationumpxEditFieldLabel.Text = 'Calibration (um/px)';

            % Create CalibrationumpxEditField
            app.CalibrationumpxEditField = uieditfield(app.ProfilePanel, 'numeric');
            app.CalibrationumpxEditField.Limits = [0 Inf];
            app.CalibrationumpxEditField.ValueChangedFcn = createCallbackFcn(app, @CalibrationumpxEditFieldValueChanged, true);
            app.CalibrationumpxEditField.Position = [337 9 70 22];
            app.CalibrationumpxEditField.Value = 1;

            % Create PeakDistancepxEditFieldLabel
            app.PeakDistancepxEditFieldLabel = uilabel(app.ProfilePanel);
            app.PeakDistancepxEditFieldLabel.HorizontalAlignment = 'center';
            app.PeakDistancepxEditFieldLabel.WordWrap = 'on';
            app.PeakDistancepxEditFieldLabel.Position = [15 4 104 32];
            app.PeakDistancepxEditFieldLabel.Text = 'Peak Distance (px)';

            % Create PeakDistancepxEditField
            app.PeakDistancepxEditField = uieditfield(app.ProfilePanel, 'numeric');
            app.PeakDistancepxEditField.Limits = [0 Inf];
            app.PeakDistancepxEditField.ValueChangedFcn = createCallbackFcn(app, @PeakDistancepxEditFieldValueChanged, true);
            app.PeakDistancepxEditField.Position = [134 9 39 22];
            app.PeakDistancepxEditField.Value = 25;

            % Create ZLineProminenceEditFieldLabel
            app.ZLineProminenceEditFieldLabel = uilabel(app.ProfilePanel);
            app.ZLineProminenceEditFieldLabel.HorizontalAlignment = 'right';
            app.ZLineProminenceEditFieldLabel.Position = [11 40 109 22];
            app.ZLineProminenceEditFieldLabel.Text = 'Z Line Prominence ';

            % Create ZLineProminenceEditField
            app.ZLineProminenceEditField = uieditfield(app.ProfilePanel, 'numeric');
            app.ZLineProminenceEditField.Limits = [0 1];
            app.ZLineProminenceEditField.ValueChangedFcn = createCallbackFcn(app, @ZLineProminenceEditFieldValueChanged, true);
            app.ZLineProminenceEditField.Position = [134 40 39 22];
            app.ZLineProminenceEditField.Value = 0.5;

            % Create SarcomerePanel
            app.SarcomerePanel = uipanel(app.MyoProfilerUIFigure);
            app.SarcomerePanel.Title = 'Sarcomere Panel';
            app.SarcomerePanel.Position = [894 16 432 580];

            % Create SarcomereMean
            app.SarcomereMean = uiaxes(app.SarcomerePanel);
            title(app.SarcomereMean, 'Average Sarcomere Intensity')
            xlabel(app.SarcomereMean, 'Location (um)')
            ylabel(app.SarcomereMean, {'Normalized Intensity'; 'Over Sarcomere'})
            app.SarcomereMean.Box = 'on';
            app.SarcomereMean.Position = [9 390 415 160];

            % Create Sarcomeres
            app.Sarcomeres = uiaxes(app.SarcomerePanel);
            title(app.Sarcomeres, 'Sarcomere Intensities')
            xlabel(app.Sarcomeres, 'Location (um)')
            ylabel(app.Sarcomeres, {'Normalized Intensity'; 'Over Sarcomere'})
            app.Sarcomeres.Box = 'on';
            app.Sarcomeres.Position = [8 209 415 160];

            % Create SummaryTable
            app.SummaryTable = uitable(app.SarcomerePanel);
            app.SummaryTable.ColumnName = {'Channel'; 'Profile Number'; 'Color'; 'SL (um)'; 'FWHM (um)'};
            app.SummaryTable.RowName = {};
            app.SummaryTable.Position = [8 12 417 179];

            % Create ZLineStainCheckBox
            app.ZLineStainCheckBox = uicheckbox(app.MyoProfilerUIFigure);
            app.ZLineStainCheckBox.Visible = 'off';
            app.ZLineStainCheckBox.Enable = 'off';
            app.ZLineStainCheckBox.Text = 'Z-Line Stain';
            app.ZLineStainCheckBox.Position = [231 -65 87 22];

            % Show the figure after all components are created
            app.MyoProfilerUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = myofibril_profiler_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.MyoProfilerUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.MyoProfilerUIFigure)
        end
    end
end