%% This script extracts temporal averages of the ridge profile
% (CC 2022 by Lukas Hauer)

close all
clearvars
clc

path_in     = '01_Files to Process';
path_out    = '02_Results';
files_id    = dir([path_in,'/*.tif']);




for i=1:numel(files_id)
    
    try
        this_path_in = [path_in,'/',files_id(i).name];
        this_path_out = [path_out,'/',files_id(i).name(1:end-4)];
        mkdir(this_path_out)
    
        process_lscm_data(this_path_in,this_path_out)
    catch
        disp(['error in ',files_id(i).name])
    end

end




function process_lscm_data(this_path_in,this_path_out)
        %% some constant extractions
    tiff_info   = imfinfo(this_path_in);
    frame_num   = numel(tiff_info)/2;
    pxl_size    = 0.2443794;    %um/pxl
    n_sioil     = pxl_size*1.4;
    d_time      = 0.63;         %s/frame
    
    lub_col = [215 0 0]./255;
    net_col = [215 215 0]./255;
        %% allocate variables
    ridge_pos       = nan(frame_num,2);
    rigdge_len      = nan(frame_num,1);
    time            = nan(frame_num,1);
    
    ridge_air_lub       = nan;
    ridge_drop_lub      = nan;
    ridge_air_net       = nan;
    ridge_drop_net      = nan;
    
        %% process images
    for i=1:frame_num
%         try
            net_idx = i*2-1;
            lub_idx = i*2;
            
            time(i) = (i-1)*d_time;
    
            lub = imread(this_path_in, lub_idx);
            net = imread(this_path_in, net_idx);
    
            [~,profile_lub] = get_interface_profile(lub,1);
            [~,profile_net] = get_interface_profile(net,2.5);
            
            profile_lub = profile_lub*n_sioil;
            profile_net = profile_net*n_sioil;
            
            [ridge_tip_y,ridge_tip_x] = max(profile_lub);
            ridge_pos(i,:) = [ridge_tip_x,ridge_tip_y];
            
                                
            [base_line_lub,dep_point_lub] = get_baseLine(profile_lub,ridge_tip_x);
            [base_line_net,dep_point_net] = get_baseLine(profile_net,ridge_tip_x);
        
            int_start = mean(dep_point_lub,dep_point_net);
            
            rigdge_len(i) = abs(int_start-ridge_tip_x)*pxl_size;
                
            ridge_profile_lub   = profile_lub - base_line_lub;
            ridge_profile_net   = profile_net - base_line_net;           
            
                %adding lub-air side
            this_ridge_air_lub  = flip(ridge_profile_lub(1:ridge_tip_x));
            ridge_air_lub       = dim_add(this_ridge_air_lub,...
                                            ridge_air_lub);
                %adding lub-drop side
            this_ridge_drop_lub  = ridge_profile_lub(ridge_tip_x:end);
            ridge_drop_lub       = dim_add(this_ridge_drop_lub,...
                                            ridge_drop_lub);           
                %adding net-air side                            
            this_ridge_air_net  = flip(ridge_profile_net(1:ridge_tip_x));
            ridge_air_net       = dim_add(this_ridge_air_net,...
                                            ridge_air_net);
                %adding net-drop side
            this_ridge_drop_net  = ridge_profile_net(ridge_tip_x:end);
            ridge_drop_net       = dim_add(this_ridge_drop_net,...
                                            ridge_drop_net);                   
    end
    %% Post Processing
    %Statistics
        %Lubricant profile 
    [ridge_air_lub_mean,...
        ridge_air_lub_std,...
        ridge_air_lub_n]        = statistics(ridge_air_lub);       
    x_air_lub                   = pxl_size*(0:length(ridge_air_lub_mean)-1)';
    
    [ridge_drop_lub_mean,...
        ridge_drop_lub_std,...
        ridge_drop_lub_n]       = statistics(ridge_drop_lub);
    x_drop_lub                  = -pxl_size*(1:length(ridge_drop_lub_mean))';       
        %Network profile
    [ridge_air_net_mean,...
        ridge_air_net_std,...
        ridge_air_net_n]        = statistics(ridge_air_net);
    x_air_net                   = pxl_size*(0:length(ridge_air_net_mean)-1)';
    
    [ridge_drop_net_mean,...
        ridge_drop_net_std,...
        ridge_drop_net_n]        = statistics(ridge_drop_net);
    x_drop_net                   = -pxl_size*(1:length(ridge_drop_net_mean))';
    
        %save data
    T_lub = table([flip(x_drop_lub);x_air_lub],...
                [flip(ridge_drop_lub_mean);ridge_air_lub_mean],...
                [flip(ridge_drop_lub_std);ridge_air_lub_std],...
                [flip(ridge_drop_lub_n);ridge_air_lub_n],...
                'VariableNames',{'x','height','heightStd','n'});
    writetable(T_lub,[this_path_out,'/Lub_ridge_height.csv']);
    
    T_net = table([flip(x_drop_net);x_air_net],...
                [flip(ridge_drop_net_mean);ridge_air_net_mean],...
                [flip(ridge_drop_net_std);ridge_air_net_std],...
                [flip(ridge_drop_net_n);ridge_air_net_n],...
                'VariableNames',{'x','height','heightStd','n'});
    writetable(T_net,[this_path_out,'/Net_ridge_height.csv']);
    
        %statistical cutoff
    ridge_air_lub_mean  = stat_cut_off(ridge_air_lub_mean,ridge_air_lub_n);
    ridge_drop_lub_mean = stat_cut_off(ridge_drop_lub_mean,ridge_drop_lub_n);
    ridge_air_net_mean  = stat_cut_off(ridge_air_net_mean,ridge_air_net_n);
    ridge_drop_net_mean = stat_cut_off(ridge_drop_net_mean,ridge_drop_net_n);
    
        %Max Ridge Points
    lub_height_mean = ridge_air_lub_mean(1);
    lub_height_std  = ridge_air_lub_std(1);
    n_lub           = ridge_drop_lub_n(1);
    net_height_mean = ridge_air_net_mean(1);
    net_height_std  = ridge_air_net_std(1);
    n_net           = ridge_air_net_n(1);
    
    idx_start   = strfind(this_path_out,'crop.lif -')+length('crop.lif -')+1;
    idx_end     = strfind(this_path_out,'um')-1;
    speed       = str2double(this_path_out(idx_start:idx_end));
    
    T = table(speed,n_lub,lub_height_mean,lub_height_std,n_net,...
                net_height_mean,net_height_std);
    writetable(T,[this_path_out,'/average_ridge_height.csv'])
    
    
    figure(1)
    yline(0,':')
    xline(0,':')
    hold on
    plot_this_error(x_air_lub,ridge_air_lub_mean,ridge_air_lub_std./std(ridge_air_lub_n),lub_col)
    plot_this_error(x_air_net,ridge_air_net_mean,ridge_air_net_std./sqrt(ridge_air_net_n),net_col)
    
    plot_this_error(x_drop_lub,ridge_drop_lub_mean,ridge_drop_lub_std./std(ridge_drop_lub_n),lub_col)
    plot_this_error(x_drop_net,ridge_drop_net_mean,ridge_drop_net_std./sqrt(ridge_drop_net_n),net_col)
    
    xlabel('$x~\mathrm{\mu m}$','Interpreter','Latex')
    ylabel('$y~\mathrm{\mu m}$','Interpreter','Latex')
    xlim([-120,130])
    ylim([-10,40])
    style_plot()
    saveas(gca,[this_path_out,'/entire_profile.png'])
    
    figure(2)
    yline(0,':')
    hold on
    plot_this_error(x_air_lub,ridge_air_lub_mean,ridge_air_lub_std./std(ridge_air_lub_n),lub_col)   
    plot_this_error(x_air_net,ridge_air_net_mean,ridge_air_net_std./sqrt(ridge_air_net_n),net_col)
    xlabel('$x~\mathrm{\mu m}$','Interpreter','Latex')
    ylabel('$y~\mathrm{\mu m}$','Interpreter','Latex')
    xlim([0,130])
    ylim([-1,35])
    style_plot()
    saveas(gca,[this_path_out,'/air_profile.png'])
    
    close all
    end


%% functions
function [base_line,dep_point] = get_baseLine(profile,ridge_x)

    profile_left = profile(1:ridge_x);
    
    profile_mode    = mode(profile_left);
    fit_band        = 0.5*std(profile_left,'omitnan');

    fit_idx = find(profile_mode-fit_band<profile_left &...
                    profile_left<profile_mode+fit_band);
                
    fit_profile = profile_left(fit_idx);
    
    base_fit    = polyfit(fit_idx,fit_profile,1);
    base_line   = base_fit(1)*(1:length(profile))+base_fit(2);
    base_line   = base_line';
    dep_point   = max(fit_idx);
end

function [A_mean,A_std,A_n] = statistics(A)
    A_mean  = mean(A,2,'omitnan');
    A_std   = std(A,0,2,'omitnan');
    A_n     = sum(1-isnan(A),2); 
end

function array = stat_cut_off(array,n_array)
        this_stat_cut_off       = mean(n_array)-std(n_array);
        this_stat_cut_off_idx   = this_stat_cut_off>n_array;
        array(this_stat_cut_off_idx) = nan;
end

function add_up_array = dim_add(new_array,old_array)
    
    [l_a,~]     = size(new_array);
    [l_b,row]   = size(old_array);
        
    if l_a < l_b
        add_up_array = nan(l_b,row+1);
    else
        add_up_array = nan(l_a,row+1);
    end
    
    add_up_array(1:l_b,1:row) = old_array;
    add_up_array(1:l_a,row+1)     = new_array;
    1;
end

function [img_bin,profile] = get_interface_profile(img,w0)

    img_filt        = imadjust(img);
    
    img_renorm      = wiener2(img_filt,[15 15]);
    img_renorm      = imadjust(imgaussfilt(img_renorm,w0*[1 3]));
    img_renorm      = image_crop(img_renorm,0.1);
    img_bin         = imbinarize(img_renorm);

    profile         = extract_profile(img_bin);
end

function profile_this = extract_profile(img)
    [~,nc]=size(img);
    profile_this = nan(nc,1);
    for i=1:nc
        try 
            profile_this(i) = find(img(:,i),1,'last');
        catch
            profile_this(i) = nan;
        end
    end
end

function img = image_crop(img,pct)
    
    [nr,~] = size(img);
    crop_up = round(pct*nr);
    crop_down = round((1-pct)*nr);
    
    img = img(crop_up:crop_down,:);

end

function plot_this_error(x,y,y_err,color)
    errorbar(x,y,y_err,'.','color',color);
%     p1 = scatter(x,y,10,color,'filled');
    hold on
    scatter(x,y,20,color,'filled', ...
             'MarkerFaceAlpha',1/2);
end

function style_plot()
    set(gca,'Color','w','XColor',[0 0 0],'YColor',[0 0 0])
    set(gcf,'Color','w')
    set(gca,'TickLabelInterpreter','latex');
    set(gcf, 'InvertHardcopy', 'off')
    set(gca,'FontSize',18)
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
end
