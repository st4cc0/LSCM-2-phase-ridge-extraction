% (CC 2022 by Lukas Hauer)
close all
clearvars
clc

path_in     = '01_Files to Process';
path_out    = '02_Results';
files_id    = dir([path_in,'/*.tif']);




for i=1:numel(files_id)
    
    this_path_in = [path_in,'/',files_id(i).name];
    this_path_out = [path_out,'/',files_id(i).name(1:end-4)];
    mkdir(this_path_out)
    mkdir([this_path_out,'/frames'])
    
    process_lscm_data(this_path_in,this_path_out)

end




function process_lscm_data(this_path_in,this_path_out)
        %% some constant extractions
    tiff_info   = imfinfo(this_path_in);
    frame_num   = numel(tiff_info)/2;
    pxl_size    = 0.2443794;    %um/pxl
    n_sioil     = 1.4;
    d_time      = 0.63;         %s/frame
    
    lub_col = [215 0 0]./255;
    net_col = [215 215 0]./255;
        %% allocate variables
    ridge_area_lub  = nan(frame_num,1);
    ridge_area_net  = nan(frame_num,1);
    ridge_pos       = nan(frame_num,2);
    ridge_hmax_lub  = nan(frame_num,1);
    ridge_hmax_net  = nan(frame_num,1);
    rigdge_len      = nan(frame_num,1);
    time            = nan(frame_num,1);
        %% process images
    for i=1:frame_num
        try
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
            
            ridge_profile_lub   = profile_lub(1:ridge_tip_x);
            ridge_profile_net   = profile_net(1:ridge_tip_x);
                                
            [base_line_lub,dep_point_lub] = get_baseLine(profile_lub,ridge_tip_x);
            [base_line_net,dep_point_net] = get_baseLine(profile_net,ridge_tip_x);
        
            int_start = mean(dep_point_lub,dep_point_net);
            
            rigdge_len(i) = abs(int_start-ridge_tip_x)*pxl_size;
        
            av_base = 0.5*(base_line_lub(1:ridge_tip_x)+base_line_net(1:ridge_tip_x));
        
            ridge_profile_lub   = ridge_profile_lub - av_base;
            ridge_profile_net   = ridge_profile_net - av_base;
    
%             [~,outs_idx_lub] = rmoutliers(ridge_profile_lub,'movmedian',10);
%             [~,outs_idx_lub] = rmoutliers(ridge_profile_net,'movmedian',10);
            
%             ridge_profile_lub(outs_idx_lub) = nan;
%             ridge_profile_net(outs_idx_net) = nan;
            
%             ridge_profile_lub = close_nan(ridge_profile_lub);
%             ridge_profile_net = close_nan(ridge_profile_net);
            
            f_pass = 1/10;
            ridge_profile_lub_lp = lowpass(ridge_profile_lub,f_pass,...
                                    1/pxl_size,'Steepness',.9999);
            ridge_profile_net_lp = lowpass(ridge_profile_net,f_pass,...
                                    1/pxl_size,'Steepness',.9999);
        
            ridge_area_lub(i)   = sum(ridge_profile_lub_lp(int_start:end))*pxl_size;
            ridge_area_net(i)   = sum(ridge_profile_net_lp(int_start:end))*pxl_size; 
            
            ridge_hmax_lub(i)   = max(ridge_profile_lub_lp)*pxl_size;
            ridge_hmax_net(i)   = max(ridge_profile_net_lp)*pxl_size;
                %% plot
                if 1 == 1
                    lub     = imadjust(lub);
                    lub     = image_crop(lub,0.1);
    
                    net     = imadjust(net);
                    net     = image_crop(net,0.1);
    
                    [nr,nc] = size(net);      
                    int_map = 0.8;
            
                    red_map = [linspace(0,int_map,256);...
                                        zeros(1,256);...
                                        zeros(1,256)]';
                    figure(1)
                    clf
                    subplot('position',[0,0.5,1,0.5])
                    imshow(lub,red_map)
                    hold on
                    plot(profile_lub/n_sioil,'r.')
                    plot(profile_net/n_sioil,'y.')
                    plot_this_point(ridge_tip_x,ridge_tip_y/n_sioil,'white',100)
                    set(gca,'ydir','normal','xdir','reverse')
                    
                    yellow_map = [linspace(0,int_map,256);...
                                    linspace(0,int_map,256);...
                                    zeros(1,256)]';
                        
                    subplot('position',[0,0,1,0.5])
                    imshow(net,yellow_map)
                    hold on
                    plot(profile_lub/n_sioil,'r.')
                    plot(profile_net/n_sioil,'y.')
    
                    x_pos = 100;
                    y_pos = 50;
    
                    set(gcf,'position',[x_pos,y_pos,nc,2*nr])
                    
                    str = sprintf('%1.1f s',time(i));
                        dim = [0.68 0.89 0.1 0.1];
                        annotation('textbox',dim,'String',str,'Color',[1 1 1],...
                            'FontName','CMU Serif','FontWeight','Bold',...
                            'FontSize',18,'LineStyle','none');
                    set(gca,'ydir','normal','xdir','reverse')
                    set(gcf, 'InvertHardcopy', 'off')
                    saveas(gca,[this_path_out,sprintf('/frames/frame_%06i.jpg',i)])
%               pause(.1)
                    
                    figure(2)
                    clf
                    subplot(2,1,1)
                    hold on
                    x_len = (0:length(ridge_profile_lub)-1)*pxl_size;
                    
                    plot_this_point(x_len,flip(ridge_profile_lub)*pxl_size,lub_col,40)
                    plot_this_point(x_len,flip(ridge_profile_net)*pxl_size,net_col,40)
                    plot(x_len,flip(ridge_profile_lub_lp)*pxl_size,'Color',[1 .3 .3])
                    plot(x_len,flip(ridge_profile_net_lp)*pxl_size,'Color',[1 1 .3]) 
                    xline(int_start*pxl_size,'k','$\mathrm{ridge~domain}$','LabelHorizontalAlignment','left',...
                        'Interpreter','Latex')
%                     yline(0,'w:')
                    yline(0,'k:')
                    
%                     set(gca,'Color','k','XColor',[1 1 1],'YColor',[1 1 1])
%                     set(gcf,'Color','k')

                    set(gca,'Color','w')
                    set(gcf,'Color','w')
                    
                    set(gca,'TickLabelInterpreter','latex');
                    ylabel('$h~\mathrm{(\mu m)}$','Interpreter','Latex')
                    xlabel('$x~\mathrm{(\mu m)}$','Interpreter','Latex')
                    axis equal
                    xlim([0,max(x_len)])
                    ylim([-5 35])
                    subplot(2,1,2)
                    hold on
                    x_len = ((int_start:length(ridge_profile_lub))-int_start)*pxl_size;
                    plot_this_point(x_len,flip(ridge_profile_lub(int_start:end))*pxl_size,lub_col,40)
                    plot_this_point(x_len,flip(ridge_profile_net(int_start:end))*pxl_size,net_col,40)
                    plot(x_len,flip(ridge_profile_lub_lp(int_start:end))*pxl_size,'Color',[1 .6 .6],'LineWidth',1.5)
                    plot(x_len,flip(ridge_profile_net_lp(int_start:end))*pxl_size,'Color',[1 1 .6],'LineWidth',1.5) 
                    yline(0,'k:')
                    
%                     set(gca,'Color','k','XColor',[1 1 1],'YColor',[1 1 1])
%                     set(gcf,'Color','k')
                    set(gca,'Color','w')
                    set(gcf,'Color','w')

                    set(gca,'TickLabelInterpreter','latex');
                    ylabel('$h~\mathrm{(\mu m)}$','Interpreter','Latex')
                    xlabel('$x~\mathrm{(\mu m)}$','Interpreter','Latex')
                    ylim([-2 35])
                    xlim([0,max(x_len)])
                    set(gcf, 'InvertHardcopy', 'off')
                    saveas(gca,[this_path_out,sprintf('/frames/ridge_%06i.jpg',i)])
                    
                end
        catch
            fprintf('error in frame %i \n',i)
        end
    end
    
    close all
    
    subplot(4,1,1)
    plot(time,ridge_area_lub,'.r')
    hold on
    plot(time,ridge_area_net,'.y')
    xlim([0 max(time)])
    xticklabels([])
    ylabel('Area (µm^2)')
    
    subplot(4,1,2)
    plot(time,ridge_hmax_lub,'.r')
    hold on
    plot(time,ridge_hmax_net,'.y')
    xlim([0 max(time)])
    xticklabels([])
    ylabel('h_{max} (µm)')
    
    subplot(4,1,3)
    plot(time,rigdge_len,'x-')
    xlim([0 max(time)])
    xticklabels([])
    ylabel('ridge length (µm)')
    
    subplot(4,1,4)
    ridge_rel_speed = (ridge_pos(2:end,1)-ridge_pos(1:end-1,1))/d_time;
    speed_time = 0.5*(time(2:end)+time(1:end-1));
    semilogy(speed_time(1:end-1),abs(ridge_rel_speed(1:end-1)),'x-')
    xlim([0 max(time)])
    xlabel('time (s)')
    ylabel('l_{ridge} (µm/s)')
    
    saveas(gca,[this_path_out,'/some_specs.jpg'])
    close all
    
    csvwrite([this_path_out,'/ridge_area_lub.csv'],ridge_area_lub)
    csvwrite([this_path_out,'/ridge_area_net.csv'],ridge_area_net)
    csvwrite([this_path_out,'/ridge_hmax_lub.csv'],ridge_hmax_lub)
    csvwrite([this_path_out,'/ridge_hmax_net.csv'],ridge_hmax_net)
    csvwrite([this_path_out,'/rigdge_len.csv'],rigdge_len)
    csvwrite([this_path_out,'/ridge_pos.csv'],ridge_pos)
    
    n_lub           = length(ridge_hmax_lub);    
    lub_height_mean = mean(ridge_hmax_lub,'omitnan');
    lub_height_std  = std(ridge_hmax_lub,'omitnan');
    
    n_net           = length(ridge_hmax_net);
    net_height_mean = mean(ridge_hmax_net,'omitnan');
    net_height_std  = std(ridge_hmax_net,'omitnan');
    
    idx_start   = strfind(this_path_out,'crop.lif -')+length('crop.lif -')+1;
    idx_end     = strfind(this_path_out,'um')-1;
    speed       = str2double(this_path_out(idx_start:idx_end));
    
    T = table(speed,n_lub,lub_height_mean,lub_height_std,n_net,...
                net_height_mean,net_height_std);
    writetable(T,[this_path_out,'/average_ridge_height.csv'])
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

function [img_bin,profile] = get_interface_profile(img,w0)

    img_filt        = imadjust(img);
    
    img_renorm = wiener2(img_filt,[15 15]);
    img_renorm = imadjust(imgaussfilt(img_renorm,w0*[1 3]));
    
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

function plot_this_point(x,y,color,size)
    p1 = plot(x,y,'.','color',color);
    scatter(x,y,size,'filled', ...
             'MarkerFaceAlpha',1/2,'MarkerFaceColor',p1.Color);
end
