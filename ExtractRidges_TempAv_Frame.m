%% T
%This script generates temporal averaged LSCM images

close all
clearvars
clc

path_in     = '01_Files to Process';
path_out    = '02_Results';
files_id    = dir([path_in,'/*.tif']);




for i=1:numel(files_id)
    
    try
        this_path_in = [path_in,'/',files_id(i).name];
        this_path_out = [path_out,'/',files_id(i).name(1:end-4),'_tempAV_LSCM'];
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
    d_time      = 0.63;         %s/frame
    
        %% allocate variables
    time            = nan(frame_num,1);
    
        %% process images
    for i=1:frame_num
%         try
            
            net_idx = i*2-1;
            lub_idx = i*2;
            
            time(i) = (i-1)*d_time;
    
            lub = imread(this_path_in, lub_idx);
            net = imread(this_path_in, net_idx);
            
            net = wiener2(net,[5 5]);
            net = imadjust(imgaussfilt(net,[1 2]));      
            
            [~,profile_lub] = get_interface_profile(lub,1);
            
            [ridge_tip_y,ridge_tip_x] = max(profile_lub);
            ridge_pos_lub = [ridge_tip_x,ridge_tip_y];
            
            if i==1
                lub_merged_array    = zeros(size(lub));
                net_merged_array    = zeros(size(lub));
                cent_old            = ridge_pos_lub;
            end
            
            [lub_merged_array,~] = fix_those_arrays(lub_merged_array,...
                                                lub,cent_old,ridge_pos_lub);
            [net_merged_array,cent_old] = fix_those_arrays(net_merged_array,...
                                                net,cent_old,ridge_pos_lub);                             
    end
    %% Post Processing
    lub_merged_array = lub_merged_array./frame_num;
    lub_merged_array = uint8(lub_merged_array);
    lub_merged_array = imadjust(lub_merged_array);
        
    net_merged_array = net_merged_array./frame_num;
    net_merged_array = uint8(net_merged_array);
    net_merged_array = imadjust(net_merged_array);
    
    
    
    plot_case = 1;
        %case 1: plot LSCM with black background
        %case 2: plot LSCM with white background
        %case 3: plot sigmoidal fit
    
    switch plot_case
        case 1
        int_map = 0.8;
        red_map = [linspace(0,int_map,256);...
                                        zeros(1,256);...
                                        zeros(1,256)]';
        subplot('position',[0,0.5,1,0.5])
        imshow(lub_merged_array,red_map)
        set(gca,'ydir','normal')
    
        yellow_map = [linspace(0,int_map,256);...
                                    linspace(0,int_map,256);...
                                    zeros(1,256)]';
                                
        subplot('position',[0,0,1,0.5],'ydir','reverse')
        imshow(net_merged_array,yellow_map)
        set(gca,'ydir','normal')
    
        imwrite(flip(flip(net_merged_array,1),2),[this_path_out,'/net_channel.png'])
        imwrite(flip(flip(lub_merged_array,1),2),[this_path_out,'/lub_channel.png'])
        
        case 2
        int_map = 0.8;
        
        red_map = [linspace(1,int_map,256);...
                   linspace(1,0,256);...
                   linspace(1,0,256)]';
        
        yellow_map = [linspace(1,int_map,256);...
                      linspace(1,int_map,256);...
                      linspace(1,0,256)]';
                  
%         r_map = meshgrid(red_map(:,1),yellow_map(:,1));
%         g_map = meshgrid(red_map(:,2),yellow_map(:,2));
%         b_map = meshgrid(red_map(:,3),yellow_map(:,3));
                                
        
        subplot('position',[0,0.5,1,0.5])
        imshow(lub_merged_array,red_map)
        set(gca,'ydir','normal')          
                  
        subplot('position',[0,0,1,0.5],'ydir','reverse')
        imshow(net_merged_array,yellow_map)
        set(gca,'ydir','normal')
    
        imwrite(net_merged_array,[this_path_out,'/net_channel.png'])
        imwrite(lub_merged_array,[this_path_out,'/lub_channel.png'])

        case 3
            
        int_map = 0.8;    
        
        [lub_sig_filt,~] = sigmoidal_filter(lub_merged_array(1:250,:));
        [net_sig_filt,~] = sigmoidal_filter(net_merged_array(1:250,:));
        
        lub_sig_filt = uint8(rescale(lub_sig_filt,1,256));
        net_sig_filt = uint8(rescale(net_sig_filt,1,256));        
        
        red_map = [linspace(1,int_map,256);...
                   linspace(1,0,256);...
                   linspace(1,0,256)]';
        
        yellow_map = [linspace(1,int_map,256);...
                      linspace(1,int_map,256);...
                      linspace(1,0,256)]';                  
                  
        subplot('position',[0,0.5,1,0.5])
        imshow(lub_sig_filt,red_map)
        set(gca,'ydir','normal')
        set(gca,'xdir','reverese')
                                
        subplot('position',[0,0,1,0.5],'ydir','reverse')
        imshow(net_sig_filt,yellow_map)
        set(gca,'ydir','normal')
    
        imwrite(net_merged_array,[this_path_out,'/net_channel.png'])
        imwrite(lub_merged_array,[this_path_out,'/lub_channel.png'])
    end
    1;
    end


%% functions

function [this_merged_array,c_out] = fix_those_arrays(this_a,this_b,this_c_a,this_c_b)
    try
        weigths = ones(size(this_a));
        [this_weigth_1,this_weigth_2,~] = ...
            merger_arrays(weigths,weigths,this_c_a,this_c_b);
        [this_merged_array_1,this_merged_array_2,c_out] = ...
            merger_arrays(this_a,this_b,this_c_a,this_c_b);
    catch

    end
    weight_array        = this_weigth_1 + this_weigth_2;
    weight_idx          = weight_array>0;
    this_merged_array   = this_merged_array_1 + this_merged_array_2;
%     this_merged_array(weight_idx) = this_merged_array(weight_idx)./weight_array(weight_idx);

    function [merged_array_1,merged_array_2,c_out] = merger_arrays(a,b,c_a,c_b)

        disp_vec        = round(c_a) - round(c_b);
        [m_a,n_a]       = size(a);
        [m_b,n_b]       = size(b);
        
        merged_array_1  = zeros(max(m_a,m_b)+abs(disp_vec(2)),max(n_a,n_b)+abs(disp_vec(1)));
        merged_array_2  = zeros(max(m_a,m_b)+abs(disp_vec(2)),max(n_a,n_b)+abs(disp_vec(1)));
    
        if disp_vec(1) >= 0 && disp_vec(2) >= 0
            disp_vec = abs(disp_vec)+1;
            merged_array_1(1:m_a,...
                                    1:n_a) = a;
            merged_array_2(disp_vec(2):disp_vec(2)+m_b-1,...
                                    disp_vec(1):disp_vec(1)+n_b-1) = b;
            c_out = c_a;
        elseif disp_vec(1) <= 0 && disp_vec(2) >= 0
            disp_vec = abs(disp_vec)+1;
            merged_array_1(1:m_a,...
                                    disp_vec(1):n_a+disp_vec(1)-1) = a;
            merged_array_2(disp_vec(2):disp_vec(2)+m_b-1,...
                                    1:n_b) = b;
            c_out = [c_b(1),c_a(2)];
        elseif disp_vec(1) >= 0 && disp_vec(2) <= 0
            disp_vec = abs(disp_vec)+1;
            merged_array_1(disp_vec(2):disp_vec(2)+m_a-1,...
                                    1:n_a) = a;
            merged_array_2(1:m_b,...
                                    disp_vec(1):disp_vec(1)+n_b-1) = b;
            c_out = [c_a(1),c_b(2)];
        elseif disp_vec(1) <= 0 && disp_vec(2) <= 0
            disp_vec = abs(disp_vec)+1;
            merged_array_1(disp_vec(2):disp_vec(2)+m_a-1,...
                                    disp_vec(1):disp_vec(1)+n_a-1) = a;
            merged_array_2(1:m_b,1:n_b) = b;
            c_out = c_b;
        end
    close all
    end
end

function [sigmoidal_filt,if_profile] = sigmoidal_filter(img)

    [nr,nc] = size(img);

    sigmoidal_filt  = zeros(nr,nc);
    if_profile      = zeros(nc,1);
    
    startpoint = [1 0.5*nr];
    lower_bounds = [0.8 1];
    upper_bounds = [100 nr];
    
    for i=1:nc
        profile = img(:,i);
        x       = 1:nr; 
        profile = rescale(profile,0,2);
        fit_type = '(1+tanh(a*(b-x)))';
%         fit_type = '(1+tanh((a-x)))';
        sigm=fit(x',profile,fit_type,'StartPoint',startpoint,...
                                    'Lower',lower_bounds,'Upper',upper_bounds);
    
        sigmoidal = 1+tanh(sigm.a*(sigm.b-x));
%         sigmoidal = 1+tanh(sigm.a-x);
        sigmoidal_filt(:,i) = sigmoidal;
    
%         if_profile(i) = sigm.b;
%         startpoint = [sigm.a sigm.b];
    end
end

function delta_mode = calibrate_base_line(profile_a,profile_b)

    profile_delta       = profile_a - profile_b;
    [delta_dist,edges]  = histcounts(profile_delta,64);
    deltas              = 0.5*(edges(2:end)+edges(1:end-1));
    [~,idx_max]         = max(delta_dist);
    delta_mode          = deltas(idx_max);

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