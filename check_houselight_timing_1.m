clear;

% data directory
datadirs{1} = 'C:\Users\User\Desktop\coco_analysis\data\CA3\2020_12_03';


% average sampling interval (s)
dt = 1/30;

for d=1:length(datadirs)
    
    datadir = datadirs{d};
    
    mkdir([datadir '/qualitycheck']);
    
    % =========================================================================
    % read the mean signal from behavioral videos/mean files (used to visualize
    % the timing of the houselight events).
    
    % make video list
    viddir = [datadir '/BehavCam_0'];
    file_pattern = [viddir '/*.avi'];
    files = dir(file_pattern);
    Nvids = length(files);
    nums = [];
    Vnums = [];
    
    if Nvids>0
        nums = [];
        for i1=1:Nvids
            temp = files(i1).name;
            temp(end-3:end) = [];
            if not(isempty(temp)) && isempty(find(temp==' '))
                nums = [nums; i1 str2num(temp)];
            end
        end
        
        % sort the file numbers
        [u,v] = sort(nums(:,2),'ascend');
        Vnums = nums(v,:);
    end
    
    % make mean signal file list
    meandir = [datadir '/qualitycheck'];
    file_pattern = [meandir '/bh_mean*.mat'];
    files = dir(file_pattern);
    Nmeans = length(files);
    nums = [];
    Mnus = [];
    
    if Nmeans>0
        nums = [];
        for i1=1:Nmeans
            temp = files(i1).name;
            temp(1:length('bh_mean')) = [];
            temp(end-3:end) = [];
            if not(isempty(temp)) && isempty(find(temp==' '))
                nums = [nums; i1 str2num(temp)];
            end
        end
        
        % sort the file numbers
        [u,v] = sort(nums(:,2),'ascend');
        Mnums = nums(v,:);
    end
    
    % merge the signals from the .mat and (if necessary).avi files
    Nfiles = max(size(Vnums,1),size(Mnums,1));
    
    mean_bh = [];
    
    for f=1:Nfiles
        
        avifile = '';
        if Nvids>0
            numv = Vnums(f,2);
            avifile = [viddir '/' num2str(numv) '.avi'];
        end
        
        meanfile = '';
        if Nmeans>0
            numm = Mnums(f,2);
            meanfile = [meandir '/bh_mean' num2str(numm) '.mat'];
        end
        
        % load the mat file if it exists, the avi file if not, or stop
        % otherwise
        if exist(meanfile,'file')==2
            disp(['Loading mean file ' meanfile '...']);
            load(meanfile);
            mean_bh = [mean_bh bh_mean_signal];
        else
            if exist(avifile,'file')==2
                msvidObj = VideoReader(avifile);
                disp(['Loading video ' avifile '...']);
                video=msvidObj.read();
                bh_mean_signal = squeeze(mean(mean(mean(video,1),2),3));
                mean_bh = [mean_bh; bh_mean_signal];
            else
                disp('Cannot find mean signal file. Stopping here.')
            end
        end
    end
    
    
    % =========================================================================
    % extract the houselight events from the schedules file
    
    schedules = [datadir '/schedules.csv'];
    
    % do some schedule.csv file cleanup
    
    % remove all the '"'
    fid = fopen(schedules,'rt');
    X = fread(fid);
    fclose(fid);
    X = char(X.');
    % replace string S1 with string S2
    Y = strrep(X, '"', '');
    
    % remove the file's header which finishes at 'Arg5_Value'
    pos_start = strfind(Y,'Arg5_Value');
    if not(isempty(pos_start))
        pos_start = pos_start + length('Arg5_Value');
        Y(1:pos_start) = [];
    end
    
    fid2 = fopen(schedules,'wt');
    fwrite(fid2,Y);
    fclose (fid2);
   
    
    % declare arrays to speed up reading files
    % 1 = time, 2 = on(1)/off(-1)/otherwise(0),
    houselight = zeros(100000,2);
    
    index = 0;
    
    % read the trajectory and the detection probability
    fid = fopen(schedules);
    tline = fgetl(fid);
    
    % important information starts after line that starts with Evnt_Time
    good_stuff = 0;
    
    while ischar(tline)
        
        tline = fgetl(fid);
        
        % make sure that we are in the relevant part of the file
        if length(tline)>length('Evnt_Time')
            % if line is the header, skip it and start reading after that
            if contains(tline,'Evnt_Time')
                good_stuff = 1;
                tline = fgetl(fid);
            else
                % if line starts with a 0, we are in the data that needs to
                % be read
                if str2num(tline(1))==0
                    good_stuff = 1;
                end
            end
        end
        
        if not(isequal(tline,-1)) && good_stuff
            
            % analyse the content of each line. "," is the separator
            separator_pos = find(tline==',');
            
            % extract 4th field
            field4 = tline(separator_pos(3)+1:separator_pos(4)-1);
            
            % extract third field
            field3 = tline(separator_pos(2)+1:separator_pos(3)-1);

            % extract first field
            field1 = tline(1:separator_pos(1)-1);

            % check if this is a houselight event
            index = index + 1;
            
            % store the time
            houselight(index,1) = str2num(field1);
            
            % store flags if this is a houselight event
            if strcmp(field4,'HouseLight #1')
                % check if this is a on/off
                if strcmp(field3,'Output On Event')
                    houselight(index,2) = 1;
                else
                    if strcmp(field3,'Output Off Event')
                        houselight(index,2) = -1;
                    else
                        disp('Strange Houselight field in schedules')
                        pause
                    end
                end
            end
            
        end
    end
    
    fclose(fid);
    
    % remove the excess
    houselight(index+1:end,:) = [];
    
    
    % =========================================================================
    
    figure('Position',[484 793 1724 550]);
    
    % make summary figure
    subplot(2,1,1)
    plot(mean_bh)
    axis tight
    xlabel('Time (frames)')
    ylabel('Mean bh video signal')
    
    title(datadir,'Interpreter','None')
    
    subplot(2,1,2)
    % plot the part of schedules that was recorded by the bh camera
    Nf = length(mean_bh);
    upper = Nf*dt;
    houselight(find(houselight(:,1)>upper),:) = [];
    plot(houselight(:,1),houselight(:,2))
    axis tight
    xlabel('Time (s)')
    ylabel('houselight from schedules')
    
    % save image
    img = getframe(gcf);
    imwrite(img.cdata, [datadir '/qualitycheck/houselight_verification.png']);

    close all
    
end















