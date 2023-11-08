classdef TimedProgressBar < handle
    %PROGRESSBAR Progress bar class for matlab loops which also works with
    %   parfor. PROGRESSBAR works by creating a transient file in your working
    %   directory, and then keeping track of the loop's progress within that
    %   file. This workaround is necessary because parfor workers cannot
    %   communicate with one another so there is no simple way to know which
    %   iterations have finished and which haven't.
    %   Meanwhile, refrain any other competing output to the command line,
    %   otherwise it will mess it. Just use the input strings to comunicate
    %   while de cycle is going on.
    %
    % METHODS:  TimedProgressBar(); constructs an object
    %               and initializes the progress monitor
    %               for a set of targetCount upcoming calculations.
    %           progress(); updates the progress inside your loop and
    %               displays an updated progress bar.
    %           stop(); deletes progressbar_(random_number).txt and finalizes
    %               the progress bar.
    %
    % EXAMPLE:
    %           targetCount = 100;
    %           barWidth= targetCount/2;
    %           p = TimedProgressBar( targetCount, barWidth, ...
    %                                 'Computing, please wait for ', ...
    %                                 ', already completed ', ...
    %                                 'Concluded in ' );
    %           parfor i=1:targetCount
    %              pause(rand);     % Replace with real code
    %              p.progress;      % Also percent = p.progress;
    %           end
    %           p.stop;             % Also percent = p.stop;
    %
    % To get percentage numbers from progress and stop methods call them like:
    %       percent = p.progress;
    %       percent = p.stop;
    %
    % date: 2014/05/27
    % author: Antonio Jose Cacho, ajsccacho@gmailcom
    % author: Stefan Doerr (previous version: ProgressBar)
    % author: Jeremy Scheff (previous version: parfor_progress)
    
    % Modified: Dustin Harvey for Echo
    
    properties ( SetAccess= protected, GetAccess= protected )
        fname
        fstartname
        format
        waitMsg
        finishTimeMsg
        startMsg
        barWidth
        textWidth
        rewindLength
        targetCount
        startTime
        cleanup
        enabled = true
        graphical
        graphicalObj
        graphicalCancel = false
    end
    
    methods
        
        function obj = TimedProgressBar( targetCount, barWidth, waitMsg, ...
                finishTimeMsg )
            
            obj.graphical = isdeployed;
            
            % Don't create progress bars on parallel workers
            if checkParallel() && ~isempty(getCurrentWorker)
                obj.enabled = false;
            end
                
            waitSz= length( waitMsg );
            finishSz= length( finishTimeMsg );
            if waitSz > finishSz
                finishTimeMsg= [ repmat( ' ', 1, waitSz - finishSz ), ...
                    finishTimeMsg ];
            else
                waitMsg= [ repmat( ' ', 1, finishSz - waitSz ), ...
                    waitMsg ];
            end
            obj.fname= obj.setFName();
            obj.fstartname = [obj.fname(1:(end-4)),'_start.txt'];
            obj.format= [ '%03d:%02d:%02.0f. %3.0f%% complete'];
            obj.textWidth= length(obj.waitMsg) + 24;
            obj.targetCount = targetCount;
            obj.barWidth = barWidth;
            obj.waitMsg = waitMsg;
            obj.finishTimeMsg = finishTimeMsg;
            obj.startTime = tic;
            obj.rewindLength = waitSz+obj.barWidth+28;
            obj.startMsg = [char(10),obj.waitMsg,...
                repmat(' ',1,obj.rewindLength-1-length(obj.waitMsg))];
            
            % Don't create progress bars if another progress bar is active
            enabledObjFname = obj.enabledObjectStore(true);
            if ~strcmp(enabledObjFname,obj.fname)
                obj.enabled = false;
            end
            
            f= fopen( obj.fname, 'w+');
            fclose(f);
            
            if obj.enabled && obj.graphical
                obj.graphicalObj = parfor_progressbar(targetCount, obj.waitMsg,...
                    'CreateCancelBtn',@(~,~)onCancelButton(obj));
                obj.graphicalObj.wbh.WindowStyle = 'modal';
            end
        end
        
        function start(obj)
            if ~obj.enabled
                return
            end
            
            disp(obj.startMsg);
        end
        
        function progress(obj)
            if ~obj.enabled
                return
            end
            
            if obj.graphical
                if obj.graphicalCancel
                    error('Execution cancelled.');
                end
                
                obj.graphicalObj.iterate(1);
                return
            end
            
            obj.updateFile();
            
            if java.io.File(obj.fstartname).isFile()
                timeElapsed = toc(obj.startTime);
                [ percent, remainingTime,]= obj.getProgress(timeElapsed);
                timedMsg= obj.getTimedMsg( remainingTime, percent );
                obj.showStatus( percent, [ obj.waitMsg timedMsg ] );
                return
            end
            
            if checkParallel()
                primaryWorker = getCurrentTask();
            else
                primaryWorker = [];
            end
            primaryWorker = isempty(primaryWorker) || primaryWorker.ID==1;
            
            if primaryWorker && (toc(obj.startTime)>3)
                obj.start();
                                
                f = fopen( obj.fstartname, 'w+');
                fclose( f );
            end
        end
        
        function stop(obj)
            if ~obj.enabled
                return
            end
            
            if obj.graphical
                obj.graphicalObj.close();
                return
            end
            
            if java.io.File(obj.fstartname).isFile()
                timeElapsed = toc(obj.startTime);
                percent= 100;
                timedMsg= obj.getTimedMsg( timeElapsed, percent );
                obj.showStatus( percent, [ obj.finishTimeMsg timedMsg ] );
                delete(obj.fstartname);
            end
            if java.io.File(obj.fname).isFile()
                delete(obj.fname);
            end
        end
        
        function delete(obj)
            obj.enabledObjectStore(false);
        end
        
    end
    
    methods  ( Static )
        
        function fname= getFName()
            fname= [  'TPBtemp_', genUniqueString ];
            fname = [ fname(1:32), '.txt' ];
            while exist( fname, 'file' )
                fname= [  'TPBtemp_', genUniqueString ];
                fname = [ fname(1:32), '.txt' ];
            end
            fname = fullfile(tempdir,fname);
        end
        
    end
    
    
    methods  ( Access= protected )
        
        function fname= setFName(obj)
            fname= obj.getFName();
            while exist( fname, 'file' )
                fname= obj.getFName();
            end
        end
        
        function [ percent, remainingTime]= getProgress(obj, timeElapsed)
            f= fopen( obj.fname, 'r' );
            [ ~, progressCount ]= fscanf( f, '%ld' );
            fclose( f );
            
            percent = progressCount * 100 / obj.targetCount;
            if percent > 0
                remainingTime= timeElapsed * ( 100 / percent - 1 );
            else
                remainingTime= timeElapsed * 200;   % initial duration estimate
            end
        end
        
        function showStatus( obj, percent, statusMsg )
            if ~obj.enabled
                return
            end
            
            x= round( percent * obj.barWidth / 100 );
            marker= '>';
            if x < obj.barWidth
                bar= [ ' [', repmat( '=', 1, x ), marker, ...
                    repmat( ' ', 1, obj.barWidth - x - 1 ), ']' ];
            else
                bar= [ ' [', repmat( '=', 1, x ), ...
                    repmat( ' ', 1, obj.barWidth - x ), ']' ];
            end
            statusLine= [ statusMsg, bar ];
            
            % console print:
            % <--------status message----------->_[<------progress bar----->]
            % Wait for 001:28:36.7, completed 38% [========>                ]
            
            cursorRewinder= repmat( char(8), 1, obj.rewindLength );
            statusLine = [repmat( ' ', 1, obj.rewindLength-length(statusLine)-1),...
                statusLine];
            disp( [ cursorRewinder, statusLine ] );
        end
        
        function timedMsg= getTimedMsg( obj, timeSec, percent )
            hh= fix( timeSec / 3600 );
            mm= fix( ( timeSec - hh * 3600 ) / 60 );
            ss= max( ( timeSec - hh * 3600 - mm * 60 ), 0.1 );
            timedMsg = sprintf( obj.format, hh, mm, ss, percent );
        end
        
        function updateFile(obj)
            if ~obj.enabled
                return
            end
            
            f = fopen(obj.fname, 'a');
            fprintf(f, '1\n');
            fclose(f);
        end
        
        function enabledObjFname = enabledObjectStore(obj,enable)
            persistent storedObjectFname
            
            if isempty(storedObjectFname) && enable
                storedObjectFname = obj.fname;
            elseif ~isempty(storedObjectFname) && strcmp(storedObjectFname,obj.fname) && ~enable
                storedObjectFname = [];
            end
            
            enabledObjFname = storedObjectFname;
        end
        
        function onCancelButton(obj)
            obj.graphicalCancel = true;
            obj.graphicalObj.close();
        end
        
    end
    
end


function parallelAvailable = checkParallel()
%CHECKPARALLEL Check if parallel is available
% [PARALLELAVAILABLE]=CHECKPARALLEL()
%
% <strong>Output</strong>
%
% parallelAvailable ([1x1] logical)
%                   - true if Parallel Computing Toolbox installed and
%                     licensed

% Check for Parallel Computing Toolbox
persistent parallelAvailable_

if isempty(parallelAvailable_)
    verInfo = ver;
    parallelAvailable_ = ismember('Parallel Computing Toolbox',{verInfo.Name}) && ...
        license('test','Distrib_Computing_Toolbox');
end

parallelAvailable = parallelAvailable_;
end

function s = genUniqueString
%GENUNIQUESTRING Generate a unique alphanumeric string
% [S]=GENUNIQUESTRING()
%
% <strong>Output</strong>
%
% s (String)
s=[char(java.util.UUID.randomUUID()) char(java.util.UUID.randomUUID())];
s=strrep(s,'-','0');
s=['v' s];
N=namelengthmax;
s=s(1:(min(N,end)));
end