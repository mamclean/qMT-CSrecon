% ---------------------------------------------------------------------
% RAW-DATA FILE READER - qMT-offset scans saved as p-files
% Ethan Macdonald, edited by Melany Mclean
% Edited - Jun 13, 2016
%
% Any scan that has been saved as a p-file can be opened and reconstructed.
% Doesn't work because Reconstruct file is missing.
% ---------------------------------------------------------------------

classdef PfileReader
    % QSM - Class object that takes meSPGR data and uses it to create
    % susceptability maps
    %
    % 
    
%     INITIALIZE STUFF

	properties (SetAccess = public)
        
        fileName
        pathObj
        
        reconstructionMode = 1 ; 
            % 1 - 3D-GRE, SPGR
            % 2 - 
      
        PfileParams
        
    end
    
    properties (SetAccess = private)
        
    end
    
    methods (Static)
        
    end
    
    methods (Access = public)
        
        
        % FUNCTIONS - READ OR RECONSTRUCT
        
        function obj = PfileReader( path_Obj, FileName )
            
            obj.fileName = FileName ;
            obj.pathObj = path_Obj ;
            
        end
        
        function obj = Reconstruct( obj, doDisplay )
            
            if( obj.reconstructionMode == 1 )
                
                obj.PfileParams = fGetPfileParams( ... % This is in the private folder
                    [obj.pathObj.root_data_pathway ...
                    obj.pathObj.subject_directories{1}], ...
                    obj.fileName ) ;
                
            end
            
            if doDisplay == 1
                % Displays dicom information
%                 Display = zeros(obj.B0field.sourceData.xLengthFGRE, ...
%                                 obj.B0field.sourceData.yLengthFGRE, 1 , ...
%                                 obj.B0field.sourceData.zLengthFGRE ) ;
% 
%                 for Slice = 1:obj.B0field.sourceData.zLengthFGRE
%                     Display(:,:,1,Slice) = obj.volB0background(:,:,Slice);
%                 end
% 
%                 figure;
%                 montage(Display, 'DisplayRange', []);
%                 title(['Background B_0 Volume']);
%                 colorbar;
%                 drawnow; snapnow;
                
            end
            
        end
        
        
    end
    
    methods (Access = private)
        
    end
    
end

