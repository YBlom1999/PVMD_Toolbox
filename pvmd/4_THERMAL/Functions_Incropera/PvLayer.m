classdef PvLayer < handle
    % PvLayer()
    % Description: 
    % This class is used to define the properties of a layer in the PV
    % module. 
    
    properties
        name %char array
        thickness %layer thickness in m.
        k %material thermal conductivity in W/m/K.
        rho %material density in kg/m3.
        cp %material specific heat capacity in J/kg/K.
        dx %horizontal spacing between nodes in the layer in m.
        dy %vertical spacing between nodes in the layer in m.
        emissivity % material emissivity. 
    end
    properties (Dependent)
      nIntNodeY % number of vertical nodes in the layer.
   end
    
    methods
        function obj = PvLayer(varargin)
            % Constructor documentation.
            % obj = PvLayer()
            % obj = PvLayer(key1,val1,key2,val2,...)
            % Valid optional keys and values:
            % key: 'name' - value: char array
            % key: 'thickness' - value: double (m)
            % key: 'thermal conductivity' - value: double (W/m/K)
            % key: 'density' - value: double (kg/m3)
            % key: 'heat capacity' - value: double (J/kg/K)
            % key: 'dx' - value: double (m)
            % key: 'nIntNodeY' - value: double (number of vertical nodes)
            % key: 'emissivity' - value: double
            %
            % Example of use:
            % To define a typical c-Si monofacial PV module:
            % pvModLayers = [PvLayer('name','Top glass','thickness',3/1000,'k',1.8,'rho',2700,'cp',750,'dx',0.01,'nIntNodeY',6,'emissivity',0.84)
            %            PvLayer('name','Top EVA','thickness',0.5/1000,'k',0.32,'rho',960,'cp',2090,'dx',0.01,'nIntNodeY',1)
            %            PvLayer('name','Cell','thickness',0.2/1000,'k',149,'rho',2330,'cp',838,'dx',0.01,'nIntNodeY',3)
            %            PvLayer('name','Bottom EVA','thickness',0.5/1000,'k',0.32,'rho',960,'cp',2090,'dx',0.01,'nIntNodeY',1)
            %            PvLayer('name','PET','thickness',0.4/1000,'k',0.56,'rho',1370,'cp',1760,'dx',0.01,'nIntNodeY',1,'emissivity',0.89)];

            obj.name = '<Not assigned>';
            obj.thickness = NaN;
            obj.k = NaN;
            obj.rho = NaN;
            obj.cp = NaN;
            obj.dx = NaN;
            obj.dy = NaN;
            obj.emissivity = NaN;
            if ~isempty(varargin)
                numargin = length(varargin);
                if ~rem(numargin,2)%there must be a value for each key
                    intNodeAux = NaN;
                    for k=1:2:numargin
                        val = varargin{k+1};
                        switch varargin{k}
                            case {'name','Name'}
                                if ischar(val)
                                    obj.name = val;
                                else
                                    error('Name value must be char.');
                                end
                            case {'thickness','Thickness'}
                                if isnumeric(val)
                                    obj.thickness = val;
                                else
                                    error('Thickness value must be numeric.');
                                end
                            case {'thermal conductivity','Thermal conductivity','k'}
                                if isnumeric(val)
                                    obj.k = val;
                                else
                                    error('Thermal conductivity value must be numeric.');
                                end
                            case {'density','Density','rho'}
                                if isnumeric(val)
                                    obj.rho = val;
                                else
                                    error('Density value must be numeric.');
                                end
                            case {'heat capacity','specific heat capacity',...
                                    'Heat capacity','Specific heat capacity','cp'}
                                if isnumeric(val)
                                    obj.cp = val;
                                else
                                    error('Specific heat capacity value must be numeric');
                                end
                            case {'dx','DX','delta x','Delta x'}
                                if isnumeric(val)
                                    obj.dx = val;
                                else
                                    error('Delta x value must be numeric.');
                                end
                            case {'dy','DY','delta y','Delta y'}
                                if isnumeric(val)
                                    obj.dy = val;
                                else
                                    error('Delta y value must be numeric.');
                                end
                            case {'nIntNodeY','internal nodes y','Internal nodes Y','internal nodes','Internal nodes'}
                                if isnumeric(val) && val == floor(val)
                                    intNodeAux = val;
                                else
                                    error('The number of internal nodes must be an integer.');
                                end
                            case {'emissivity','Emissivity','epsilon','Epsilon'}
                                if isnumeric(val) && val <=1
                                    obj.emissivity = val;
                                else
                                    error('The emvissivity must be a a number between 0 and 1.');
                                end
                            otherwise
                                try
                                    error([varargin{k},' is an invalid input key.']);
                                catch %if you can't concatenate then the key is not a char
                                    error('One of the keys is not a char.');
                                end
                        end
                    end
                else
                    error('Invalid input argument format.');
                end
                if ~isnan(intNodeAux)
                    obj.dy = obj.thickness/(intNodeAux+1);
                end
            end
        end
        
        
        function retVal = get.nIntNodeY(obj)
            retVal = obj.thickness / obj.dy - 1;
%             if floor(retVal) ~= retVal
%                 error(['Delta y (dy) value in ',obj.name,' layer is invalid']);
%             end
        end
    end
end

