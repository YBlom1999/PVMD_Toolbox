classdef ThermalNode < handle
    
    properties
        position
        type %according to Juan's nummbering
        topLayer
        botLayer
        
    end
    
    methods
        function obj = ThermalNode()
            %UNTITLED10 Construct an instance of this class
            %   Detailed explanation goes here
            obj.position = [NaN NaN];
            obj.topLayer = PvLayer();
            obj.botLayer = PvLayer();
            obj.type = NaN;
        end
    end
end

