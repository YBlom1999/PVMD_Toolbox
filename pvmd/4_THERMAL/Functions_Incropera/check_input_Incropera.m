function Valid_input_Incropera = check_input_Incropera(thermal_input)
%Valid_input_Incropera checks if the input for Incropera is valid
Valid_input_Incropera = true;
Nlayers = thermal_input.Nlayers;
AssignmentLayers = thermal_input.assignment_layers;
if Nlayers ~= size(thermal_input.layers,2) || max(AssignmentLayers) > Nlayers
    Valid_input_Incropera = false;
end
end