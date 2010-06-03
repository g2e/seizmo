function [mod]=taupvmod(file)

% split file into its parts
[path,name,type]=fileparts(file);

% Import the TauP package
import edu.sc.seis.TauP.*

% Create a VelocityModel object
mod=VelocityModel();

% Establish that the model is of type modelFileType
mod.setFileType(type);

% Read in the file
mod.readVelocityFile(file);

% Set the model name
mod.setModelName(name);

% Test the model
if ~mod.validate()
    clear mod
end

end
