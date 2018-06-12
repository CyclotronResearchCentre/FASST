function out = isfieldRecursive(varargin)
%allFieldsExist = isfieldRecursive(myStructure,'fieldOfTopLevel','fieldOfFirstField','FieldOfThatNextField')
%This function has two main purposes:
%1) To check whether a field of a field of field (etc.) of a structure exists. Can be used to check as many levels deep as
%desired. Provide the first input as a structure and any number of subsequent inputs as strings
%containing field names to be checked (or cell arrays of these, see below).
%Example
% myStructure.right.directory = 'd';
%allFieldsExist = isfieldRecursive(specimen,'right','directory')
% allFieldsExist =
%      1
%2) To test multiple branches of multiple fields of structures
%recursivley until all are found to be present. To check whether more than one field exists at a
%certain level, pass a cell array of strings and fields will be searched for at that level.
%Example:
% myStructure.calibration.left.fc = 1;
% myStructure.calibration.right.fc = 1;
% myStructure.calibration.centre.fc = 1;
%allFieldsExist = isfieldRecursive(myStructure,'calibration',{'left','right','centre'})
% allFieldsExist =
%      1
%allFieldsExist = isfieldRecursive(myStructure,'calibration',{'left','right','centre'},'fc')
% allFieldsExist =
%      1
%allFieldsExist = isfieldRecursive(myStructure,'calibration',{'left','right','centre','blah'},'fc')
% allFieldsExist =
%      0
%allFieldsExist = isfieldRecursive(myStructure,'calibration',{'left','right','centre'},{'fc','kc'})
% allFieldsExist =
%      0
% myStructure.calibration.left.kc = 1;
%allFieldsExist = isfieldRecursive(myStructure,'calibration',{'left','right','centre'},{'fc','kc'})
% allFieldsExist =
%      0
% myStructure.calibration.right.kc = 1;
% myStructure.calibration.centre.kc = 1;
%allFieldsExist = isfieldRecursive(myStructure,'calibration',{'left','right','centre'},{'fc','kc'})
% allFieldsExist =
%      1
%Author: M Arthington
%Date: 2010/07/06


out = true;
while length(varargin)>=2 && out
	if iscell(varargin{2})
		for i=1:length(varargin{2})
			if length(varargin)>2
				out = out && isfieldRecursive(varargin{1},varargin{2}{i},varargin{3:end});
			else
				out = out && isfieldRecursive(varargin{1},varargin{2}{i});
			end
		end
		if out
			varargin{1} = varargin{1}.(varargin{2}{1});%Select the first cell so that the while loop can exit, even though this will already have been checked
			varargin = {varargin{1} varargin{3:end}};
		end
	else
		if out && isfield(varargin{1},varargin{2})
			varargin{1} = varargin{1}.(varargin{2});
			varargin = {varargin{1} varargin{3:end}};
		else
			out = false;
		end
	end
end