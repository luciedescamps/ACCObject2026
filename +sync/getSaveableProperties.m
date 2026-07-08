function propNames = getSaveableProperties(v)
% GETSAVEABLEPROPERTIES - returns a cell array of an object's property
% names that are mutable and persistent. Takes either an object or a class
% name as input

if isobject(v)
    % Get metaclass from obj
    mc = metaclass(v);
else
    mc = meta.class.fromName(v);
end

% Get metaproperties
mps = mc.PropertyList;

mps = mps(~[mps.Constant] & ~[mps.Transient] & ~[mps.Dependent]);

propNames = {mps.Name};

end